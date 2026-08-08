import sys
import os
import argparse
from .RunCmdsMP import logger
from .__version__ import version
import logging
from collections import OrderedDict

# 禁用 fontTools 的所有 INFO 级别日志
logging.getLogger('fontTools').setLevel(logging.WARNING)


bindir = os.path.dirname(os.path.realpath(__file__))

def args_common(parser):
	pass


def args_dotplot(parser):
	from .dot_plotter import dotplot_args
	args = dotplot_args(parser)
	return dict(args=args)


def func_dotplot(**kargs):
	from .dot_plotter import xmain as dot_plotter
	dot_plotter(**kargs)

def args_depth(parser):
	from .ploidy_plotter import ploidy_args
	args = ploidy_args(parser)
	return dict(args=args)

def func_depth(**kargs):
	from .ploidy_plotter import xmain as ploidy_plotter
	ploidy_plotter(**kargs)

def args_retention(parser):
	from .retention_plotter import retention_args
	retention_args(parser)

def func_retention(**kargs):
	from .retention_plotter import xmain
	xmain(**kargs)

def args_filter(parser):
	parser.add_argument('-s', '-synteny', required=True,  type=str,  nargs='*',
						dest='collinearities',  metavar='FILE',
						help="Collinearity files output from MCscanX, WGDI, or MCscan/JCVI. \
[required]")
	parser.add_argument('-o', '-orthology', required=True,  type=str,  nargs='*',
						dest='orthologs',  metavar='FOLDER/FILE',
						help="Orthologues output from OrthoFinder (folder), or OrthoMCL (file). \
[required]")
	parser.add_argument('-c', '-cutoff',  type=float, default=0.6,
						dest='min_ratio',  metavar='FLOAT',
						help="Cutoff (lower limit) of Orthology Index (rataining blocks > this cutoff) [default=%(default)s]")
	parser.add_argument('-u', '-upper',  type=float, default=1,
						dest='max_ratio',  metavar='FLOAT',
						help="Upper limit of Orthology Index (rataining blocks <= this cutoff) [default=%(default)s]")
	parser.add_argument('-n', '-min_n',  type=int, default=0,
						dest='min_n',  metavar='INT',
						help="Minimum gene number in a block [default=%(default)s]")
	parser.add_argument('-g', '-gff', type=str, nargs='+',
						dest='gff',  metavar='FILE',
						help="Gff file. [required for `-d`]")
	parser.add_argument('-d', '-min_dist',  type=int, default=None,
						dest='min_dist',  metavar='INT',
						help="Minimum distance to remove a tandem repeated block [default=%(default)s]")
	parser.add_argument('-stat', default=None, metavar='FILE',
						dest='out_stats', type=str,
						help="Output stats by species pairs. [default=%(default)s]")
	parser.add_argument('-oo', default=False,
						dest='output_orthology', action='store_true',
						help="Output retained orthology instead of synteny. [default=%(default)s]")


def func_filter(**kargs):
	from .mcscan import identify_orthologous_blocks
	identify_orthologous_blocks(**kargs)

def _add_shared_hog_args(parser, og_required=True, s_required=True,
						 t_required=True, group=None):
	"""Shared arguments for prune and hog subcommands.

	If group (an argparse argument group) is given, arguments are added to
	it instead of the parser directly.
	"""
	target = group if group is not None else parser
	target.add_argument('-og', '-orthogroup', required=og_required, type=str,
						dest='ogfile', metavar='FILE',
						help='Orthogroup file (MCL format)')
	target.add_argument('-s', '-synteny', required=s_required, type=str, nargs='+',
						dest='orthfiles', metavar='FILE',
						help='Ortholog/Collinearity files')
	target.add_argument('-t', '-sptree', required=t_required, type=str,
						dest='sptreefile', metavar='FILE',
						help='Species tree file (Newick)')
	target.add_argument('-inparalog', action='store_true', default=False,
						dest='paralog',
						help='Include paralogs in `-s` input. Do NOT '
							 'enable if the input paralogs is not terminal inparalogs')
	target.add_argument('--min-child-species', type=int, default=2,
						dest='min_child_species', metavar='INT',
						help=argparse.SUPPRESS)  # default 2, hidden
	target.add_argument('--cross-speciation', action='store_true', default=True,
						dest='cross_speciation',
						help=argparse.SUPPRESS)  # default True, hidden
	target.add_argument('--drop-no-cross', action='store_true', default=False,
						dest='drop_no_cross',
						help='Drop (instead of merge) HOGs whose genes do not span all child branches')


def args_hog(parser):
	g_build = parser.add_argument_group('HOG building')
	_add_shared_hog_args(parser, group=g_build)
	g_out = parser.add_argument_group('Output')
	g_out.add_argument('-pre', '-prefix', type=str, default='HOGs',
					   dest='outpre', metavar='FILE',
					   help='Output prefix [default=%(default)s]')
	g_out.add_argument('--out-stats', action='store_true', default=False,
					   help='Output copy-number statistics TSV (<prefix>.stats.tsv)')
	g_out.add_argument('--bar-plot', action='store_true', default=False,
					   help='Output bar chart of copy-number distribution (<prefix>.bar.pdf/.png)')
	g_out.add_argument('--tree-plot', action='store_true', default=False,
					   help='Output species tree with copy-number pie charts at nodes (<prefix>.tree.pdf/.png)')
	g_out.add_argument('--max-copies', type=int, default=5,
					   help='Max copy number to track in stats/plot [default=%(default)s]')

def func_hog(**kargs):
	from .hog import xmain as hog_main
	hog_main(**kargs)
def args_paralog(parser):
	g_hog = parser.add_argument_group('HOG rebuilding')
	_add_shared_hog_args(parser, og_required=False, s_required=False,
						 t_required=False, group=g_hog)
	g_hog.add_argument('--write-hog', action='store_true', default=True,
					   dest='write_hog',
					   help=argparse.SUPPRESS)  # default True, hidden

	g_in = parser.add_argument_group('Prebuilt inputs (skip rebuilding)')
	g_in.add_argument('--hog', type=str, default=None,
					  dest='hog_tsv', metavar='FILE',
					  help='Load HOGs from existing HOGs.tsv (skip rebuilding); '
						   'use the same species tree as when it was built')
	g_in.add_argument('--paralog', type=str, default=None,
					  dest='inparalog', metavar='FILE',
					  help='Load paralog pairs from existing paralog.tsv '
						   '(first 3 columns only)')

	g_ix = parser.add_argument_group('Paralog indexing')
	g_ix.add_argument('-ss', '--self-synteny', type=str, nargs='+', default=None,
					  dest='self_synteny', metavar='FILE',
					  help='Synteny file for paralog indexing (default: same as -s; '
						   'required if -s is not given)')
	g_ix.add_argument('-n', '-min_n', type=int, default=0,
					  dest='min_n', metavar='INT',
					  help='Minimum gene number in a block [default=%(default)s]')
	g_ix.add_argument('-g', '-gff', type=str, nargs='+',
					  dest='gff', metavar='FILE',
					  help='Gff file. [required for `-d`]')
	g_ix.add_argument('-d', '-min_dist', type=int, default=None,
					  dest='min_dist', metavar='INT',
					  help='Minimum distance to remove a tandem repeated block [default=None]')
	g_ix.add_argument('--pi-cutoff', type=float, default=0.05,
					  dest='pi_cutoff', metavar='FLOAT',
					  help='BPI cutoff: blocks whose best-branch BPI falls below '
						   'this are assigned to the root (BPI = paralog_pairs / '
						   'block_gene_pairs) [default=%(default)s]')
	g_ix.add_argument('--nodes', metavar='NODE', nargs='+', type=str, default=None,
					  dest='nodes',
					  help='Tree nodes to report paralogs for (default: all)')
	g_ix.add_argument('--species', metavar='SPECIES', nargs='+', type=str, default=None,
					  dest='species',
					  help='Species to report paralogs for (default: all)')

	g_out = parser.add_argument_group('Output')
	g_out.add_argument('-pre', '-prefix', type=str,
					   dest='prefix', metavar='PREFIX', default='paralog_index',
					   help='Output prefix [default=%(default)s]')
	g_out.add_argument('--no-index', action='store_true', default=False,
					   dest='no_index',
					   help='Disable paralog indexing, output raw paralog pairs')
	g_out.add_argument('--tree-plot', action='store_true', default=False,
					   dest='tree_plot',
					   help='Output species tree with syntenic/non-syntenic paralog pies '
							'(<prefix>.tree.pdf/.png)')
	g_out.add_argument('--heatmap', action='store_true', default=False,
					   dest='heatmap',
					   help='Output block x branch PI heatmap '
							'(<prefix>.heatmap.tsv/.pdf/.png); -t gives '
							'tree-ordered columns')
	g_out.add_argument('--cluster', action='store_true', default=False,
					   dest='heatmap_cluster',
					   help=argparse.SUPPRESS)
	g_out.add_argument('--scale', type=str, default=None,
					   choices=['linear', 'log', 'log2', 'log10', 'sqrt'],
					   dest='heatmap_scale',
					   help='Scale heatmap row height by block gene count '
							'(default: uniform rows) [choices: %(choices)s]')

def func_paralog(**kargs):
	from .paralog import xmain as paralog_main
	paralog_main(**kargs)


def args_detandem(parser):
	parser.add_argument('-og', '-orthogroup', required=True, type=str,
						dest='ogfile', metavar='FILE',
						help='Orthogroup file (MCL format). [required]')
	parser.add_argument('-g', '-gff', required=True, type=str, nargs='+',
						dest='gfffile', metavar='FILE',
						help='GFF annotation file. [required]')
	parser.add_argument('-s', '-synteny', type=str, nargs='*', default=None,
						dest='orthfiles', metavar='FILE',
						help='Ortholog/Collinearity files for selecting the '
							 'gene with highest ortholog degree per cluster. [optional]')
	parser.add_argument('-d', '-dist', type=int, default=200,
						dest='tandem_dist', metavar='INT',
						help='Maximum index distance for tandem duplication. '
						'[default=%(default)s]')
	parser.add_argument('--out-tandem', type=str, default=None,
						dest='out_tandem', metavar='FILE',
						help='Output tandem clusters (retained gene marked with *)')

def func_detandem(**kargs):
	from .detandem import Detandem
	Detandem(**kargs).run()

def args_evaluate(parser):
	from .evaluator import evaluate_args
	evaluate_args(parser)


def args_ksplot(parser):
	from .ks_plotter import ksplot_args
	ksplot_args(parser)

def args_prune(parser):
	_add_shared_hog_args(parser)
	parser.add_argument('-o', '--output', type=str,
						dest='output', metavar='FILE', default=None,
						help='Output file path [default: stdout]')
	parser.add_argument('-restore-gene','--restore-gene', action='store_true', 
						dest='restore_gene',default=False,
						help='When a species has no genes retained, restore the first gene of that species (default: False)')
	parser.add_argument('-restore-log','--restore-log', type=str,
						dest='restore_log',
						help='Log file path for recording restored genes')


def func_prune(**kargs):
	import sys
	from .cluster_filter_from_hog import process_og_with_hog
	hog_args = {
		'ogfile': kargs['ogfile'],
		'orthfiles': kargs['orthfiles'],
		'sptreefile': kargs['sptreefile'],
		'paralog': kargs['paralog']
	}
	output = kargs.get('output') or sys.stdout
	process_og_with_hog(kargs['ogfile'], hog_args, output, kargs['restore_gene'], kargs['restore_log'])

def func_evaluate(**kargs):
	from .evaluator import evaluate_main
	evaluate_main(**kargs)


def func_ksplot(**kargs):
	from .ks_plotter import xmain as ksplot_main
	ksplot_main(**kargs)

def args_sim(parser):
	try:
		from .evolution_simulator_ak import sim_args
	except ImportError as e:
		logger.warning('Cannot register sim: {}'.format(e))
		return
	sim_args(parser)

def func_sim(**kargs):
	from .evolution_simulator_ak import xmain as sim_main
	sim_main(**kargs)

def args_rak(parser):
	parser.add_argument('-og', '-orthogroup', required=True, type=str,
						dest='ogfile', metavar='FILE',
						help='Orthogroup file [required]')
	parser.add_argument('-s', '-synteny', required=True, type=str,
						dest='orthfiles', metavar='FILE', nargs='*',
						help='Ortholog/Collinearity files [required]')
	parser.add_argument('-t', '-sptree', required=True, type=str,
						dest='sptreefile', metavar='FILE',
						help='Species tree file [required]')
	parser.add_argument('-g', '-gff', required=True, type=str, nargs='+',
						dest='gfffile', metavar='FILE',
						help='GFF annotation file [required]')
	parser.add_argument('-prefix', type=str, default='AKR',
						dest='outpre', metavar='FILE',
						help='Output prefix [default=%(default)s]')
	parser.add_argument('-paralog', default=False,
						dest='paralog', action='store_true',
						help='Include paralogs [default=%(default)s]')
	parser.add_argument('-rounds', type=int, default=3,
						dest='rounds', metavar='INT',
						help='Optimization rounds [default=%(default)s]')
	parser.add_argument('-chrlist', type=str, default=None,
						dest='chrom_list', metavar='FILE',
						help='Chromosome list file (one per line) to retain [default=%(default)s]')
	parser.add_argument('-mingenes', type=int, default=0,
						dest='min_genes', metavar='INT',
						help='Minimum gene number per chromosome to retain [default=%(default)s]')
	parser.add_argument('--gfa-debug', action='store_true', default=False,
						dest='gfa_debug',
						help='Output GFA files at each pipeline phase for debugging')

def func_rak(**kargs):
	from .AK import AKR
	akr = AKR(**kargs)
	akr.run()

def args_rakeval(parser):
	from .eval_ak import add_eval_args
	add_eval_args(parser)

def func_rakeval(**kargs):
	from .eval_ak import eval_main
	eval_main(**kargs)

def args_cluster(parser):
	parser.add_argument('-s', '-synteny', required=True,  type=str,  nargs='*',
						dest='collinearities',  metavar='FILE',
						help="Collinearity files from `filter` sub-command. [required]")
	parser.add_argument('-o', '-orthology', type=str,  nargs='*',
						dest='orthologs',  metavar='FOLDER/FILE', default=None,
						help="Orthologues output from OrthoFinder (folder), or OrthoMCL (file). \
This will use Orthology Index as weight for MCL [default=%(default)s]")
	parser.add_argument('-I', '-inflation', type=float, default=1.5, metavar='FLOAT',
						dest='inflation',
						help="Inflation for MCL (varying this parameter affects granularity) \
[default=%(default)s]")
	parser.add_argument('-outgroup', type=str, default=None, metavar='TAXON/FILE',
						dest='outgroup', nargs='*',
						help="Outgroups to exclude from orthogroups (prior to `-ingroup`) \
[default=%(default)s]")
	parser.add_argument('-ingroup', type=str, default=None, metavar='TAXON/FILE',
						dest='ingroup', nargs='*',
						help="Ingroups that are only included [default=%(default)s]")
	parser.add_argument('-prefix', type=str, default='cluster',
						dest='outpre',
						help="Output prefix [default=%(default)s]")
	parser.add_argument('-m','-method', type=str, default='mcl',
						dest='method',
						help="cluster method (mcl, comp) [default=%(default)s]")

def func_cluster(**kargs):
	from .mcscan import cluster_by_mcl
	cluster_by_mcl(**kargs)


def args_outgroup(parser):
	parser.add_argument('-s', '-synteny', required=True,  type=str,  nargs='*',
						dest='collinearities',  metavar='FILE',
						help="Collinearity files from `filter` sub-command. [required]")
	parser.add_argument('-og', '-orthogroup', required=True,  type=str,
						dest='orthogroup',  metavar='FILE',
						help="Orthogroups output from `cluster` sub-command. [required]")
	parser.add_argument('-outgroup', type=str, required=True, metavar='TAXON',
						dest='outgroup', nargs='*',
						help="Outgroups to include to orthogroups [required]")
	parser.add_argument('-cutoff',  type=float, default=0.2,
						dest='min_ratio',  metavar='FLOAT',
						help="Cutoff (lower limit) of links to outgroup genes [default=%(default)s]")


def func_outgroup(**kargs):
	from .mcscan import cluster_add_outgroup
	cluster_add_outgroup(**kargs)


def args_phylo_common(parser):
	parser.add_argument('-mc', '-max_copies', type=float, default=6,
						dest='max_copies', metavar='INT',
						help="To limit a common maximum copy number for every species. [default=%(default)s]")
	parser.add_argument('-sc', '-singlecopy', default=None,
						dest='singlecopy', action='store_true',
						help="Only retrieve singlecopy genes (=`-max_copies 1`). [default=%(default)s]")
	parser.add_argument('-ss', '-spsd', type=str, default=None,
						dest='spsd',  metavar='FILE',
						help="To limit a specific copy number for each species (format: 'TAXON<tab>NUMBER'). [default=%(default)s]")
	parser.add_argument('-fmt', type=str, default='orthomcl',
						dest='source', choices=['orthomcl', 'orthofinder', 'mcscanx'],
						help="Format of `-orthogroup` input. [default=%(default)s]")
	parser.add_argument('-tree_tool', type=str, default='iqtree',
						dest='tree_tool',  metavar='STR',
						choices=['iqtree', 'fasttree'], 
						help="Tree building tool: iqtree or fasttree. [default=%(default)s]")

def args_stats(parser):
	parser.add_argument('-og', '-orthogroup', required=True,  type=str,
						dest='input',  metavar='FILE',
						help="Orthogroups output from `cluster` or `outgroup` sub-commands. [required]")
	parser.add_argument('-mm', '-max_missing', type=float, default=0.4,
						dest='max_taxa_missing', metavar='FLOAT',
						help="To allow maximum ratio of missing species. [default=%(default)s]")
	args_phylo_common(parser)


def func_stats(**kargs):
	from .mcscan import orthomcl_stats
	orthomcl_stats(**kargs)


def args_phylo(parser):
	import uuid
	uid = uuid.uuid1()
	default_tmpdir = './tmp-{}'.format(uid)

	parser.add_argument('-og', '-orthogroup', required=True,  type=str,
						dest='input',  metavar='FILE',
						help="Orthogroups output from `cluster` or `outgroup` sub-commands. \
[required]")
	parser.add_argument('-pep', required=True,  type=str,
						dest='pep',  metavar='PEP FILE',
						help="Protein fasta file. [required]")
	parser.add_argument('-cds', type=str, default=None,
						dest='cds',  metavar='CDS FILE',
						help="CDS fasta file. [default=%(default)s]")
	parser.add_argument('-both', default=False,
						dest='both', action='store_true',
						help="To use both CDS and PEP to build gene trees (only valid when `-cds` is true). \
[default: %(default)s]")
	parser.add_argument('-root', '-outgroup', type=str, metavar='TAXON',
						dest='root', nargs='*', default=None,
						help="Outgroups to root gene trees [default=%(default)s]")
	parser.add_argument('-pre', '-prefix', type=str, default='sog',
						dest='suffix', metavar='STR',
						help="Output prefix. [default=%(default)s]")
	parser.add_argument('-mm', '-max_missing', type=float, default=0.4,
						dest='max_taxa_missing', metavar='FLOAT',
						help="To allow maximum ratio of missing species. [default=%(default)s]")
	args_phylo_common(parser)

	parser.add_argument('-aligner', type=str, default='muscle',
						dest='aligner', metavar='STR', choices=['muscle', 'mafft'], 
						help="Aligner: muscle (v5/v3) or mafft. [default=%(default)s]")

	parser.add_argument('-only_aln', default=False,
						dest='onlyaln', action='store_true',
						help="Only aligning sequences, to skip trimal and iqtree. [default: %(default)s]")
	parser.add_argument('-concat', default=False,
						dest='concat', action='store_true',
						help="To concatenate alignments for tools such as IQTREE \
(valid when `-singlecopy` is true). [default: %(default)s]")
	parser.add_argument('-trimal_opts', type=str, default='-automated1',
						dest='trimal_opts',  metavar='STR',
						help="TrimAl options. [default='%(default)s']")
	parser.add_argument('-iqtree_opts', type=str, default='',
						dest='iqtree_opts',  metavar='STR',
						help="IQ-TREE options. [default='%(default)s']")
	parser.add_argument('-fast', default=False,
						dest='fast', action='store_true',
						help="Speedup IQ-TREE by restricting model set (JTT for PEP and GTR for CDS). \
[default: %(default)s]")
	parser.add_argument('-p', '-ncpu', type=int, default=20,
						dest='ncpu', metavar='INT',
						help="Number of processors. [default=%(default)s]")
	parser.add_argument('-tmp', '-tmpdir', type=str, default=default_tmpdir,
						dest='tmpdir',  metavar='FOLDER',
						help="Temporary folder. [default=%(default)s]")
	parser.add_argument('-clean', default=False,
						dest='clean', action='store_true',
						help="Cleanup temporary folder. [default: %(default)s]")
	# overwrite


def func_phylo(**kargs):
	from .mcscan import orthomcl_to_astral
	orthomcl_to_astral(**kargs)


# ── Subcommand grouping for help display ──────────────────────────────

class GroupedHelpFormatter(argparse.RawDescriptionHelpFormatter):
	"""argparse HelpFormatter that groups subcommands under categorized headers."""

	def _format_action(self, action):
		# Intercept subparser actions that carry group metadata
		if hasattr(action, '_cmd_groups'):
			return self._format_grouped_subparsers(action)
		return super()._format_action(action)

	def _format_grouped_subparsers(self, action):
		# --- header (same logic as _format_action) ---
		help_position = min(self._action_max_length + 2,
							self._max_help_position)
		help_width = max(self._width - help_position, 11)
		action_width = help_position - self._current_indent - 2
		action_header = self._format_action_invocation(action)

		if not action.help:
			tup = self._current_indent, '', action_header
			action_header = '%*s%s\n' % tup
			indent_first = 0
		elif len(action_header) <= action_width:
			tup = self._current_indent, '', action_width, action_header
			action_header = '%*s%-*s  ' % tup
			indent_first = 0
		else:
			tup = self._current_indent, '', action_header
			action_header = '%*s%s\n' % tup
			indent_first = help_position

		parts = [action_header]

		if action.help and action.help.strip():
			help_text = self._expand_help(action)
			if help_text:
				help_lines = self._split_lines(help_text, help_width)
				parts.append('%*s%s\n' % (indent_first, '', help_lines[0]))
				for line in help_lines[1:]:
					parts.append('%*s%s\n' % (help_position, '', line))
		elif not action_header.endswith('\n'):
			parts.append('\n')

		# --- grouped command listing ---
		flat = [(n, h) for names in action._cmd_groups.values() for n, h in names]
		w = max(len(n) for n, _ in flat)
		for group_name, cmds in action._cmd_groups.items():
			parts.append('  ' + group_name + ':\n')
			for name, help_text in cmds:
				parts.append('    {:<{}}  {}\n'.format(name, w, help_text))
			parts.append('\n')

		return self._join_parts(parts)


CMD_GROUPS = OrderedDict([
	('Visualization', [
		('dotplot', 'Generate Ks/OI/subgenome/ancestor-colored dot plots with versatile functions.'),
		('depth',   'Generate mutiple bar charts for synteny depth (indicator of relative ploidy).'),
		('retention', 'Plot gene retention rate (syntenic/total) along chromosomes by sliding window.'),
		('ksplot',  'Plot mutiple Ks distributions: histogram, density, and ridge plots.'),
		('evaluate', 'Evaluate and compare synteny: synteny decay, fractionation rate, OI, and gene copy-number.'),
	]),
	('Syntenic Orthogroups', [
		('filter',   'Filter synteny by Orthology Index to generate orthologous synteny.'),
		('cluster',  'Cluster orthologous synteny into syntenic orthogroups (SOGs).'),
		('outgroup', 'Add outgroups to SOGs.'),
		('detandem', 'Remove tandem duplicate genes from SOGs.'),
	]),
	('Hierarchical Orthologous Groups', [
		('hog',      'Split Hierarchical Orthologous Groups (HOGs) from SOGs using synteny and species tree.'),
		('paralog',  'Output paralogous gene pairs and synteny blocks per branch based on HOGs.'),
		('prune',    'Prune SOGs to single-copy per species based on HOGs.'),
	]),
	('Phylogenomics', [
		('phylo', 'Reconstruct gene trees from SOGs.'),
		('stats', 'Make summary of SOGs for phylogeny.'),
	]),
	('Karyotype Evolution [experimental]', [
		('rak', 'Reconstruct ancestral karyotypes based on HOGs and telomere-centric model. [experimental]'),
		('sim', 'Simulate chromosome rearrangement evolution.'),
		('rakeval', 'Evaluate karyotype reconstruction against simulation truth.'),
	]),
])

# args_* function lookup for grouped subparser creation
_ARGS_FN = {
	'dotplot': args_dotplot, 'depth': args_depth, 'retention': args_retention,
	'ksplot': args_ksplot,
	'evaluate': args_evaluate,
	'filter': args_filter, 'cluster': args_cluster, 'outgroup': args_outgroup,
	'hog': args_hog, 'detandem': args_detandem,
	'paralog': args_paralog,
	'prune': args_prune,
	'phylo': args_phylo, 'stats': args_stats,
	'rak': args_rak, 'sim': args_sim,
	'rakeval': args_rakeval,
}


def makeArgs():
	parser = argparse.ArgumentParser(
		formatter_class=GroupedHelpFormatter,
		description='Play with Orthology Index and orthologous synteny.',
	)
	parser.add_argument(
		'-v', '--version',
		action='version',
		version=version
	)
	# subcommands
	subparsers = parser.add_subparsers(help='sub-command help')
	subparsers._cmd_groups = CMD_GROUPS
	for group_name, cmds in CMD_GROUPS.items():
		for cmd_name, help_text in cmds:
			p = subparsers.add_parser(cmd_name, help=help_text)
			_ARGS_FN[cmd_name](p)

	if len(sys.argv) == 1:
		parser.print_help(sys.stderr) 
		sys.exit(1) 
	args = parser.parse_args()
	return args


FUNC = {
	'dotplot': func_dotplot,
	'filter': func_filter,
	'cluster': func_cluster,
	'outgroup': func_outgroup,
	'phylo': func_phylo,
	'stats': func_stats,
	'depth': func_depth,
	'retention': func_retention,
	'hog': func_hog,
	'detandem': func_detandem,
	'rak': func_rak,
	'sim': func_sim,
	'rakeval': func_rakeval,
	'ksplot': func_ksplot,
	'evaluate': func_evaluate,
	'paralog': func_paralog,
	'prune': func_prune,
}


def main():
	args = makeArgs()  # options
	logger.info('Command: {}'.format(' '.join(sys.argv)))
	logger.info('Version: {}'.format(version))
	logger.info('Arguments: {}'.format(args.__dict__))
	key = sys.argv[1]
	func = FUNC[key]  # functions
	func(**args.__dict__)  # execute
	logger.info('Completed\n')


if __name__ == '__main__':
	main()
