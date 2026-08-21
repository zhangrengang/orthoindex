# coding: utf-8
import sys
import os
import glob
import re
import itertools
from Bio import SeqIO
from Bio import Phylo
import numpy as np
try:
	from xopen import xopen as open
except ImportError:
	from .small_tools import open_file as open
from .small_tools import mkdirs, flattern
from collections import Counter, OrderedDict
from lazy_property import LazyWritableProperty as lazyproperty
from .RunCmdsMP import run_cmd, run_job, logger


def catAln(inALNs, outALN, allow_missing=True, idmap=None):
	'''concatenate alignments'''
	if len(inALNs) == 0:
		raise ValueError('Zero inputs')
	if allow_missing:
		species = set([])
		for inALN in inALNs:
			if not os.path.exists(inALN):
				continue
			for rc in SeqIO.parse(inALN, 'fasta'):
				if idmap is None:
					sp, g = gene_format_common(rc.id)
				else:
					try:
						sp = idmap[rc.id]
					except KeyError:
						continue
				species.add(sp)
	names = [os.path.basename(aln).split('.')[0] for aln in inALNs]
	d_seqs = {}
	lens = []
	nm = 0
	for inALN in inALNs:
		if not os.path.exists(inALN):
			nm += 1
			if nm < 10:
				logger.warn('{} not exists'.format(inALN))
			elif nm == 10:
				logger.warn('...')
			continue
		ntax = 0
		sps = set([])
		for rc in SeqIO.parse(inALN, 'fasta'):
			ntax += 1
			if idmap is None:
				sp, g = gene_format_common(rc.id)
			else:
				try:
					sp = idmap[rc.id]
				except KeyError as e:
					logger.warn(
						'{} in {} not exists in idmap'.format(rc.id, inALN))
					seq = ''
					continue
			if sp in sps:
				logger.warn('{} repeats in {}; ignored'.format(rc.id, inALN))
				continue
			seq = str(rc.seq)
			try:
				d_seqs[sp] += [seq]
			except KeyError:
				d_seqs[sp] = [seq]
			sps.add(sp)
		lens += [len(seq)]
		if allow_missing:
			for sp in species-sps:
				ntax += 1
				seq = ''.join(['-']*len(seq))
				try:
					d_seqs[sp] += [seq]
				except KeyError:
					d_seqs[sp] = [seq]
	xlens = ','.join(map(str, lens))
	names = ','.join(names)
	if nm > 0:
		logger.warn('{} alignment files missing'.format(nm))
	description = 'taxa:{} genes:{} sites:{} blocks:{} names:{}'.format(
		ntax, len(lens), sum(lens), xlens, names)
	logger.info('Concatenating {} genes of {} taxa ({} sites)'.format(
		len(lens), ntax, sum(lens)))
	for sp, seqs in list(d_seqs.items()):
		seqs = ''.join(seqs)
		print('>{} {}\n{}'.format(sp, description, seqs), file=outALN)


class Group():  # Parse groups.tsv and iterate to return each row.
	def __init__(self, inGrp):
		self.inGrp = inGrp

	def __iter__(self):
		return self.parse()

	def parse(self):
		i = 0
		for line in open(self.inGrp):
			i += 1
			if i == 1:
				self.species = line.strip().split('\t')
				continue
			yield GroupRecord(line, self.species)


class GroupRecord(object):  # Parse each line
	def __init__(self, line, species):
		line = line.strip('\n\r').split('\t')
		self.ogid = self.id = line[0]
		self.genes = [genes.split(', ') for genes in line[1:]]
		self.species = species
		self.genes = [self._strip(genes) for genes in self.genes]
		self.raw_genes = line[1:]
		self.nsp = len(self.species)

	def __iter__(self):
		return iter(self.genes)

	def get_group(self, sps=None):
		'''Get all genes for the specified set of species (default to all).'''
		if isinstance(sps, str):
			sps = {sps}
		else:  # list, set
			pass
		for sp, genes in zip(self.species, self.genes):
			if sps is not None and sp not in set(sps):
				continue
			if isinstance(genes, str):
				genes = [genes]
			for gene in genes:
				yield gene

	@property
	def counts(self):
		# self.genes]
		return [len(genes) for genes in list(self.spdict.values())]

	@property
	def spdict(self):
		return OrderedDict(list(zip(self.species, self.genes)))

	@property
	def counter(self):
		return {sp: len(genes) for sp, genes in sorted(self.spdict.items())}

	@property
	def singlecopy_ratio(self):
		try:
			return 1.0*len(self.singlecopy_dict)/self.nsp
		except ZeroDivisionError:
			return 0

	@property
	def singlecopy_dict(self):
		return {genes[0]: sp for sp, genes in list(self.spdict.items()) if len(genes) == 1}

	@property
	def singlecopy_genes(self):
		return list(self.singlecopy_dict.keys())

	def _strip(self, values):
		return [v for v in values if v]


def copy_number_density(OFdir, outPrefix='orthogroups'):
	datafile = '{}.density.data'.format(outPrefix)
	d_count = {}
	for group in OrthoFinder(OFdir).orthogroups:
		for sp, count in list(group.counter.items()):
			try:
				d_count[sp][count] += 1
			except KeyError:
				try:
					d_count[sp][count] = 1
				except KeyError:
					d_count[sp] = {count: 1}
	with open(datafile, 'w') as f:
		line = ['species', 'gene_count', 'frequency']
		print('\t'.join(line), file=f)
		for sp, counter in sorted(d_count.items()):
			for count, freq in sorted(counter.items()):
				if count < 1:
					continue
				line = [sp, count, freq]
				line = list(map(str, line))
				print('\t'.join(line), file=f)
	rsrc = outPrefix + '.density.r'
	outfig = outPrefix + '.density.pdf'
	xlim = 6
	with open(rsrc, 'w') as f:
		print('''datafile = '{datafile}'
data = read.table(datafile, head=T)
library(ggplot2)
p <- ggplot(data, aes(x=gene_count, y=frequency, color=species)) + \
	geom_line() + \
	scale_x_continuous(breaks=seq(1, {xlim}, 1), limits=c(1, {xlim}))
ggsave('{outfig}', p, width=12, height=7)
'''.format(datafile=datafile, outfig=outfig, xlim=xlim), file=f)
	cmd = 'Rscript {}'.format(rsrc)
	os.system(cmd)


class OrthoMCLGroup():  # Parse groups.txt and iterate to return each line.
	def __init__(self, inGrp, sps=None, min_genes=2):
		self.inGrp = inGrp
		self.sps = sps
		self.min_genes = min_genes

	def __iter__(self):
		return self.parse()

	def parse(self):
		for line in open(self.inGrp):
			rc = OrthoMCLGroupRecord(line, sps=self.sps)
			if len(rc.genes) < self.min_genes:
				continue
			yield rc

	def get_species(self):
		species = set([])
		for rc in self:
			species = species | set(rc.species)
		return species

def mcl2count(mclgroup, fout=sys.stdout):
	species = set([])
	counters = []
	for rc in OrthoMCLGroup(mclgroup):
		counter = rc.counter
		species = species| set(counter)
		counters += [(rc.ogid, counter)]
	species = sorted(species)
	line = ['ID'] + species
	fout.write('\t'.join(line) + '\n')
	for ogid, counter in counters:
		count = [counter.get(sp, 0) for sp in species]
		line = [ogid] + list(map(str, count))
		fout.write('\t'.join(line) + '\n')

class OrthoMCLGroupRecord(GroupRecord):  # Parse each line.
	def __init__(self, line=None, genes=None, ogid=None, nsp=None, sps=None):
		if line is not None:
			line = line.strip().split()
			self.ogid = line[0].strip(':')
			self.genes = line[1:]
		if genes is not None:
			self.genes = genes
		if ogid is not None:
			self.ogid = ogid
		self.species = self.get_species()
		if sps is not None:
			self.genes = list(self.get_group(sps=sps))
			self.species = self.get_species()
		self.nsp = nsp if nsp else len(self.species)

	def __len__(self):
		return len(self.genes)

	def get_species(self):
		return [gene.split('|')[0] for gene in self.genes]

	@property
	def spdict(self):
		d = OrderedDict()
		for sp, gene in zip(self.species, self.genes):
			try:
				d[sp] += [gene]
			except KeyError:
				d[sp] = [gene]
		return d

	@property
	def mean_copies(self):
		copies = [len(genes) for sp, genes in list(self.spdict.items())]
		return 1.0 * sum(copies) / len(copies)

	def write(self, fout):
		print('{}: {}'.format(self.ogid, ' '.join(self.genes)), file=fout)

	def classify(self):
		d_count = self.counter
		for sp, count in list(d_count.items()):
			if len(d_count) == 1:
				if count > 1:
					catlog = 2, 'Uniuqe OGs'
				elif count == 1:
					catlog = 3, 'Unassigned genes'
			else:
				if count > 1:
					catlog = 1, 'Multiple-copy OGs'
				elif count == 1:
					catlog = 0, 'Single-copy OGs'
			yield sp, catlog, count


def parse_species(species, orthofinder=None):
	if isinstance(species, str):
		if exists_and_size_gt_zero(species):
			species = [line.strip().split()[0] for line in open(species)]
		else:
			species = species.split()
	elif orthofinder is not None:
		species = orthofinder.Species
	return species


def to_hybpiper(ResultsDir, cdsSeq=None, outOGSeq=None, species=None,
				min_singlecopy=0.7, only_stats=False):
	species = parse_species(species)
	Orthogroups = OrthoFinder(ResultsDir).Orthogroups
	if not only_stats:
		d_seqs = seq2dict(cdsSeq)
	ratios = []
	i = 0
	for og in Group(Orthogroups):
		if species:
			og.genes = [og.spdict[sp] for sp in species]
			og.species = species
		if only_stats:
			ratios += [og.singlecopy_ratio]
			continue
		if not og.singlecopy_ratio >= min_singlecopy:
			continue
		i += 1
		for gene, sp in list(og.singlecopy_dict.items()):
			rc = d_seqs[gene]
			rc.id = '{}-{}'.format(sp.replace('-', '_'), og.id)
			SeqIO.write(rc, outOGSeq, 'fasta')
	if only_stats:
		print('{}\t{}'.format('total', len(ratios)), file=sys.stderr)
		for cutoff in range(50, 105, 5):
			cutoff = cutoff/1e2
			ratios = [x for x in ratios if x >= cutoff]
			print('>={}\t{}'.format(cutoff, len(ratios)), file=sys.stderr)
		return
	print('{} OGs'.format(i), file=sys.stderr)


def venn(ResultsDir, outTsv, species=None):
	species = parse_species(species)
	result = OrthoFinder(ResultsDir)
	d_groups = {}
	for group in result.get_orthogroups(sps=species):  # get_orthologs_cluster
		ogid = group.ogid
		for sp in set(group.species):
			try:
				d_groups[sp] += [ogid]
			except KeyError:
				d_groups[sp] = [ogid]
	for sp, ogs in sorted(d_groups.items()):
		line = [sp] + ogs
		print('\t'.join(line), file=outTsv)


def pan_stats(OFdir, species):
	result = OrthoFinder(OFdir)
	species = parse_species(species, result)
	d_groups = {}
	for group in result.get_mcl_orthogroups(species):
		ogid = group.ogid
		for sp in set(group.species):
			try:
				d_groups[sp].add(ogid)
			except KeyError:
				d_groups[sp] = {ogid}
	d_groups = OrderedDict([(sp, d_groups[sp]) for sp in species])
	_flower_plot(d_groups, outfig='Flower.plot.pdf')

	groups = list(d_groups.values())
	data = []
	for i in range(1, len(species)+1):
		insects, dispens = [], []
		print(i, end=' ')
		for xgroups in comb(groups, i):
			insect, dispen = insect_groups(xgroups)
			insects += [len(insect)]
			dispens += [len(dispen)+len(insect)]
		print(len(insects))
		data += [[i] + stats_data(insects) + stats_data(dispens)]
	plot_pan(data, outfig='Core-Pan.plot.pdf')


def _flower_plot(d_groups, outfig):
	import matplotlib.pyplot as plt
	from flower_plot import flower_plot
	from Stairway_plot import create_colors
	colors = create_colors(len(d_groups)+2)
	insect, _ = insect_groups(list(d_groups.values()))
	n_core = len(insect)
	data = OrderedDict()
	for (sp, ids), color in zip(list(d_groups.items()), colors[2:]):
		others = [_ids for _sp, _ids in list(d_groups.items()) if _sp != sp]
		_ins, _dis = insect_groups(others)
		others = _ins | _dis
		shell = len(ids & others)
		unique = len(ids - others)
		data[sp] = dict(color=color, shell=shell, unique=unique)

	fig, ax = plt.subplots(
		subplot_kw={'aspect': 'equal'}, figsize=(9, 9), dpi=300)
	flower_plot(data, n_core=n_core,
				core_color=colors[0], shell_color=colors[1], alpha=0.3, ax=ax)
	plt.savefig(outfig, )


def plot_pan(data, outfig, alpha=0.2):
	import matplotlib.pyplot as plt
	data = np.array(data)
	x = data[:, 0]
	cc, cd = ('#1f77b4', '#ff7f0e',)
	yc, ycL, ycU = data[:, 1], data[:, 2], data[:, 3]
	yd, ydL, ydU = data[:, 4], data[:, 5], data[:, 6]
	plt.figure(figsize=(8, 6))
	plt.plot(x, yc, 'o-', label='Core', color=cc, )
	plt.fill_between(x, ycL, ycU, color=cc,  alpha=alpha)
	plt.plot(x, yd, 's-', label='Pan', color=cd, )
	plt.fill_between(x, ydL, ydU, color=cd,  alpha=alpha)
	plt.legend(loc='best')
	plt.xlabel('Genome number')
	plt.ylabel('Family number')
	plt.xticks(x, list(map(int, x)))
	plt.savefig(outfig)


def stats_data(data):
	return [np.median(data), np.percentile(data, 2.5), np.percentile(data, 97.5)]


def comb(groups, N, MAX_COMB=1000):
	i = 0
	for xgrp in itertools.combinations(groups, N):
		yield list(xgrp)
		i += 1
		if i > MAX_COMB:
			break


def insect_groups(groups):
	for i, ids in enumerate(groups):
		if i == 0:
			insect = ids
			continue
		insect = insect & ids
	dispens = set([])
	for ids in groups:
		dispens = dispens | (ids - insect)
	return insect, dispens


def to_astral(ResultsDir, pepSeq, outTrees, species=None,
			  tmpdir='/io/tmp/share', min_singlecopy=0.7):
	from RunCmdsMP import run_job
	tmpdir = '{}/to_astral.{}'.format(tmpdir, os.getpid())
	logger.info('change tmpdir to {}'.format(tmpdir))
	if not os.path.exists(tmpdir):
		mkdirs(tmpdir)
	species = parse_species(species)
	result = OrthoFinder(ResultsDir)
	if species is None:
		species = result.Species
	print(species, file=sys.stderr)
	d_seqs = seq2dict(pepSeq)
	cmd_list = []
	treefiles = []
	for og in result.get_orthogroups(sps=species):
		d_singlecopy = {genes[0]: sp for sp, genes in list(
			og.spdict.items()) if len(genes) == 1}
		singlecopy_ratio = 1.0*len(d_singlecopy) / len(species)
		if not singlecopy_ratio >= min_singlecopy:
			continue
		outSeq = '{}/{}.pep'.format(tmpdir, og.ogid)
		f = open(outSeq, 'w')
		for gene, sp in list(d_singlecopy.items()):
			try:
				rc = d_seqs[gene]
			except KeyError:
				continue
			rc.id = sp
			SeqIO.write(rc, f, 'fasta')
		f.close()
		alnSeq = outSeq + '.aln'
		alnTrim = alnSeq + '.trimal'
		treefile = alnTrim + '.treefile'
		cmd = '. ~/.bashrc; [ ! -s {} ]'.format(treefile)
		cmds = [cmd]
		cmd = 'mafft --auto {} > {}'.format(outSeq, alnSeq)
		cmds += [cmd]
		cmd = 'trimal -automated1 -in {} -out {} &> /dev/null'.format(
			alnSeq, alnTrim)
		cmds += [cmd]
		cmd = 'iqtree -s {} -bb 1000 -mset JTT -nt 1 &> /dev/null'.format(
			alnTrim)
		cmds += [cmd]
		cmds = ' && '.join(cmds)
		cmd_list += [cmds]
		treefiles += [treefile]
	cmd_list2 = []
	bin = 10
	for i in range(0, len(cmd_list), bin):
		cmd_list2 += ['\n'.join(cmd_list[i:i+bin])]
	cmd_file = '{}/cmds.list'.format(tmpdir)
	run_job(cmd_file, cmd_list=cmd_list2, tc_tasks=20)
	for treefile in treefiles:
		if not os.path.exists(treefile):
			logger.warn('{} do not exist, check the log file `{}`'.format(
				treefile, treefile.replace('.treefile', '.log')))
			continue
		for line in open(treefile):
			outTrees.write(line)


def retrieve_orthologs(ResultsDir, collinearity, gff, fout=sys.stdout):
	import networkx as nx
	from mcscan import Collinearity
	G = nx.Graph()  # OG Graph
	for group in OrthoFinder(ResultsDir).get_orthogroups():
		for g1, g2 in itertools.combinations(group.genes, 2):
			G.add_edge(g1, g2)

	orthologs = set([])
	for rc in Collinearity(collinearity, gff):
		d_chrom = rc.d_chrom  # chrom: [g1,g2,...]
		d_gene = rc.d_gene		# gene_id: gene
		chr1, chr2 = rc.chr1, rc.chr2
		genes1, genes2 = rc.genes1, rc.genes2  # [g1,g2,g3,...]
		last_g1, last_g2 = genes1[0], genes2[0]
		for g1, g2 in zip(genes1[1:], genes2[1:]):
			g1_idx = last_g1.index, g1.index
			g2_idx = last_g2.index, g2.index
			g1_start, g1_end = min(g1_idx), max(g1_idx)
			g2_start, g2_end = min(g2_idx), max(g2_idx)
			if g1_end - g1_start == 1 or g2_end-g2_start == 1:  # continuous
				continue
			inner_genes1 = d_chrom[chr1][g1_start+1:g1_end]
			inner_genes2 = d_chrom[chr2][g2_start+1:g2_end]
			for inner_g1, inner_g2 in itertools.product(inner_genes1, inner_genes2):
				inner_gid1, inner_gid2 = inner_g1.id, inner_g2.id
				if G.has_edge(inner_gid1, inner_gid2):
					orthologs.add((inner_gid1, inner_gid2))
			last_g1, last_g2 = g1, g2
	for ortholog in orthologs:
		print('\t'.join(ortholog), file=fout)
	print(len(orthologs), 'pairs retrieved.', file=sys.stderr)


def retrieve_allele(ResultsDir, collinearity, gff, fout=sys.stdout, min_block=10, win_size=20, diff_sp=True, min_score=100, sps=None):
	from mcscan import Collinearity
	import networkx as nx
	result = OrthoFinder(ResultsDir)
	sps = parse_species(sps, result)
	G = nx.Graph()  # OG Graph
	for group in result.get_orthogroups(sps):
		for g1, g2 in itertools.combinations(group.genes, 2):
			G.add_edge(g1, g2, OG=group.ogid)

	blastG = nx.Graph()  # blast Graph
	best_hits = {}
	for line in result.get_blast(sps):
		g1, g2 = line[:2]
		key = (g1, g2)
		score = float(line[-1])
		if key not in best_hits or (key in best_hits and best_hits[key] < score):
			best_hits[key] = score
	for (g1, g2), score in list(best_hits.items()):
		if blastG.has_edge(g1, g2):
			blastG.edge[g1][g2]['weight'] + 1
		else:
			blastG.add_edge(g1, g2, OG='none', weight=1)
	altG = nx.Graph() 	# allel Graph
	_sps = set([])
	blocks = []
	for rc in Collinearity(collinearity, gff):
		if rc.N < min_block:
			continue
		if not (rc.species1 in set(sps) and rc.species2 in set(sps)):
			continue
		if diff_sp and rc.species1 == rc.species2:
			continue
		_sps = _sps | set([rc.species1, rc.species2])
		blocks += [(rc.score, rc.N, rc.genes1, rc.genes2, rc.chr1,
					rc.chr2, rc.species1, rc.species2, rc.Alignment)]

	sps = _sps
	d_chrom = rc.d_chrom  # chrom: [g1,g2,...]
	d_genes = rc.d_gene		# gene.id: gene
	d_syn = {}
	blocks = sorted(blocks, reverse=1)
	for score, N, genes1, genes2, chr1, chr2, sp1, sp2, Alignment in blocks:
		key = (sp1, sp2)
		if key not in d_syn:
			d_syn[key] = set([])
		for g1, g2 in zip(genes1, genes2):
			id1, id2 = g1.id, g2.id
			if not id1 in d_syn[key] and not id2 in d_syn[key]:
				og = G.edge[g1.id][g2.id]['OG']
				altG.add_edge(
					g1.id, g2.id, source='MCScanX-{}'.format(Alignment), OG=og, weight=4)
				d_syn[key] = d_syn[key] | {g1.id, g2.id}
	# OG

	def mean_dist(inner_g1, inner_g2, g1_idx, g2_idx):
		import math
		return math.sqrt(abs(inner_g1.index - g1_idx) * abs(inner_g2.index - g2_idx))
	for score, N, genes1, genes2, chr1, chr2, sp1, sp2, Alignment in blocks:
		key = (sp1, sp2)
		for g1, g2 in zip(genes1, genes2):
			g1_idx = g1.index
			g2_idx = g2.index
			g1_start, g1_end = max(0, g1_idx-win_size/2), g1_idx+win_size/2
			g2_start, g2_end = max(0, g2_idx-win_size/2), g2_idx+win_size/2
			inner_genes1 = d_chrom[chr1][g1_start:g1_end+1]
			inner_genes2 = d_chrom[chr2][g2_start:g2_end+1]
			inner_pairs = list(itertools.product(inner_genes1, inner_genes2))
			inner_pairs = sorted(inner_pairs, key=lambda x: mean_dist(
				x[0], x[1], g1_idx, g2_idx))
			for inner_g1, inner_g2 in inner_pairs:
				id1, id2 = (inner_g1.id, inner_g2.id)
				if G.has_edge(id1, id2) and not altG.has_edge(id1, id2) and \
						not id1 in d_syn[key] and not id2 in d_syn[key]:  # reduce network
					og = G.edge[id1][id2]['OG']
					altG.add_edge(id1, id2, source='OG', OG=og, weight=2)
					d_syn[key] = d_syn[key] | {id1, id2}
	# blast
	for score, N, genes1, genes2, chr1, chr2, Alignment in blocks:
		for g1, g2 in zip(genes1, genes2):
			g1_idx = g1.index
			g2_idx = g2.index
			g1_start, g1_end = max(0, g1_idx-win_size/2), g1_idx+win_size/2
			g2_start, g2_end = max(0, g2_idx-win_size/2), g2_idx+win_size/2
			inner_genes1 = d_chrom[chr1][g1_start:g1_end+1]
			inner_genes2 = d_chrom[chr2][g2_start:g2_end+1]
			for inner_g1, inner_g2 in itertools.product(inner_genes1, inner_genes2):
				id1, id2 = (inner_g1.id, inner_g2.id)
				if blastG.has_edge(id1, id2) and not altG.has_edge(id1, id2) and not altG.has_node(id1) and not altG.has_node(id2):
					weight = blastG.edge[id1][id2]['weight'] * 0.5
					og = 'none'
					altG.add_edge(id1, id2, source='blast',
								  OG=og, weight=weight)
	# resolve repeats
	d_degrees = altG.degree(weight='weight')
	print(sorted(sps), file=sys.stderr)
	lines = []
	for cmpt in nx.connected_components(altG):
		sg = altG.subgraph(cmpt)
		grp = OrthoMCLGroupRecord(genes=cmpt)
		sp_dict = grp.spdict
		primary_genes = []
		alter_genes = []
		for sp in sorted(sps):
			genes = sp_dict.get(sp, [])
			g_pri = max(genes, key=lambda x: d_degrees[x]) if genes else '-'
			g_alt = set(genes) - set([g_pri])
			g_alt = [d_genes[x].raw_gene for x in g_alt]
			g_alt = ','.join(g_alt) if g_alt else '-'
			primary_genes += [g_pri]
			alter_genes += [g_alt]
		attrs = []
		for pri_g1, pri_g2 in itertools.combinations(primary_genes, 2):
			try:
				attr0 = altG.edge[pri_g1][pri_g2]
				attr = '{source}:{OG}'.format(**attr0)
				_sps = [d_genes[pri_g1].species, d_genes[pri_g2].species]
				line = [pri_g1, pri_g2] + _sps + [attr0['source'], attr0['OG']]
				print('\t'.join(line), file=sys.stderr)
			except KeyError:
				attr = '-'
			attrs += [attr]
		idxes = [d_genes[g].index for g in primary_genes if g != '-']
		gidx = round(1.0 * sum(idxes) / len(idxes), 2)
		idxes = [d_genes[g].ichr for g in primary_genes if g != '-']
		cidx = 1.0 * sum(idxes) / len(idxes)
		primary_genes = [d_genes[g].raw_gene if g !=
						 '-' else '-' for g in primary_genes]
		line = [cidx, gidx] + primary_genes + alter_genes + attrs
		lines += [line]
	# title
	col1 = col2 = sorted(sps)
	col3 = ['{}-{}'.format(sp1, sp2)
			for sp1, sp2 in itertools.combinations(sorted(sps), 2)]
	line = ['']*2 + ['# primary alleles'] + ['']*(len(col1)-1) + ['# secondary alleles'] + ['']*(len(col2)-1) + \
					['# sources of primary gene pairs'] + ['']*(len(col3)-1)
	print('\t'.join(line), file=fout)
	line = ['chrom', 'idx'] + col1 + col2 + col3
	print('\t'.join(line), file=fout)
	for line in sorted(lines):
		print('\t'.join(map(str, line)), file=fout)


def lazy_orthologs(ortholog):
	if os.path.isdir(ortholog+'/species_to_species_orthologs'):
		parser = SonicParanoid(ortholog)
	elif glob.glob('{}/table.*-*'.format(ortholog)):
		parser = SonicParanoid(ortholog, inparanoid=True)
	else:
		parser = OrthoFinder(ortholog)
	return parser


class SonicParanoid:
	def __init__(self, ResultsDir, inparanoid=False):
		self.ResultsDir = ResultsDir
		self.inparanoid = inparanoid
		if inparanoid:
			self.Orthologues = glob.glob('{}/table.*-*'.format(ResultsDir))
		else:
			self.Orthogroups = '{}/ortholog_groups/ortholog_groups.tsv'.format(
				ResultsDir)
			self.Orthologues = glob.glob(
				'{}/species_to_species_orthologs/*/*'.format(ResultsDir))

	def get_homologs(self, sps=None, sp1=None, sp2=None, **kargs):
		'''Obtain pairwise orthologs'''
		if sp1 is not None or sp2 is not None:
			sps = {sp1, sp2}
		if sps is not None:
			orthoFiles = []
			for sp1, sp2 in itertools.permutations(sps, 2):
				if self.inparanoid:
					orthoFile = '{}/table.{}-{}'.format(
						self.ResultsDir, sp1, sp2)
				else:
					orthoFile = '{}/species_to_species_orthologs/{}/{}-{}'.format(
						self.ResultsDir, sp1, sp1, sp2)
				orthoFiles += [orthoFile]
		else:
			orthoFiles = self.Orthologues
		ortho_pairs = set([])
		idx = 0
		for orthoFile in orthoFiles:
			idx += 1
			for line in SonicParanoidOrthologs(orthoFile):
				Genes_1 = line.OrthoA
				Genes_2 = line.OrthoB
				for g1, g2 in itertools.product(Genes_1, Genes_2):
					if (g2, g1) in ortho_pairs:
						continue
					ortho_pairs.add((g1, g2))  # orthologs
				for Genes in [Genes_1, Genes_2]:
					if len(Genes) < 2:
						continue
					for g1, g2 in itertools.combinations(Genes, 2):
						if (g2, g1) in ortho_pairs:
							continue
						ortho_pairs.add((g1, g2))  # inparalogs
		return ortho_pairs


def test_sonicparanoid(ResultsDir):
	for g1, g2 in SonicParanoid(ResultsDir).get_homologs():
		print(g1, g2)


class SonicParanoidOrthologs:
	def __init__(self, orthologs):
		self.orthologs = orthologs

	def __iter__(self):
		return self._parse()

	def _parse(self):
		i = 0
		for line in open(self.orthologs):
			i += 1
			if i == 1:
				continue
			yield SonicParanoidOrthologsLine(line)


class SonicParanoidOrthologsLine:
	def __init__(self, line):
		Size, Relations, OrthoA, OrthoB = line.strip().split('\t')
		self.OrthoA = self.parse_ortho(OrthoA)
		self.OrthoB = self.parse_ortho(OrthoB)

	def parse_ortho(self, Ortho):
		return [val for i, val in enumerate(Ortho.split()) if i % 2 == 0]


class OrthoFinder:
	def __init__(self, ResultsDir):
		self.ResultsDir = ResultsDir
		self.WorkingDirectory = '{}/WorkingDirectory/'.format(ResultsDir)
		self.SpeciesTreeAlignment_path = \
			'{}/MultipleSequenceAlignments/SpeciesTreeAlignment.fa'.format(
				ResultsDir)
		self.SpeciesIDs = '{}/WorkingDirectory/SpeciesIDs.txt'.format(
			ResultsDir)
		self.SpeciesTree_rooted = \
			'{}/Species_Tree/SpeciesTree_rooted_node_labels.txt'.format(
				ResultsDir)
		self.Orthologues = glob.glob(
			'{}/Orthologues/*/*__v__*.tsv'.format(ResultsDir))
		self.Duplications = '{}/Gene_Duplication_Events/Duplications.tsv'.format(
			ResultsDir)
		self.Orthogroups = '{}/Orthogroups/Orthogroups.tsv'.format(ResultsDir)
		self.MCLOrthogroups = '{}/Orthogroups/Orthogroups.txt'.format(
			ResultsDir)
		self.SequenceIDs = '{}/WorkingDirectory/SequenceIDs.txt'.format(
			ResultsDir)

	@property
	def SpeciesTreeAlignment(self):
		msa = self.SpeciesTreeAlignment_path
		if os.path.exists(msa):
			return msa
		# v3: no SpeciesTreeAlignment.fa; concatenate from Orthogroups_for_concatenated_alignment.txt
		og_file = '{}/Species_Tree/Orthogroups_for_concatenated_alignment.txt'.format(
			self.ResultsDir)
		if not os.path.exists(og_file):
			return msa  # will raise FileNotFoundError downstream
		aln_dir = '{}/MultipleSequenceAlignments/'.format(self.ResultsDir)
		og_ids = [l.strip() for l in open(og_file) if l.strip()]
		aln_files = ['{}{}.fa'.format(aln_dir, og) for og in og_ids]
		with open(msa, 'w') as out:
			catAln(aln_files, out)
		return msa

	@property
	def orthogroups(self):
		return Group(self.Orthogroups)

	def get_orthogroups(self, sps=None):
		for group in Group(self.Orthogroups):
			genes = list(group.get_group(sps=sps))
			yield OrthoMCLGroupRecord(ogid=group.ogid, genes=genes)

	def get_mcl_orthogroups(self, sps=None):
		for group in OrthoMCLGroup(self.MCLOrthogroups):
			genes = list(group.get_group(sps=sps))
			yield OrthoMCLGroupRecord(ogid=group.ogid, genes=genes)

	@property
	def OGDict(self):
		'''dictionary: gene id -> OG id '''
		d = {}
		for rc in Group(self.Orthogroups):
			for genes in rc.genes:
				for gene in genes:
					d[gene] = rc.ogid
		return d

	@property
	def GenesDict(self):
		'''OG id -> genes'''
		d = {}
		for group in Group(self.Orthogroups):
			d[group.ogid] = list(group.get_group())
		return d

	def get_oglogs(self):
		'''Get all orthologous relationships from self.Orthogroups.'''
		for rc in Group(self.Orthogroups):
			genes = list(rc.get_group())
			for g1, g2 in itertools.combinations(genes, 2):
				yield g1, g2

	@property
	def Single_Copy_Orthologue(self):
		'''Single_Copy OG'''
		Single_Copy_OGs = [line.strip() for line in
						   open('{}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt'.format(self.ResultsDir))]
		return Single_Copy_OGs

	def formatdir(self, xdir):
		return '{}/{}'.format(self.ResultsDir, xdir)

	def checkdir(self, xdir):
		fdir = self.formatdir(xdir)
		tgz = '{}.tgz'.format(fdir.rstrip('/'))
		tgz2 = '{}.tgz'.format(xdir.rstrip('/'))
		if os.path.exists(tgz) and not os.path.exists(fdir):
			cmd = 'cd {} && tar xzf {}'.format(self.ResultsDir, tgz2)
			print(cmd, file=sys.stderr)
			run_cmd(cmd)

	def Single_Copy_Codon_Align(self, cdsSeqs, tmpdir='/tmp'):
		'''Generate codon alignments for single-copy orthogroups.'''
		Single_Copy_OGs = self.Single_Copy_Orthologue
		singlecopy_seqdir = self.formatdir('Single_Copy_Orthologue_Sequences')
		self.checkdir('Single_Copy_Orthologue_Sequences')
		d_cds = seq2dict(cdsSeqs)
		ALNs = []
		for OG in Single_Copy_OGs:
			pepSeq = '{}/{}.fa'.format(singlecopy_seqdir, OG)
			cdsSeq = '{}/{}.cds'.format(tmpdir, OG)
			f = open(cdsSeq, 'w')
			for rc in SeqIO.parse(pepSeq, 'fasta'):
				SeqIO.write(d_cds[rc.id], f, 'fasta')
			f.close()
			pepAln = '{}/{}.pep.aln'.format(tmpdir, OG)
			cmd = 'mafft --auto {} > {} 2> /dev/null'.format(pepSeq, pepAln)
			os.system(cmd)
			pepTrim = pepAln + '.trimal'

			cdsAln = cdsSeq + '.aln'
			cmd = 'pal2nal.pl -output fasta {} {} > {} 2> /dev/null'.format(
				pepAln, cdsSeq, cdsAln)
			os.system(cmd)
			if os.path.getsize(cdsSeq) > 0 and os.path.getsize(cdsAln) == 0:
				print('Error in CMDS `{}`'.format(cmd), file=sys.stderr)
				continue
			cdsTrim = cdsAln + '.trimal'
			cmd = 'trimal -gt 0.8 -in {} -out {} &> /dev/null'.format(
				cdsAln, cdsTrim)
			os.system(cmd)
			ALNs += [cdsTrim]
		return ALNs

	def Single_Copy_Pep_Align(self):
		return ['{}/MultipleSequenceAlignments/{}.fa'.format(self.ResultsDir, OG)
				for OG in self.Single_Copy_Orthologue]

	@property
	def SequenceIDdict(self):
		'''Mapping of new IDs to original IDs for protein sequences.'''
		d = {}
		for line in open(self.SequenceIDs):
			if line.startswith('#'):
				continue
			id, geneName = line.strip().split(': ')
			d[id] = geneName
		return d

	def get_SequenceDict(self):
		'''Mapping relationships between protein sequence IDs and records.'''
		self.Sequences = glob.glob(
			'{}/WorkingDirectory/Species*.fa'.format(self.ResultsDir))
		d_seq_ids = self.SequenceIDdict
		d_seqs = {}
		for Sequence in self.Sequences:
			for rc in SeqIO.parse(Sequence, 'fasta'):
				rc.id = d_seq_ids[rc.id]
				d_seqs[rc.id] = rc
		return d_seqs

	def get_Blast(self, sps=None, fout=sys.stdout):
		for _, _, temp in self.get_blast(sps=sps):
			print('\t'.join(temp), file=fout)

	def to_wgdi(self, species=None, outdir='wgdi', split=True):
		sps = parse_species(species)
		d_handle = {}
		for sp1, sp2, temp in self.get_blast(sps=sps):
			if (sp1, sp2) in d_handle:
				pass
			else:
				out = '{}/{}-{}.blast'.format(
					outdir, sp1, sp2) if split else '{}/{}.blast'.format(outdir, 'all')
				if not split and d_handle:
					d_handle[(sp1, sp2)] = list(d_handle.values())[0]
				else:
					d_handle[(sp1, sp2)] = open(out, 'w')
			handle = d_handle[(sp1, sp2)]
			print('\t'.join(temp), file=handle)
		for handle in list(d_handle.values()):
			handle.close()

	def get_blast(self, sps=None):
		d_sp0 = self.SpeciesIDdict
		d_sp = self.reverse_SpeciesIDdict
		d_seq = self.SequenceIDdict
		if sps is None:
			spIds = list(d_sp.values())
		else:
			nofound = [sp for sp in sps if sp not in d_sp]
			if nofound:
				print('{} is not found'.format(nofound), file=sys.stderr)
			spIds = [d_sp.get(sp) for sp in sps if sp in d_sp]
		spIds = sorted(spIds)
		for sp1, sp2 in itertools.product(spIds, spIds):
			blast = '{}/Blast{}_{}.txt.gz'.format(
				self.WorkingDirectory, sp1, sp2)
			for line in open(blast):
				temp = line.strip().split('\t')
				temp[0] = d_seq[temp[0]]
				temp[1] = d_seq[temp[1]]
				yield d_sp0[sp1], d_sp0[sp2], temp

	@property
	def SpeciesIDdict(self):
		'''Mapping relationship between new species IDs and names.'''
		d = {}
		for line in open(self.SpeciesIDs):
			if line.startswith('#'):
				continue
			temp = line.strip().split()
			id, spName = temp[:2]
			spName = '.'.join(spName.split('.')[:-1])
			id = id.strip(':')
			d[id] = spName
		return d

	@lazyproperty
	def Species(self):
		'''List of species names.'''
		species = []
		for line in open(self.SpeciesIDs):
			if line.startswith('#'):
				continue
			temp = line.strip().split()
			id, spName = temp[:2]
			spName = '.'.join(spName.split('.')[:-1])
			species += [spName]
		return species

	@property
	def reverse_SpeciesIDdict(self):
		'''Dictionary of species names and new IDs.'''
		return dict([(spName, id) for id, spName in list(self.SpeciesIDdict.items())])

	def spName2Id(self, *sps):
		'''Get the new ID of the specified species.'''
		d_sp = self.reverse_SpeciesIDdict
		return [d_sp[sp] for sp in sps]

	def spId2Name(self, *sps):
		'''Get the name of the specified species ID.'''
		d_sp = self.SpeciesIDdict
		return [d_sp[sp] for sp in sps]

	def get_species(self, genes, sep='|'):
		'''Retrieve species name from gene name, inherited from the OrthoMCL format.'''
		return [gene.split(sep)[0] for gene in genes]

	def get_species_specific(self, sp):
		'''Get all genes that are unique to the specified species.'''
		i = 0
		d_genes = {}
		all_sp_genes = []
		for line in open(self.Orthogroups):
			i += 1
			temp = line.rstrip('\r\n').split('\t')
			if i == 1:
				spidx = temp.index(sp)
				continue
			od_id = temp[0]
			sp_genes = temp[spidx]
			temp.pop(spidx)
			genes = temp[1:]
			sp_genes = sp_genes.split(', ')
			all_sp_genes += sp_genes
			if not (sp_genes and set(genes) == set([''])):
				continue
			for gene in sp_genes:
				d_genes[gene] = od_id
		all_sp_genes = set(all_sp_genes)
		sp_id0 = self.reverse_SpeciesIDdict[sp]
		for seq_id, gene in list(self.SequenceIDdict.items()):
			sp_id, sid = seq_id.split('_')
			if sp_id == sp_id0 and gene not in all_sp_genes:
				d_genes[gene] = None
		return d_genes

	def classify_genes(self, sps=None, fout=sys.stdout):
		if sps is None:
			sps = set(self.Species)
		d_cls = OrderedDict([(sp, {}) for sp in sps])
		catlogs = set([])
		for og in OrthoMCLGroup(self.MCLOrthogroups):
			for sp, catlog, count in og.classify():
				if not sp in sps:
					continue
				try:
					d_cls[sp][catlog] += count
				except KeyError:
					d_cls[sp][catlog] = count
				catlogs.add(catlog)
		catlogs = sorted(catlogs)
		line = ['Species'] + [cl for _, cl in catlogs]
		print('\t'.join(line), file=fout)
		for sp in sps:
			count = [d_cls[sp].get(catlog, 0) for catlog in catlogs]
			line = [sp] + count
			line = list(map(str, line))
			print('\t'.join(line), file=fout)

	def get_species_specific2(self, sp, ex_sps=[]):
		'''Get all genes that are unique to the specified species.'''
		if not ex_sps:
			ex_sps = set(self.Species) - set([sp])
		d_specific = {}
		ogs, genes = 0, 0
		mogs, mgenes = 0, 0
		for og in OrthoMCLGroup(self.MCLOrthogroups):
			d_count = og.counter
			d_genes = og.spdict
			sp_ngene = d_count.get(sp, 0)
			ex_sps_ngene = [d_count.get(ex_sp, 0) for ex_sp in ex_sps]
			# specific: >0 in this and all 0 in others
			if sp_ngene > 0 and not any(ex_sps_ngene):
				ogs += 1
				genes += sp_ngene
				if sp_ngene > 1:
					mogs += 1
					mgenes += sp_ngene
				for gene in d_genes[sp]:
					d_specific[gene] = og.ogid
		print('exclude species: {}\ntotal specific OGs: {}\ntotal specific genes:{}\n\
specific multi-copy OGs: {}\nspecific multi-copy genes: {}'.format(
			len(ex_sps), ogs, genes, mogs, mgenes), file=sys.stderr)
		return d_specific

	def count_og(self, sps):
		ogs, genes = 0, 0
		for og in Group(self.Orthogroups):
			d_count = og.counter
			count = [d_count[sp] for sp in sps]
			if any(count):
				ogs += 1
				genes += sum(count)
		return ogs, genes

	def get_blast_files(self, sp1, sp2, byName=True):
		'''Retrieve the BLAST files inter- or intra-species.'''
		if byName:
			sp1, sp2 = self.spName2Id(sp1, sp2)
		blast_files = []
		if sp1 == sp2:
			SPs = [(sp1, sp2)]
		else:
			SPs = [(sp1, sp2), (sp2, sp1)]
		for sp1, sp2 in SPs:
			blast_files += ['{}/WorkingDirectory/Blast{}_{}.txt.gz'.format(
				self.ResultsDir, sp1, sp2)]
		return blast_files

	def get_paralogs(self, sp=None, byName=True, min_support=0.5):
		'''Retrieve paralogs from Gene_Duplication_Events/'''
		if sp is not None and not byName:
			sp, = self.spId2Name(sp)
		Duplications = self.Duplications
		para_pairs = set([])
		i = 0
		for line in open(Duplications):
			i += 1
			if i == 1:
				continue
			temp = line.strip().split('\t')
			Orthogroup, Species_Trer_Node, Gene_Tree_Node, Support, Type, Genes_1, Genes_2 = temp
			Support = float(Support)
			if Support < min_support:
				continue
			Genes_1 = Genes_1.split(', ')
			Genes_2 = Genes_2.split(', ')
			Genes_1 = list(map(gene_format_p, Genes_1))
			Genes_2 = list(map(gene_format_p, Genes_2))
			for (sp1, g1), (sp2, g2) in itertools.product(Genes_1, Genes_2):
				if sp1 == sp2:
					if sp is not None and sp != sp1:
						continue
					if (g2, g1) in para_pairs:
						continue
					para_pairs.add((g1, g2))
		return para_pairs

	def get_paralogs2(self, sps=None, sp=None, byName=True):
		'''Retrieve paralogs from Orthologues/'''
		if len(sps) == 1 and sp is None:
			sp = sps[0]
		if sp:
			if not byName:
				sp, = self.spId2Name([sp])
			orthoFiles = glob.glob(
				'{}/Orthologues/Orthologues_{}/*.tsv'.format(self.ResultsDir, sp))
		elif sps:
			orthoFiles = []
			for sp1, sp2 in itertools.permutations(sps, 2):
				if not byName:
					sp1, sp2 = self.spId2Name(sp1, sp2)
				orthoFiles += ['{}/Orthologues/Orthologues_{}/{}__v__{}.tsv'.format(
					self.ResultsDir, sp1, sp1, sp2)]
		else:
			orthoFiles = self.Orthologues
		para_pairs = set([])
		idx = 0
		for orthoFile in orthoFiles:
			idx += 1
			i = 0
			for line in open(orthoFile):
				i += 1
				if i == 1:
					continue
				temp = line.strip().split('\t')
				Orthogroup, Genes_1, Genes_2 = temp
				Genes_1 = Genes_1.split(', ')
				Genes_2 = Genes_2.split(', ')
				Genes_1 = list(map(gene_format_o, Genes_1))
				Genes_2 = list(map(gene_format_o, Genes_2))
				for Genes in [Genes_1, Genes_2]:
					if len(Genes) < 2:
						continue
					for (sp1, g1), (sp2, g2) in itertools.combinations(Genes, 2):
						assert sp1 == sp2
						if sp is not None and sp != sp1:
							continue
						if (g2, g1) in para_pairs:
							continue
						para_pairs.add((g1, g2))
		return para_pairs

	def get_paralogs3(self):
		'''Retrieve paralogs from Orthogroups/Orthogroups.tsv'''
		para_pairs = set([])
		for rc in Group(self.Orthogroups):
			for genes in rc.genes:
				for g1, g2 in itertools.combinations(genes, 2):
					para_pairs.add((g1, g2))
		return para_pairs

	def get_orthologs(self, sps=None, sp1=None, sp2=None, byName=True):
		'''Retrieve orthologous pairs'''
		if sps is not None:
			orthoFiles = []
			for sp1, sp2 in itertools.permutations(sps, 2):
				if not byName:
					sp1, sp2 = self.spId2Name(sp1, sp2)
				orthoFiles += ['{}/Orthologues/Orthologues_{}/{}__v__{}.tsv'.format(
					self.ResultsDir, sp1, sp1, sp2)]
		elif sp1 is not None or sp2 is not None:
			if not byName:
				sp1, sp2 = self.spId2Name(sp1, sp2)
			orthoFiles = ['{}/Orthologues/Orthologues_{}/{}__v__{}.tsv'.format(
				self.ResultsDir, sp1, sp1, sp2),
				'{}/Orthologues/Orthologues_{}/{}__v__{}.tsv'.format(
				self.ResultsDir, sp2, sp2, sp1),]
		else:
			orthoFiles = self.Orthologues
		ortho_pairs = set([])
		idx = 0
		for orthoFile in orthoFiles:
			idx += 1
			i = 0
			for line in open(orthoFile):
				i += 1
				if i == 1:
					continue
				temp = line.strip().split('\t')
				try:
					Orthogroup, Genes_1, Genes_2 = temp
				except ValueError:
					print('Error:', orthoFile, temp, file=sys.stderr)
					continue
				Genes_1 = Genes_1.split(', ')
				Genes_2 = Genes_2.split(', ')
				Genes_1 = list(map(gene_format_o, Genes_1))
				Genes_2 = list(map(gene_format_o, Genes_2))
				for (sp1, g1), (sp2, g2) in itertools.product(Genes_1, Genes_2):
					assert sp1 != sp2
					if (g2, g1) in ortho_pairs:
						continue
					ortho_pairs.add((g1, g2))
		return ortho_pairs

	def get_homologs(self, sps=None, sp1=None, sp2=None, byName=True, **kargs):
		'''Retrieve orthologs and inparalogs'''
		if sp1 is not None or sp2 is not None:
			sps = {sp1, sp2}
		if sps is not None:
			orthoFiles = []
			for sp1, sp2 in itertools.permutations(sps, 2):
				if not byName:
					sp1, sp2 = self.spId2Name(sp1, sp2)
				orthoFiles += ['{}/Orthologues/Orthologues_{}/{}__v__{}.tsv'.format(
					self.ResultsDir, sp1, sp1, sp2)]
		else:
			orthoFiles = self.Orthologues
		ortho_pairs = set([])
		idx = 0
		for orthoFile in orthoFiles:
			idx += 1
			i = 0
			for line in open(orthoFile):
				i += 1
				if i == 1:
					continue
				temp = line.strip().split('\t')
				try:
					Orthogroup, Genes_1, Genes_2 = temp
				except ValueError:
					print('Error:', orthoFile, temp, file=sys.stderr)
					continue
				Genes_1 = Genes_1.split(', ')
				Genes_2 = Genes_2.split(', ')
				Genes_1 = list(map(gene_format_o, Genes_1))
				Genes_2 = list(map(gene_format_o, Genes_2))
				for (sp1, g1), (sp2, g2) in itertools.product(Genes_1, Genes_2):
					if sp1 is not None and sp2 is not None:
						assert sp1 != sp2
					if (g2, g1) in ortho_pairs:
						continue
					ortho_pairs.add((g1, g2))  # orthologs
				for Genes in [Genes_1, Genes_2]:
					if len(Genes) < 2:
						continue
					for (sp1, g1), (sp2, g2) in itertools.combinations(Genes, 2):
						if sp1 is not None and sp2 is not None:
							assert sp1 == sp2
						if sps is not None and sp1 not in sps:
							continue
						if (g2, g1) in ortho_pairs:
							continue
						ortho_pairs.add((g1, g2))  # inparalogs
		return ortho_pairs

	def get_orthologs_cluster(self, **kargs):
		'''Retrieve the clusters of orthologs, which are a subset of OGs.'''
		import networkx as nx
		G = nx.Graph()
		for g1, g2 in self.get_orthologs(**kargs):
			G.add_edge(g1, g2)
		i = 0
		for cmpt in nx.connected_components(G):
			i += 1
			ogid = 'OG_{}'.format(i)
			yield OrthoMCLGroupRecord(genes=cmpt, ogid=ogid)

	def get_singlecopy_orthologs(self, **kargs):
		for cmpt in self.get_orthologs_cluster(**kargs):
			if self.is_singlecopy(cmpt):
				yield cmpt

	def is_singlecopy(self, cmpt):
		sps = [gene_format_o(g)[0] for g in cmpt]
		if len(sps) == len(set(sps)):
			return True
		else:
			return False

	def get_singlecopy_orthologs2(self, sps, singlecopy_ratio=1):
		'''Retrieve the single-copy orthologs for the target species set.'''
		nsp = len(sps)
		for og in self.get_orthologs_cluster(sps=sps):
			og.nsp = nsp
			if og.singlecopy_ratio >= singlecopy_ratio:
				yield og

	def is_singlecopy2(self, cmpt, sps):
		if sps is None:
			sps = self.Species
		sps2 = [gene_format_o(g)[0] for g in cmpt]
		return len(sps2) == len(set(sps2)) and len(sps2) == len(sps)

	def convert_seq_id(self, seqfile, out_seqfile,  fmt='fasta'):
		'''convert species ID to species NAME in sequence file'''
		d_species = self.SpeciesIDdict
		if out_seqfile == self.SpeciesTreeAlignment:
			raise ValueError('stop to write {}'.format(out_seqfile))
		f = open(out_seqfile, 'w')
		for rc in SeqIO.parse(seqfile, fmt):
			rc.id = d_species[rc.id]
			rc.description += ' {} sites'.format(len(rc.seq))
			SeqIO.write(rc, f, fmt)

	def get_root(self, treefile, fmt='newick'):
		'''Retrieve the root from the tree file.'''
		tree = Phylo.read(treefile, fmt)
		root = tree.root
		if root.name == 'N0':
			for clade in root.clades:
				if not clade.name == 'N1':
					return clade.name
		return root.name

	def re_root(self, treefile, root, out_treefile, fmt='newick'):
		'''Reroot the tree using Phylo.'''
		tree = Phylo.read(treefile, fmt)
		for clade in tree.find_clades(root):
			if root == clade.name:
				root = clade
				break
		else:
			print('root {} is not found'.format(root), file=sys.stderr)
		tree.root_with_outgroup(root)
		Phylo.write(tree, out_treefile, fmt)

	def get_aln_len(self, alnfile, fmt='fasta'):
		'''Get the length of the alignment.'''
		for rc in SeqIO.parse(alnfile, fmt):
			return len(rc.seq)


def count_og(OFdir, species):
	species = parse_species(species)
	result = OrthoFinder(OFdir)
	ogs, genes = result.count_og(species)
	print('number of species: {}\nnumber of OGs: {}\nnumber of genes: {}'.format(
		len(species), ogs, genes))


def single_copy_stats(OFdir, species=None, ):
	species = parse_species(species)
	result = OrthoFinder(OFdir)
	if species is None:
		sps = species = result.Species
	else:
		sps = species
	groups = [group.singlecopy_ratio for group in result.get_singlecopy_orthologs2(
		sps=sps, singlecopy_ratio=0)]
	for percent in range(0, 105, 5):
		min_ratio = percent/1e2
		groups = [ratio for ratio in groups if ratio >= min_ratio]
		line = [percent, len(groups)]
		line = list(map(str, line))
		print('\t'.join(line), file=sys.stdout)


def single_copy_cds_align2(OFdir, pepSeqs, cdsSeqs, outALN, species=None,
						   singlecopy_ratio=1, tmpdir='./tmp', ncpu=20, mode='grid'):
	'''Retrieve the codon alignments of single-copy orthologs for the target species set.'''
	if not os.path.exists(tmpdir):
		os.mkdir(tmpdir)
	species = parse_species(species)
	result = OrthoFinder(OFdir)
	if species is None:
		sps = species = result.Species
	else:
		sps = species
	d_pep = seq2dict(pepSeqs)
	d_cds = seq2dict(cdsSeqs)
	cmd_list = []
	alnfiles = []
	for group in result.get_singlecopy_orthologs2(sps=sps, singlecopy_ratio=singlecopy_ratio):
		og = group.ogid
		outPep = '{}/{}.pep'.format(tmpdir, og)
		outCds = '{}/{}.cds'.format(tmpdir, og)
		f1 = open(outPep, 'w')
		f2 = open(outCds, 'w')
		for gene in group.singlecopy_genes:
			rc = d_pep[gene]
			rc.id, g = gene_format_common(gene)  # pep id -> taxon
			SeqIO.write(rc, f1, 'fasta')
			try:
				rc = d_cds[gene]
			except KeyError as e:
				print(e, file=sys.stderr)
				continue
			rc.id, g = gene_format_common(gene)  # cds
			SeqIO.write(rc, f2, 'fasta')
		f1.close()
		f2.close()

		cmds = []
		alnPep = outPep + '.aln'
		cmd = '. ~/.bashrc; mafft --auto {} > {} '.format(outPep, alnPep)
		cmds += [cmd]
		alnPepTrim = alnPep + '.trimal'
		cmd = 'trimal -in {} -out {} -gt 0.8'.format(alnPep, alnPepTrim)
		alnCds = outCds + '.aln'
		cmd = 'pal2nal.pl -output fasta {} {} > {} 2> /dev/null'.format(
			alnPep, outCds, alnCds)
		cmds += [cmd]
		alnCodon = alnCds + '.trimal'
		cmd = 'trimal -in {} -out {} -gt 0.8'.format(alnCds, alnCodon)
		cmds += [cmd]
		cmds = ' && '.join(cmds)
		cmd_list += [cmds]
		alnfiles += [alnCodon]
	cmd_file = '{}/cmds.list'.format(tmpdir)
	run_job(cmd_file, cmd_list=cmd_list,
			tc_tasks=ncpu, mode=mode, fail_exit=False)
	alnfiles2 = [alnfile for alnfile in alnfiles if os.path.exists(alnfile)]
	non_exists = set(alnfiles) - set(alnfiles2)
	if non_exists:
		print('{} not exists, check please.'.format(
			non_exists), file=sys.stderr)
	catAln(alnfiles2, outALN)


def collinear_og_trees(OFdir, collinearity, taxon, seqfile, tmpdir='/io/tmp/share', min_N=10):
	from mcscan import Collinearity
	tmpdir = '/io/tmp/share/collinear.15867/'
	mkdirs(tmpdir)
	result = OrthoFinder(OFdir)
	d_seqs = seq2dict(seqfile)
	d_ogdict = result.OGDict
	# Retrieve the OGs containing paralogous genesfor the specified species.
	ogs = []
	for rc in Collinearity(collinearity):
		if rc.N < min_N:
			continue
		if not rc.species1 == rc.species2 == taxon:
			continue
		for g1, g2 in rc.pairs:
			og1, og2 = d_ogdict[g1], d_ogdict[g2]
			assert og1 == og2
			ogs += [og1]
	# Construct trees for OGs.
	d_genedict = result.GenesDict
	cmd_list = []
	for og in set(ogs):
		genes = d_genedict[og]
		outSeq = '{}/{}.pep'.format(tmpdir, og)
		sps = []
		for gene in genes:
			sp, g = gene_format_common(gene)
			sps += [sp]
		mono = {'Ananas_comosus', 'Brachypodium_distachyon',
				'Musa_acuminata', 'Oryza_sativa'}
		eudi = {'Aquilegia_coerulea', 'Arabidopsis_thaliana', 'Coffea_canephora',
				'Nelumbo_nucifera', 'Populus_trichocarpa', 'Solanum_lycopersicum',
				'Tetracendron_sinense', 'Theobroma_cacao',
				'Trochodendron_aralioides'}
		magn = {'Cinnamomum_kanehirae',
				'Liriodendron_chinense', 'Litsea_cubeba'}
		sps = set(sps)
		if not (sps & mono and sps & eudi and sps & magn):
			continue
		f = open(outSeq, 'w')
		for gene in genes:
			sp, g = gene_format_common(gene)
			if sp in {'Piper_nigrum'}:
				continue
			rc = d_seqs[gene]
			SeqIO.write(rc, f, 'fasta')
		f.close()
		alnSeq = outSeq + '.aln'
		alnTrim = alnSeq + '.trimal'
		treefile = alnTrim + '.treefile'
		cmd = '[ ! -s {} ]'.format(treefile)
		cmds = [cmd]
		cmd = 'mafft --auto {} > {} 2> /dev/null'.format(outSeq, alnSeq)
		cmds += [cmd]
		cmd = 'trimal -automated1 -in {} -out {} &> /dev/null'.format(
			alnSeq, alnTrim)
		cmds += [cmd]
		cmd = 'iqtree -s {} -bb 1000 -mset JTT -nt 1 &> /dev/null'.format(
			alnTrim)
		cmds += [cmd]
		cmds = ' && '.join(cmds)
		cmd_list += [cmds]
		print(treefile)
	cmd_file = '{}/cmds.list'.format(tmpdir)


def get_all_homo_pairs(OFdir, outpairs, species=None,):
	species = parse_species(species)
	sps = species
	result = OrthoFinder(OFdir)
	for group in result.get_orthologs_cluster(sps=sps):
		for pair in itertools.combinations(group, 2):
			print('{}\t{}'.format(*pair), file=outpairs)


def to_codeml(OFdir, outDir, pepSeq, cdsSeq, species=None, min_species=0.7, min_seqs=4, singlecopy=False):
	'''Prepare positive selection files without requiring single-copy; 
	therefore, gene trees are also used.'''
	print(vars(), file=sys.stderr)
	if not os.path.exists(outDir):
		os.mkdir(outDir)
	species = parse_species(species)
	result = OrthoFinder(OFdir)
	if species is None:
		sps = None
		species = result.Species
	else:
		sps = set(species)
	d_pep = seq2dict(pepSeq)
	d_cds = seq2dict(cdsSeq)
	cmd_list = []
	seqfiles, treefiles = [], []
	outGroup = '{}/groups.txt'.format(outDir)
	f = open(outGroup, 'w')
	for group in result.get_orthologs_cluster(sps=sps):
		if len(group) < min_seqs:
			continue
		g_species = result.get_species(group)
		if 1.0*len(set(g_species)) / len(set(species)) < min_species:
			continue
		is_singlecopy = len(species) == len(g_species) == len(set(g_species))
		if singlecopy and not is_singlecopy:
			continue
		og = group.ogid
		print('{}: {}'.format(og, ' '.join(sorted(group))), file=f)

		# prepare sequences
		outPep = '{}/{}.pep'.format(outDir, og)
		outCds = '{}/{}.cds'.format(outDir, og)
		f1 = open(outPep, 'w')
		f2 = open(outCds, 'w')
		for gene in group:
			rc = d_pep[gene]
			rc.id = format_id_for_iqtree(gene)  # pep id	-> taxon
			SeqIO.write(rc, f1, 'fasta')
			rc = d_cds[gene]
			rc.id = format_id_for_iqtree(gene)  # cds
			SeqIO.write(rc, f2, 'fasta')
		f1.close()
		f2.close()

		# prepare commands
		cmds = []
		alnPep = outPep + '.aln'
		# align
		cmd = 'mafft --auto {} > {} 2> /dev/null'.format(outPep, alnPep)
		cmds += [cmd]
		seqfile = alnCds = outCds + '.aln'  # seqfile for codeml
		cmd = 'pal2nal.pl -output fasta {} {} | cut -f1 -d " " >{} 2> /dev/null'.format(
			alnPep, outCds, alnCds)
		cmds += [cmd]
		alnTrim = alnPep + '.trimal'
		cmd = 'trimal -automated1 -in {} -out {} &> /dev/null'.format(
			alnPep, alnTrim)
		cmds += [cmd]
		cmd = 'iqtree -s {} -mset JTT -nt 1 &> /dev/null'.format(alnTrim)
		cmds += [cmd]
		iqtreefile = alnTrim + '.treefile'
		treefile = unrooted_treefile = alnTrim + '.tre'  # treefile for codeml
		cmd = 'nw_topology {}  > {}'.format(
			iqtreefile, unrooted_treefile)  # nw_reroot -d -
		cmds += [cmd]
		cmds = ' && '.join(cmds)
		cmd_list += [cmds]
		seqfiles += [seqfile]
		treefiles += [treefile]

	f.close()
	# run
	cmd_list2 = []
	bin = 40
	for i in range(0, len(cmd_list), bin):
		cmd_list2 += ['\n'.join(cmd_list[i:i+bin])]
	cmd_file = '{}/cmds.list'.format(outDir)
	run_job(cmd_file, cmd_list=cmd_list2, tc_tasks=50)
	non_seqfiles = [
		seqfile for seqfile in seqfiles if not exists_and_size_gt_zero(seqfile)]
	non_treefiles = [
		treefile for treefile in treefiles if not exists_and_size_gt_zero(treefile)]
	for non_exists in [non_seqfiles, non_treefiles]:
		if non_exists:
			print('{} not exists, check.'.format(non_exists), file=sys.stderr)


def format_id_for_iqtree(id):
	return re.compile(r'[^\w]+').sub('_', id)


def format_for_iqtree(inSeq, outSeq):
	for rc in SeqIO.parse(inSeq, 'fasta'):
		sp, _ = gene_format_common(rc.id)
		rc.id = format_id_for_iqtree(rc.id)
		print('{}\t{}'.format(rc.id, sp), file=sys.stderr)
		SeqIO.write(rc, outSeq, 'fasta')


def exists_and_size_gt_zero(FILE):
	return os.path.exists(FILE) and os.path.getsize(FILE) > 0


def to_paml(OFdir, outDir, cdsSeq, species=None, singlecopy=True):
	'''Generate files required for PAML positive selection analysis; 
	select common genes; 
	for multi-copy genes, choose the one with the shortest branch length; 
	if singlecopy=True, then only completely single-copy genes are selected.'''
	print(vars(), file=sys.stderr)
	if not os.path.exists(outDir):
		os.mkdir(outDir)
	species = parse_species(species)
	result = OrthoFinder(OFdir)
	if species is None:
		sps = None
		species = result.Species
	else:
		sps = species
	d_seqs = result.get_SequenceDict()
	i, j, s = 0, 0, 0
	groups = []
	for genes in result.get_orthologs_cluster(sps=sps):
		i += 1
		g_species = result.get_species(genes)
		if not (set(species) & set(g_species) == set(species)):
			continue
		j += 1
		is_singlecopy = len(species) == len(g_species)
		if singlecopy and not is_singlecopy:
			continue
		if is_singlecopy:  # single copy
			s += 1
			groups += [genes]
			continue
		groups += [select_genes_bytree(genes, d_seqs, j)]
	print('total {} groups, shared {}, single copy {}'.format(i, j, s))
	i = 0
	d_cds = seq2dict(cdsSeq)
	outGroup = '{}/groups.txt'.format(outDir)
	f = open(outGroup, 'w')
	for group in groups:
		i += 1
		og = 'OG_{}'.format(i)
		print('{}: {}'.format(og, ' '.join(group)), file=f)
		outSeq = '{}/{}.pep'.format(outDir, og)
		outCds = '{}/{}.cds'.format(outDir, og)
		f1 = open(outSeq, 'w')
		f2 = open(outCds, 'w')
		for gene in group:
			rc = d_seqs[gene]
			rc.id, g = gene_format_common(gene)  # pep
			SeqIO.write(rc, f1, 'fasta')
			rc = d_cds[gene]
			rc.id, g = gene_format_common(gene)  # cds
			SeqIO.write(rc, f2, 'fasta')
		f1.close()
		f2.close()

		alnSeq = outSeq + '.aln'
		cmd = 'mafft --auto {} > {} 2> /dev/null'.format(outSeq, alnSeq)
		os.system(cmd)
		cdsSeq = outCds + '.aln'
		cmd = 'pal2nal.pl -output fasta {} {} > {} 2> /dev/null'.format(
			alnSeq, outCds, cdsSeq)
		os.system(cmd)
		cds_aln = cdsSeq + '.trimal'
		cmd = 'trimal -in {} -out {} -gt 0.8 -phylip'.format(cdsSeq, cds_aln)
		os.system(cmd)
	f.close()


def seq2dict(seq):
	return dict([(format_of_id(rc.id), rc)for rc in SeqIO.parse(seq, 'fasta')])


def format_of_id(id):
	return re.compile(r'[\(\)]').sub('_', id)


def select_genes_bytree(genes, d_seqs, i):
	'''Select one gene from multi-copy genes based on the shortest branch length in the tree.'''
	outSeq = '/io/tmp/share/xxx{}.fa'.format(i)
	with open(outSeq, 'w') as f:
		for gene in genes:
			rc = d_seqs[gene]
			rc.id = rc.id.replace('|', '-')
			SeqIO.write(rc, f, 'fasta')
	alnSeq = outSeq + '.aln'
	print("if [ $SGE_TASK_ID -eq {} ]; then".format(i))
	cmd = 'mafft --auto {} > {} 2> /dev/null'.format(outSeq, alnSeq)
	print(cmd)
	cmd = 'iqtree -s {} -pre {} -nt AUTO &> /dev/null'.format(alnSeq, alnSeq)
	print(cmd)
	print('fi')
	treefile = alnSeq + '.treefile'
	d_dist = {}
	tree = Phylo.read(treefile, 'newick')
	for gene1 in genes:
		sp1 = gene1.split('|')[0]
		gene0 = gene1
		gene1 = gene1.replace('|', '-')
		dist = 0
		for gene2 in genes:
			sp2 = gene2.split('|')[0]
			gene2 = gene2.replace('|', '-')
			if sp1 == sp2:
				continue
			dist += tree.distance(gene1, gene2)
		d_dist[gene0] = dist
	node0 = min(genes, key=lambda x: d_dist[x])
	sp0 = node0.split('|')[0]
	nodes = [node0]
	d_dist = {}
	node0 = node0.replace('|', '-')
	for gene1 in genes:
		sp1 = gene1.split('|')[0]
		gene0 = gene1
		if sp1 == sp0:
			continue
		gene1 = gene1.replace('|', '-')
		dist = tree.distance(gene1, node0)
		try:
			d_dist[sp1][gene0] = dist
		except KeyError:
			d_dist[sp1] = {gene0: dist}
	for sp, d_sp_dist in list(d_dist.items()):
		node = min(list(d_sp_dist.keys()), key=lambda x: d_sp_dist[x])
		nodes += [node]
	print(sorted(genes), sorted(nodes))
	return nodes


def single_copy_cds_align(OFdir, cdsSeqs, outALN, tmpdir='./tmp'):
	'''Concatenate the codon alignments of single-copy genes.'''
	if not os.path.exists(tmpdir):
		os.mkdir(tmpdir)

	result = OrthoFinder(OFdir)
	inALNs = result.Single_Copy_Codon_Align(cdsSeqs, tmpdir=tmpdir)
	catAln(inALNs, outALN)


def single_copy_pep_align(OFdir, outALN, ):
	result = OrthoFinder(OFdir)
	inALNs = result.Single_Copy_Pep_Align()
	catAln(inALNs, outALN)


def cafe_count(OFdir, outCount):
	'''Generate the gene family count file for CAFE.'''
	def _count(genes):
		return len([v for v in genes.split(', ') if v.strip()])
	result = OrthoFinder(OFdir)
	ogfile = result.Orthogroups
	i = 0
	for line in open(ogfile):
		i += 1
		temp = line.rstrip('\r\n').split('\t')
		if i == 1:
			species = temp[1:]
			line = ['Desc', 'Family ID'] + species
			print('\t'.join(line), file=outCount)
			continue
		ogId = temp[0]
		group = temp[1:]
		count = [_count(genes) for genes in group]
		line = ['(null)', ogId] + count
		line = list(map(str, line))
		print('\t'.join(line), file=outCount)


def to_cafe(OFdir, outCount, species=None):
	'''Retrieve the gene counts for the target species set, for CAFE.'''
	species = parse_species(species)
	result = OrthoFinder(OFdir)
	i = 0
	d_reason = {}
	for group in result.orthogroups:
		i += 1
		if i == 1:
			if species is None:
				species = group.species
			line = ['Desc', 'Family ID'] + species
			print('\t'.join(line), file=outCount)
		counter = group.counter
		count = [counter[sp] for sp in species]
		if not cafe_filter(count, d_reason):
			continue
		line = ['(null)', group.ogid] + count
		line = list(map(str, line))
		print('\t'.join(line), file=outCount)
	print('total {total} families, {left} left, removed: {too_high_std} std>100,\
{too_high_max} max>100, {too_many_zero} none_zero<5'.format(
		total=i, **d_reason), file=sys.stderr)


def upsetplot(OFdir, outprefix, species=None, min_count=1):
	species = parse_species(species)
	print(species)
	result = OrthoFinder(OFdir)
	outCountfile = outprefix + '.set'
	outStatsfile = outprefix + '.stats'
	outCount = open(outCountfile, 'w')
	outStats = open(outStatsfile, 'w')
	i = 0
	lines = []
	d_stats = OrderedDict()
	for group in result.orthogroups:
		i += 1
		if i == 1:
			if species is None:
				species = group.species
			line = species
			print('\t'.join(line), file=outCount)
			print('\t'.join(line+['OG_number', 'gene_number']), file=outStats)
			continue
		counter = group.counter
		count = [counter[sp] for sp in species]
		gene_sum = sum(count)
		line = [min(1, v) for v in count]
		line = tuple(line)
		if sum(line) == 0:
			continue
		try:
			d_stats[line][0] += 1
			d_stats[line][1] += gene_sum
		except KeyError:
			d_stats[line] = [1, gene_sum]
		lines += [line]

	for line, count in list(d_stats.items()):
		line = list(line) + count
		line = list(map(str, line))
		print('\t'.join(line), file=outStats)
	outStats.close()

	lines_counter = Counter(lines)
	low_lines = {line for line, count in list(
		lines_counter.items()) if count < min_count}
	for line in lines:
		if line in low_lines:
			continue
		line = list(map(str, line))
		print('\t'.join(line), file=outCount)
	outCount.close()

	sets = ', '.join(['"{}"'.format(sp) for sp in species])
	nsets = len(species)
	nintersects = 0
	for i in range(len(species)):
		nintersects += len(list(itertools.combinations(species, i+1)))

	outfig = outprefix + '.pdf'
	rfile = outprefix + '.r'
	# plot
	rsrc = r'''
setfile = '{}'
outfig = '{}'
data = read.table(setfile, head=T, sep='\t')
library(UpSetR)
pdf(outfig, width=10, height=7)
upset(data, nsets={}, nintersects={},
dev.off()
'''.format(outCountfile, outfig, nsets, nintersects)
	with open(rfile, 'w') as f:
		print(rsrc, file=f)
	cmd = 'Rscript {}'.format(rfile)
	os.system(cmd)


def cafe_filter(fam_size, d_reason={}, max_sd=100, none_zero=5, min_sp=0.5):
	'''Filter out OGs with too large SD, maximum value, or too few non-null values.'''
	import numpy as np
	for reason in ['left', 'too_high_std', 'too_high_max', 'too_many_zero']:
		if reason not in d_reason:
			d_reason[reason] = 0
	clade_count = sum(1 for sp_count in fam_size if sp_count >= 1)
	if clade_count >= none_zero or 1.0*clade_count/len(fam_size) >= min_sp:
		std = np.std(fam_size)
		if std <= max_sd and max(fam_size) <= max_sd:
			d_reason['left'] += 1
			return True
		else:
			if std >= max_sd:
				d_reason['too_high_std'] += 1
			if max(fam_size) >= max_sd:
				d_reason['too_high_max'] += 1
	else:
		d_reason['too_many_zero'] += 1
	return False


def species_specific_genes(OFdir, sp, outTsv, ex_sps=[]):
	'''Species-specific genes'''
	result = OrthoFinder(OFdir)
	print('\t'.join(['gene', 'group']), file=outTsv)
	d_genes = result.get_species_specific2(sp, ex_sps)
	for gene, group in sorted(d_genes.items()):
		print('\t'.join([gene.split('|', 1)[1], str(group)]), file=outTsv)


def bootstrap_species_tree(OFdir, outdir, bootstrap=1000, iqtree_options='-mset JTT'):
	'''Reconstruct  tree using IQ-TREE.'''
	result = OrthoFinder(OFdir)
	msa = result.SpeciesTreeAlignment
	treefile = result.SpeciesTree_rooted
	new_msa = '{}/singlecopy_aligned.faa'.format(outdir)
	cmd = '[ -s {0} ] && trimal -in {0} -automated1 > {1}'.format(msa, new_msa)
	run_cmd(cmd)
	prefix = new_msa
	try:
		root = result.get_root(result.SpeciesTree_rooted)
	except IOError:
		print('re-root with {', file=sys.stderr)
		root = ''
	if root and not re.compile(r'N\d').match(root):
		iqtree_options += ' -o {}'.format(root)
	cmd = "iqtree -s {} -pre {} -bb {} -bnni -nt AUTO  {} > /dev/null".format(
		new_msa, prefix, bootstrap, iqtree_options)
	print('running cmd: {}'.format(cmd), file=sys.stderr)
	os.system(cmd)
	new_treefile = '{}.treefile'.format(prefix)
	new_treefile_rooted = '{}.rooted.tre'.format(prefix)
	print('re-root with {}'.format([root]), file=sys.stderr)
	cmd = 'nw_reroot {} {} > {}'.format(
		new_treefile, root, new_treefile_rooted)
	os.system(cmd)
	return new_treefile_rooted


def singlecopy_tree(OFdir, outdir, bootstrap=1000, iqtree_options='-mset JTT'):
	'''Construct tree using single-copy genes.'''
	result = OrthoFinder(OFdir)
	new_msa = '{}/singlecopy_aligned.faa'.format(outdir)
	prefix = new_msa
	root = result.get_root(result.SpeciesTree_rooted)
	with open(new_msa, 'w') as outALN:
		single_copy_pep_align(OFdir, outALN)
	if root or not re.compile(r'N\d').match(root):
		iqtree_options += ' -o {}'.format(root)
	cmd = "iqtree -s {} -pre {} -bb {} -bnni -nt AUTO  {} > /dev/null".format(
		new_msa, prefix, bootstrap, iqtree_options)
	os.system(cmd)
	new_treefile = '{}.treefile'.format(prefix)
	new_treefile_rooted = '{}.rooted.tre'.format(prefix)
	cmd = 'nw_reroot {} {} > {}'.format(
		new_treefile, root, new_treefile_rooted)
	os.system(cmd)
	return new_treefile_rooted


def get_singlecopy_orthologs(OFdir, outHomo, **kargs):
	result = OrthoFinder(OFdir)
	for genes in result.get_singlecopy_orthologs(**kargs):
		print('\t'.join(genes), file=outHomo)


def get_orthologs(OFdir, outHomo, **kargs):
	result = OrthoFinder(OFdir)
	for g1, g2 in result.get_orthologs(**kargs):
		line = [g1, g2]
		print('\t'.join(line), file=outHomo)
	return
	orthoFiles = result.Orthologues
	idx = 0
	for orthoFile in orthoFiles:
		idx += 1
		i = 0
		for line in open(orthoFile):
			i += 1
			if i == 1:
				continue
			temp = line.strip().split('\t')
			Orthogroup, Genes_1, Genes_2 = temp
			Genes_1 = Genes_1.split(', ')
			Genes_2 = Genes_2.split(', ')
			Genes_1 = list(map(gene_format_o, Genes_1))
			Genes_2 = list(map(gene_format_o, Genes_2))
			info = '%s.%s.%s' % (Orthogroup, idx, i-1)
			for (sp1, g1), (sp2, g2) in itertools.product(Genes_1, Genes_2):
				if sp1 != sp2:
					line = [g1, g2, info]
					print('\t'.join(line), file=outHomo)


def get_oglogs(OFdir, outHomo):
	result = OrthoFinder(OFdir)
	for g1, g2 in result.get_oglogs():
		line = [g1, g2]
		print('\t'.join(line), file=outHomo)


def gene_format_o(gene):
	try:
		sp, g = gene.split('|', 1)
	except ValueError:
		sp, g = None, gene
	return (sp, gene)


def get_paralogs(OFdir, outHomo, min_support=0.5):
	result = OrthoFinder(OFdir)
	Duplications = result.Duplications
	i = 0
	for line in open(Duplications):
		i += 1
		if i == 1:
			continue
		temp = line.strip().split('\t')
		Orthogroup, Species_Trer_Node, Gene_Tree_Node, Support, Type, Genes_1, Genes_2 \
			= temp
		Support = float(Support)
		if Support < min_support:
			continue
		Genes_1 = Genes_1.split(', ')
		Genes_2 = Genes_2.split(', ')
		Genes_1 = list(map(gene_format_p, Genes_1))
		Genes_2 = list(map(gene_format_p, Genes_2))
		info = [Orthogroup, Species_Trer_Node, Gene_Tree_Node, Support, Type]
		info = list(map(str, info))
		info = '_'.join(info)
		for (sp1, g1), (sp2, g2) in itertools.product(Genes_1, Genes_2):
			if sp1 == sp2:
				line = [g1, g2, info]
				print('\t'.join(line), file=outHomo)


def gene_format_p(gene):
	sp, g = gene.split('|', 1)
	sp = sp[:int(len(sp)/2)]
	g = sp + '|' + g
	return (sp, g)


def gene_format_common(gene):
	if not '|' in gene:
		sp, g = None, gene
		return (sp, g)
	sp, g = gene.split('|', 1)
	sp1 = sp[:len(sp)//2]
	sp2 = sp[len(sp)//2+1:]
	if sp1 == sp2:
		sp = sp1
		g = sp + '|' + g
	return (sp, g)


def MCScanX_transposed(OFdir, tsp, cspp, spmap, gff, datadir='data', outdir='result'):
	'''MCScanX_transposed pipeline'''
	if not os.path.exists(datadir):
		os.mkdir(datadir)
	d_spmap = spmap2dict(spmap)
	result = OrthoFinder(OFdir)
	t_abr = d_spmap[tsp]
	t_gff = '%s/%s.gff' % (datadir, t_abr)
	checkpoint = t_gff + '.ok'
	if not os.path.exists(checkpoint):
		prepare_gff(gff, t_gff, [tsp], d_spmap)
		os.mknod(checkpoint)
	else:
		print('checkpoint: %s exists, skipping prepare %s' %
			  (checkpoint, t_gff))

	t_blast = '%s/%s.blast' % (datadir, t_abr)
	d_geneIds = {}
	checkpoint = t_blast + '.ok'
	if not os.path.exists(checkpoint):
		d_geneIds = result.SequenceIDdict
		print(list(d_geneIds.items())[:10])
		blast_files = result.get_blast_files(tsp, tsp)
		pairs = result.get_paralogs(tsp)
		outHomo = '%s/%s.homology' % (datadir, t_abr)
		write_homo(pairs, outHomo)
		prepare_blast(tsp, tsp, pairs, d_geneIds, blast_files, t_blast)
		os.mknod(checkpoint)
	else:
		print('checkpoint: %s exists, skipping prepare %s' %
			  (checkpoint, t_blast))
	# prepare c-t.gff, c-t.blast
	for csp in cspp:
		c_abr = d_spmap[csp]
		c_gff = '%s/%s_%s.gff' % (datadir, t_abr, c_abr)
		checkpoint = c_gff + '.ok'
		if not os.path.exists(checkpoint):
			prepare_gff(gff, c_gff, [tsp, csp], d_spmap)
			os.mknod(checkpoint)
		else:
			print('checkpoint: %s exists, skipping prepare %s' %
				  (checkpoint, c_gff))
		c_blast = '%s/%s_%s.blast' % (datadir, t_abr, c_abr)
		checkpoint = c_blast + '.ok'
		if not os.path.exists(checkpoint):
			if not d_geneIds:
				d_geneIds = result.SequenceIDdict
			blast_files = result.get_blast_files(tsp, csp)
			pairs = result.get_orthologs(tsp, csp)
			outHomo = '%s/%s_%s.homology' % (datadir, t_abr, c_abr)
			write_homo(pairs, outHomo)
			prepare_blast(tsp, csp, pairs, d_geneIds, blast_files, c_blast)
			os.mknod(checkpoint)
		else:
			print('checkpoint: %s exists, skipping prepare %s' %
				  (checkpoint, c_blast))
	# run
	c_abrs = ','.join([d_spmap[csp] for csp in cspp])
	suffix = '_'.join([d_spmap[sp] for sp in [tsp]+cspp])
	outdir += '.' + suffix
	log = 'run_%s.log' % (suffix,)
	cmd = 'MCScanX_h-transposed.pl -i %s -t %s -c %s -o %s -x %s &> %s' % (
		datadir, t_abr, c_abrs, outdir, len(cspp), log)
	print('CMD: %s' % cmd)
	os.system(cmd)
	outcount = 'run_%s.count.xls' % (suffix,)
	count_mcscan(t_abr, outdir, outcount)


def count_mcscan(t_abr, outdir, outcount):
	def _count(inFile):
		return sum([1 for line in open(inFile)]) - 1

	def _out_pairs(inFile, Type, outf):
		i = 0
		for line in open(inFile):
			i += 1
			if i == 1:
				continue
			temp = line.split()
			g1, g2 = temp[0], temp[2]
			line = [g1, g2, Type]
			print('\t'.join(line), file=outf)
	pairType = outdir + '.pair.class'
	f1 = open(pairType, 'w')
	f = open(outcount, 'w')
	line = ['mode', 'genes', 'pairs']
	print('\t'.join(line), file=f)
	for type in ['proximal', 'segmental', 'tandem', 'transposed']:
		genes = '%s/%s.%s.genes' % (outdir, t_abr, type)
		pairs = '%s/%s.%s.pairs' % (outdir, t_abr, type)
		_out_pairs(pairs, type, f1)
		gene_num = _count(genes)
		pair_num = _count(pairs)
		line = [type, gene_num, pair_num]
		line = list(map(str, line))
		print('\t'.join(line), file=f)
	for enope in ['after', 'between']:
		genesx, pairsx = glob.glob('{}/{}.transposed_{}_*.genes'.format(outdir, t_abr, enope)), \
			glob.glob('{}/{}.transposed_{}_*.pairs'.format(outdir, t_abr, enope))
		for genes, pairs in zip(genesx, pairsx):
			type1 = genes.split('/')[-1].split('.')[-2]
			type2 = pairs.split('/')[-1].split('.')[-2]
			assert type1 == type2
			_out_pairs(pairs, type1, f1)
			gene_num = _count(genes)
			pair_num = _count(pairs)
			line = [type1, gene_num, pair_num]
			line = list(map(str, line))
			print('\t'.join(line), file=f)
	f1.close()
	f.close()


def write_homo(pairs, outHomo):
	f = open(outHomo, 'w')
	for line in pairs:
		print('\t'.join(line), file=f)
	f.close()


def prepare_blast(sp01, sp02, pairs, d_geneIds, blast_files, outblast):
	print('extract blast of {} from {}'.format([sp01, sp02], blast_files))
	sp0x = sorted([sp01, sp02])
	gene_pair_set = pairs

	d_blast = {}
	for blasts in blast_files:
		for line in open(blasts):
			temp = line.strip().split('\t')
			g1, g2 = temp[:2]
			g1, g2 = d_geneIds[g1], d_geneIds[g2]
			if (g1, g2) in gene_pair_set or (g2, g1) in gene_pair_set:
				temp[:2] = g1, g2
				bscore = float(temp[11])
				if (g1, g2) in d_blast and d_blast[(g1, g2)][1] < bscore:
					d_blast[(g1, g2)] = [temp, bscore]
				elif (g1, g2) not in d_blast:
					d_blast[(g1, g2)] = [temp, bscore]
	f = open(outblast, 'w')
	for key, (temp, bscore) in list(d_blast.items()):
		print('\t'.join(temp), file=f)
	f.close()


def prepare_gff(inGff, outGff, spp, d_spmap):
	print('extract gff of {} from {}'.format(spp, inGff))
	spp = set(spp)
	f = open(outGff, 'w')
	for line in open(inGff):
		temp = line.strip().split()
		chr = temp[0]
		gene = temp[1]
		if d_spmap[chr] in spp:
			f.write(line)
	f.close()


def spmap2dict(spmap):
	d = {}
	for line in open(spmap):
		temp = line.strip().split()
		d[temp[0]] = temp[2]
	return d


def get_unique_logs(OFdir, outPrefix=''):
	result = OrthoFinder(OFdir)
	out_orth = '{}orthologs.txt'.format(outPrefix)
	out_para = '{}inparalogs.txt'.format(outPrefix)
	out_para2 = '{}inparalogs2.txt'.format(outPrefix)
	out_para3 = '{}inparalogs3.txt'.format(outPrefix)
	with open(out_orth, 'w') as f:
		for g1, g2 in result.get_orthologs():
			print('\t'.join([g1, g2]), file=f)
	with open(out_para, 'w') as f:
		for g1, g2 in result.get_paralogs():
			print('\t'.join([g1, g2]), file=f)
	with open(out_para2, 'w') as f:
		for g1, g2 in result.get_paralogs2():
			print('\t'.join([g1, g2]), file=f)
	with open(out_para3, 'w') as f:
		for g1, g2 in result.get_paralogs3():
			print('\t'.join([g1, g2]), file=f)

def get_paralogs2(OFdir, sps=None, out_para=sys.stdout):
	result = OrthoFinder(OFdir)
	for g1, g2 in result.get_paralogs2(sps=sps):
		print('\t'.join([g1, g2]), file=out_para)

def aln2beast(inALN, outNEX, partify=True):
	import re
	i = 0
	ntax = sum([1 for rc in SeqIO.parse(inALN, 'fasta')])
	for rc in SeqIO.parse(inALN, 'fasta'):
		i += 1
		if i == 1:
			datatype = guess_seqtype(rc.seq)
			desc = rc.description
			try:
				genes = re.compile(r'genes:(\d+)').search(desc).groups()[0]
			except:
				genes = 0
			sites = len(rc.seq)
			try:
				blocks = re.compile(r'blocks:(.*)').search(desc).groups()[0]
				partitions = re.compile(r'(\d+)').findall(blocks)
				partitions = list(map(int, partitions))
			except:
				partitions = []
			assert len(partitions) == int(genes)
			print('''#NEXUS
begin data;
dimensions ntax={ntax} nchar={nchar};
format datatype={datatype} interleave=no gap=-;
matrix'''.format(ntax=ntax, nchar=sites, datatype=datatype), file=outNEX)
		print('{id}\t{seq}'.format(id=rc.id, seq=rc.seq), file=outNEX)
	print(''';
end;''', file=outNEX)
	if partify and partitions:
		print('begin assumptions;', file=outNEX)
		last_start = 1
		i = 0
		for partition in partitions:
			i += 1
			end = last_start + partition-1
			print('	charset part{part} = {start}-{end};'.format(part=i,
				  start=last_start, end=end), file=outNEX)
			last_start = end + 1
		assert end == int(sites)
		print('end;', file=outNEX)


def guess_seqtype(seq, gap='-'):
	char_count = Counter(seq.upper())
	print(char_count, file=sys.stderr)
	nACTG = sum([char_count.get(char, 0) for char in 'ACTG'])
	nACUG = sum([char_count.get(char, 0) for char in 'ACUG'])
	gap = char_count.get('-', 0)
	if 1e2*nACTG / (len(seq)-gap) > 80:
		return 'dna'
	elif 1e2*nACUG / (len(seq)-gap) > 80:
		return 'rna'
	else:
		return 'protein'


def og2gene(OFdir, oglist, outsv, species=None):
	species = parse_species(species)
	result = OrthoFinder(OFdir)
	oglist = {line.strip().split()[0] for line in open(oglist)}
	line = ['gene', 'OG']
	print('\t'.join(line), file=outsv)
	for group in OrthoFinder(OFdir).orthogroups:
		og = group.ogid
		if not og in oglist:
			continue
		for gene in group.get_group(sps=species):
			gene = gene.split('|', 1)[-1]
			line = [gene, og]
			print('\t'.join(line), file=outsv)


def classify_genes(OFdir, outsv=sys.stdout, species=None):
	OrthoFinder(OFdir).classify_genes(sps=parse_species(species), fout=outsv)


def add_og(OFdir, richfile, outrichfile):
	result = OrthoFinder(OFdir)
	d_genes = {}
	for group in OrthoFinder(OFdir).orthogroups:
		og = group.ogid
		for gene in group.get_group():
			gene = gene.split('|', 1)[-1]
			d_genes[gene] = og

	for i, line in enumerate(open(richfile)):
		temp = line.strip().split('\t')
		if i == 0:
			line = temp + ['OG', 'numOG']
		else:
			genes = temp[-1].split(', ')
			ogs = [d_genes[gene] for gene in genes]
			line = temp + [', '.join(ogs), len(set(ogs))]
		line = list(map(str, line))
		print('\t'.join(line), file=outrichfile)


def get_Blast(OFdir, species, outblast):
	species = parse_species(species)
	result = OrthoFinder(OFdir)
	result.get_Blast(species, fout=outblast)


def parse_key_opts(args):
	d = {}
	pops = []
	for i, arg in enumerate(args):
		kv = arg.split('=', 1)
		if len(kv) != 2:
			continue
		pops += [i]
		key, val = kv
		val = tr_numeric(val)
		d[key] = val
	for i in sorted(pops, reverse=1):
		args.pop(i)
	return d


def tr_numeric(val):
	try:
		return int(val)
	except:
		try:
			return float(val)
		except:
			return val


def main():
	print('CMD: {}'.format(' '.join(sys.argv)), file=sys.stderr)
	subcommand = sys.argv[1]
	kargs = parse_key_opts(sys.argv)
	# Reconstruct the tree with bootstrap support.
	if subcommand == 'reTree':
		OFdir = sys.argv[2]
		outdir = sys.argv[3]
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		bootstrap_species_tree(OFdir=OFdir, outdir=outdir)
	# Construct the tree using only single-copy genes.
	elif subcommand == 'reTree2':
		OFdir = sys.argv[2]
		outdir = sys.argv[3]
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		singlecopy_tree(OFdir=OFdir, outdir=outdir)
	elif subcommand == 'orthologs':
		OFdir = sys.argv[2]
		try:
			sp1, sp2 = sys.argv[3:5]
		except:
			sp1, sp2 = None, None
		outHomo = sys.stdout
		get_orthologs(OFdir, outHomo, sp1=sp1, sp2=sp2)
	elif subcommand == 'paralogs1':
		OFdir = sys.argv[2]
		outHomo = sys.stdout
		get_paralogs(OFdir, outHomo)
	elif subcommand == 'paralogs': # 2
		OFdir = sys.argv[2]
		sps = sys.argv[3:]
		get_paralogs2(OFdir, sps)
	elif subcommand == 'uniqlogs':
		OFdir = sys.argv[2]
		get_unique_logs(OFdir)
	# Retrieve all homologous gene pairs for each OG
	elif subcommand == 'oglogs':
		OFdir = sys.argv[2]
		outHomo = sys.stdout
		get_oglogs(OFdir, outHomo)
	# prepare files for cafe
	elif subcommand == 'to_cafe':
		OFdir = sys.argv[2]
		outTsv = sys.stdout
		try:
			species = sys.argv[3]
		except IndexError:
			species = None
		to_cafe(OFdir, outTsv, species=species)
	elif subcommand == 'upset':
		OFdir = sys.argv[2]
		outprefix = 'upset'
		try:
			species = sys.argv[3]
		except IndexError:
			species = None
		try:
			min_count = int(sys.argv[4])
		except IndexError:
			min_count = 1
		upsetplot(OFdir, outprefix, species=species, min_count=min_count)
	elif subcommand == 'MCScanX_transposed':
		OFdir = sys.argv[2]
		tsp, cspp, spmap, gff = sys.argv[3:7]
		cspp = cspp.split(',')
		MCScanX_transposed(OFdir, tsp, cspp, spmap, gff)
	elif subcommand == 'single_copy_cds':   # single_copy_cds alignments
		OFdir = sys.argv[2]
		cdsSeqs = sys.argv[3]
		outALN = sys.stdout
		single_copy_cds_align(OFdir, cdsSeqs, outALN)
	elif subcommand == 'single_copy_cds2':
		OFdir = sys.argv[2]
		pepSeqs = sys.argv[3]
		cdsSeqs = sys.argv[4]
		outALN = sys.stdout
		try:
			species = sys.argv[5]
		except IndexError:
			species = None
		single_copy_cds_align2(OFdir, pepSeqs, cdsSeqs,
							   outALN, species=species, **kargs)
	elif subcommand == 'single_copy_stats':
		OFdir = sys.argv[2]
		try:
			species = sys.argv[3]
		except IndexError:
			species = None
		single_copy_stats(OFdir, species)
	elif subcommand == 'single_copy_pep':  # single_copy_pep alignments
		OFdir = sys.argv[2]
		outALN = sys.stdout
		single_copy_pep_align(OFdir, outALN)
	elif subcommand == 'catAln':  # cocnat alignments
		inALNs = sys.argv[2:]
		outALN = sys.stdout
		catAln(inALNs, outALN)
	elif subcommand == 'to_paml':
		OFdir = sys.argv[2]
		outDir, cdsSeq = sys.argv[3:5]
		try:
			species = sys.argv[5]
		except IndexError:
			species = None
		singlecopy = False if len(sys.argv) > 6 else True  # default: True
		to_paml(OFdir, outDir, cdsSeq, species=species, singlecopy=singlecopy)
	elif subcommand == 'species_specific':
		OFdir = sys.argv[2]
		sp = sys.argv[3]
		ex_sps = sys.argv[4:]
		outTsv = sys.stdout
		species_specific_genes(OFdir, sp, outTsv, ex_sps)
	# Convert the alignment file to a NEX file for BEAST input.
	elif subcommand == 'aln2beast':
		inALN = sys.argv[2]
		outNEX = sys.stdout
		aln2beast(inALN, outNEX)
	# Prepare single-copy OGs for two species.
	elif subcommand == 'singlecopy_orthologs':
		OFdir = sys.argv[2]
		try:
			sp1, sp2 = sys.argv[3:5]
		except:
			sp1, sp2 = None, None
		outHomo = sys.stdout
		print(OFdir, sp1, sp2, file=sys.stderr)
		get_singlecopy_orthologs(OFdir, outHomo, sp1=sp1, sp2=sp2)
	elif subcommand == 'to_astral':
		OFdir = sys.argv[2]
		pepSeq = sys.argv[3]
		outTrees = sys.stdout
		try:
			species = sys.argv[4]
		except IndexError:
			species = None
		try:
			cutoff = float(sys.argv[5])
		except IndexError:
			cutoff = 0.7
		to_astral(OFdir, pepSeq, outTrees,
				  species=species, min_singlecopy=cutoff)
	elif subcommand == 'to_hybpiper':
		OFdir = sys.argv[2]
		cdsSeq = sys.argv[3]
		outOGSeq = sys.stdout
		try:
			species = sys.argv[4]
		except IndexError:
			species = None
		try:
			cutoff = float(sys.argv[5])
		except IndexError:
			cutoff = 0.7
		to_hybpiper(OFdir, cdsSeq, outOGSeq,
					species=species, min_singlecopy=cutoff)
	elif subcommand == 'singlecopy_stats':
		OFdir = sys.argv[2]
		try:
			species = sys.argv[3]
		except IndexError:
			species = None
		to_hybpiper(OFdir, species=species, only_stats=True)
	elif subcommand == 'venn':
		OFdir = sys.argv[2]
		outTsv = sys.stdout
		species = sys.argv[3]
		venn(OFdir, outTsv, species=species)
	elif subcommand == 'pan':
		OFdir = sys.argv[2]
		try:
			species = sys.argv[3]
		except IndexError:
			species = None
		pan_stats(OFdir, species)
	elif subcommand == 'cn_density':
		OFdir = sys.argv[2]
		copy_number_density(OFdir)
	elif subcommand == 'to_codeml':  # Positive Selection
		OFdir = sys.argv[2]
		outDir, pepSeq, cdsSeq = sys.argv[3:6]
		try:
			species = sys.argv[6]
		except IndexError:
			species = None
		to_codeml(OFdir, outDir, pepSeq, cdsSeq, species=species)
	# List the genes of each OG for enrichment analysis.
	elif subcommand == 'og2gene':
		OFdir = sys.argv[2]
		oglist = sys.argv[3]
		try:
			species = sys.argv[4]
		except IndexError:
			species = None
		outsv = sys.stdout
		og2gene(OFdir, oglist, outsv, species=species)
	# Gene classification: single-copy, multi-copy, etc.
	elif subcommand == 'classify_genes':
		OFdir = sys.argv[2]
		try:
			species = sys.argv[3]
		except IndexError:
			species = None
		outsv = sys.stdout
		classify_genes(OFdir, outsv, species)
	# Obtain all homologous gene pairs for the target species set.
	elif subcommand == 'homo_pairs':
		OFdir = sys.argv[2]
		outpairs = sys.stdout
		try:
			species = sys.argv[3]
		except IndexError:
			species = None
		get_all_homo_pairs(OFdir, outpairs, species=species)
	# Recover false-negative orthologous genes based on synteny.
	elif subcommand == 'retrieve_orthologs':
		OFdir = sys.argv[2]
		collinearity, gff = sys.argv[3:5]
		retrieve_orthologs(OFdir, collinearity, gff, fout=sys.stdout)
	# Recover false-negative orthologous genes in the syntenic flanking regions.
	elif subcommand == 'retrieve_allele':
		OFdir = sys.argv[2]
		collinearity, gff = sys.argv[3:5]
		try:
			sps = sys.argv[5]
		except IndexError:
			sps = None
		retrieve_allele(OFdir, collinearity, gff, fout=sys.stdout, sps=sps)
	elif subcommand == 'enrich_add_og':
		OFdir = sys.argv[2]
		richfile = sys.argv[3]
		outrichfile = sys.stdout
		add_og(OFdir, richfile, outrichfile)
	elif subcommand == 'format_for_iqtree':
		inSeq = sys.argv[2]
		outSeq = sys.stdout
		format_for_iqtree(inSeq, outSeq)
	elif subcommand == 'collinear_og_trees':
		OFdir = sys.argv[2]
		collinearity, taxon, seqfile = sys.argv[3:6]
		collinear_og_trees(OFdir, collinearity, taxon, seqfile)
	elif subcommand == 'get_Blast':
		OFdir = sys.argv[2]
		species = sys.argv[3]
		outblast = sys.stdout
		get_Blast(OFdir, species, outblast)
	elif subcommand == 'count_og':
		OFdir = sys.argv[2]
		species = sys.argv[3]
		count_og(OFdir, species)
	elif subcommand == 'to_wgdi':
		try:
			OFdir = sys.argv[2]
		except IndexError:
			OFdir = 'OrthoFinder/OrthoFinder/Results_*/'
#		try:
#			species = sys.argv[-1]
#		except IndexError:
#			species = None
#		OrthoFinder(OFdir).to_wgdi(parse_species(species), **kargs)
		OrthoFinder(OFdir).to_wgdi(**kargs)
	elif subcommand == 'test_spd':
		test_sonicparanoid()
	elif subcommand == 'mcl2count':
		mcl2count(sys.argv[2])
	else:
		raise ValueError('Unknown command: {}'.format(subcommand))


if __name__ == '__main__':
	main()
