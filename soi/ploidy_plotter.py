# coding: utf-8
import sys
import os
import argparse
from math import sqrt, ceil
import networkx as nx
import numpy as np
import itertools
import matplotlib.pyplot as plt
from .mcscan import XCollinearity, XGff
import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams['pdf.fonttype'] = 42


def add_ploidy_opts(parser):
	parser.add_argument('--window_size', metavar='INT', type=int, default=50,
						help="window_size. [default=%(default)s]")
	parser.add_argument('--window_step', metavar='INT', type=int, default=10,
						help="window_step. [default=%(default)s]")
	parser.add_argument('--min_block', metavar='INT', type=int, default=None,
						help="min genes for a block. [default=%(default)s]")
	parser.add_argument('--max_distance', metavar='INT', type=int, default=20,
						help="max distance from anchor genes. [default=%(default)s]")
	parser.add_argument('--min_overlap', metavar='FLOAT', type=float, default=0.4,
						help="min overlap for covering a reference window. [default=%(default)s]")
	parser.add_argument('--output_depth', metavar='FILE', type=str, default=None,
						help="output depth data to a file. [default=%(default)s]")
	parser.add_argument('--max_ploidy', metavar='INT', type=int, default=10,
                        help="upper limit for x axis. [default=%(default)s]")
	parser.add_argument('--color', metavar='COLOR', type=str, default=None,
						help="bar fill color. [default=%(default)s]")
	parser.add_argument('--edgecolor', metavar='COLOR', type=str, default=None,
						help="bar edge color. [default=%(default)s]")
	parser.add_argument('--as-proportion', action='store_true', default=False,
						help="show y-axis as proportion (0-1) instead of count. [default=%(default)s]")

def ploidy_args(parser):
	# parser.add_argument('-s', '--collinearity', metavar='INPUT_BLOCK_FILE', type=str,
						# required=True, help="the blocks (*.collinearity, output of MCSCANX)")
	# parser.add_argument('-g', '--gff', metavar='INPUT_GENE_GFF_FILE', type=str,
						# required=True, help="the annotation gff file (one of MCSCANX input)")
	parser.add_argument('-s', '-synteny', metavar='FILE', type=str, required=True, nargs='+',
						dest='collinearity',
						help="syntenic block file (*.collinearity, output of MCSCANX/WGDI)[required]")
	parser.add_argument('-g', '-gff', metavar='FILE', type=str, required=True, nargs='+',
						dest='gff',
						help="gene annotation gff file (*.gff, one of MCSCANX/WGDI input)[required]")
	parser.add_argument('-r', '-ref', metavar='SPECIES', type=str, required=True, dest='ref',
						help="reference species [required]")
	parser.add_argument('-q', '-qry', metavar='SPECIES', nargs='+', type=str, dest='qry',
						required=True, help="query species [required]")
	parser.add_argument('-pre', '-prefix', metavar='PREFIX', type=str,
						dest='output', default=None, help="output prefix")
	parser.add_argument('--format', metavar='figure file out format', action='append',
						default=['pdf', 'png'], help="default=%(default)s")
	parser.add_argument('--nrow', metavar='nrow', type=int, default=None,
						help="number of rows. default=%(default)s")
	parser.add_argument('--min_same_block', type=int, default=25,
						help=argparse.SUPPRESS) #"min gene number in a block on the same chromosome. default=%(default)s")
	add_ploidy_opts(parser)

def makeArgparse():
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		)
	ploidy_args(parser)
	args = parser.parse_args()
	return args

class Args:
	def __init__(self):
		pass


def xmain(**kargs):
	args = Args()
	for k, v in kargs.items():
		setattr(args, k, v)
	return main(args)

def main(args):
#	args = makeArgparse()
	sps = [args.ref] + args.qry
	if args.output is None:
		sps = [x[:2] for x in sps]
		args.output = 'dp.' + '-'.join(sps) + '_' + str(args.window_size)
	if args.nrow is None:
		args.nrow = int(ceil(sqrt(len(args.qry))))
	if args.window_step is None:
		args.window_step = args.window_size / 5
	if args.min_overlap is None:
		args.min_overlap = args.window_size / 2.5
	elif args.min_overlap <= 1:
		args.min_overlap = args.min_overlap*args.window_size
	args.ncol = int(ceil(1e0*len(args.qry) / args.nrow))
	args.outfigs = [args.output+'.'+fmt for fmt in args.format]
	# suptitle = 'Reference: ' + args.ref
	# xlabel = 'Relative Ploidy'.format(args.window_size)
	# args.suptitle = '{} ({})'.format(xlabel, suptitle)
	args.titles = args.qry
	print('{} x {} figure'.format(args.nrow, args.ncol), file=sys.stderr)
#	print(args.__dict__)

	plot_fold(**args.__dict__)


def plot_fold(collinearity, gff, ref, qry, **kargs):
	d_ortholog_graph, d_paralog_graph = parse_collinearity(
		collinearity, ref, qry, **kargs)
	d_coord_path, d_coord_graph = parse_gff(gff, [ref]+qry)
	data = []
	for sp in qry:
		d_fold = get_ploidy(d_coord_path[ref], d_coord_graph[ref],
							d_coord_graph[sp], d_ortholog_graph[sp],
						#	d_paralog_graph[sp], 
							**kargs)
		data += [np.array(sorted(d_fold.items()))]
		#print(sp, sorted(d_fold.items()))
	return plot_bars(data, ref=ref, **kargs)


def plot_bars(data, titles, ax=None, outfigs=None, nrow=1, ncol=1, fontsize=10, 
			  suptitle=None, max_ploidy=10, color='white', edgecolor='black',
			  ylabel='Number of windows', xlabel='Synteny depth', 
			  output_depth=None, mode='w', ref=None,
			  as_proportion=False, **kargs):
	if output_depth:
		save_depth_table(data, titles, output_depth=output_depth, 
			mode=mode, max_ploidy=max_ploidy, ref=ref)
	# bars show depth >= 1 only; depth-0 windows reported in the table
	data = [np.array([row for row in arr if row[0] >= 1]) if len(arr) else arr
			for arr in data]
	if as_proportion:
		ylabel = ylabel.replace('Number', 'Proportion', 1) if ylabel.startswith('Number') else 'Proportion of ' + ylabel
	if ax is None:
		if nrow*ncol == 1:
			ax = plt.subplot(111)
			ax = [ax]
		else:
			fig, ax = plt.subplots(
				nrow, ncol, sharex=True, figsize=(10*ncol/2, 8*ncol/2))
			cells = list(itertools.product(
				list(range(nrow)), list(range(ncol))))
			ax = [ax[cell] for cell in cells]
	else:
		ax = [ax]
	tick_label = list(range(0, max_ploidy+1))
	for i, (dat, title, sax) in enumerate(zip(data, titles, ax)):
		try:
			x = dat[:, 0]
			y = dat[:, 1].astype(float)
			if as_proportion and y.sum() > 0:
				y = y / y.sum()
			sax.bar(x, y, align='center', color=color, edgecolor=edgecolor)
		except IndexError:
			pass
		if title is not None:
			sax.set_title(title)
		sax.set_xlim(0, max_ploidy)
		if xlabel is not None and i >= (nrow-1)*ncol:
			sax.set_xlabel(xlabel, fontsize=fontsize)
		if ylabel is not None and i % ncol == 0:
			sax.set_ylabel(ylabel, fontsize=fontsize)
	plt.xticks(tick_label)
	if suptitle is not None:
		plt.suptitle(suptitle)
	if outfigs is not None:
		for outfig in outfigs:
			plt.savefig(outfig)
	else:
		return ax
def save_depth_table(data, titles, ref=None, output_depth=None, mode='w', max_ploidy=10):
    """
    每行一个物种，列为不同深度 (1 to max_ploidy)
    """
    # 1. 构建表头: Species, 0, 1, 2, 3, ..., 10+
    header = ['Reference', "Query"] + [str(i) for i in range(0, max_ploidy)] + [f"{max_ploidy}+"]
    
    rows = [header]
    
    # 2. 填充每个物种的数据
    for i, arr in enumerate(data):
        species = titles[i]
        # 初始化当前物种的计数器
        counts = {p: 0 for p in range(0, max_ploidy + 1)}
        
        for depth, count in arr:
            if depth >= max_ploidy:
                counts[max_ploidy] += count
            else:
                counts[depth] += count
        
        # 构造当前行：物种名 + 各深度的计数
        row = [str(ref), species] + [str(counts[p]) for p in range(0, max_ploidy + 1)]
        rows.append(row)

    # 3. 拼接为 TSV 文本
    output_text = "\n".join(["\t".join(row) for row in rows]) + "\n"

    # 4. 输出
    if output_depth and output_depth != 'stdout':
        with open(output_depth, mode=mode, encoding='utf-8') as f:
            f.write(output_text)
    else:
        print(output_text)

    return rows

def parse_collinearity(collinearity, ref, qry, min_block=10, min_same_block=25, **kargs):
	d_ortholog_graph = {}
	d_paralog_graph = {}
	for sp in qry:
		d_ortholog_graph[sp] = nx.Graph()
		d_paralog_graph[sp] = nx.Graph()

	qry = set(qry)
	for rc in XCollinearity(collinearity):
		if rc.chr1 == rc.chr2 and rc.N < min_same_block:
			continue
		if min_block is not None and rc.N < min_block:
			continue
		sp1, sp2 = rc.species
		if sp1 == sp2 and sp1 in qry:
#			d_paralog_graph[sp1].add_edges_from(rc.pairs)
			d_ortholog_graph[sp1].add_edges_from(rc.pairs)
			continue
		elif (sp1 == ref and sp2 in qry):
			d_ortholog_graph[sp2].add_edges_from(rc.pairs)
		elif (sp2 == ref and sp1 in qry):
			d_ortholog_graph[sp1].add_edges_from(rc.pairs)
	return d_ortholog_graph, d_paralog_graph


def parse_gff(gff, sps):
	d_coord_graph = {}
	for sp in sps:
		d_coord_graph[sp] = nx.Graph()
	sps = set(sps)
	d_gff = {}
	for line in XGff(gff):
		if not line.species in sps:
			continue
		key = (line.species, line.chrom)
		try:
			d_gff[key] += [line]
		except KeyError:
			d_gff[key] = [line]

	d_coord_path = {}
	for (sp, chrom), lines in list(d_gff.items()):
		lines = sorted(lines, key=lambda x: (x.start, -x.end))
		genes = [line.gene for line in lines]
		try:
			d_coord_path[sp] += [genes]
		except KeyError:
			d_coord_path[sp] = [genes]
		for i, line in enumerate(lines):
			d_coord_graph[sp].add_node(line.gene, chrom=chrom, index=i)
		for i, line in enumerate(lines[1:]):
			edge = (lines[i].gene, line.gene)
			d_coord_graph[sp].add_edge(*edge)
	return d_coord_path, d_coord_graph


def get_ploidy(ref_coord_paths, ref_coord_graph, qry_coord_graph, rq_ortholog_graph,
			   window_size=20, window_step=10, **kargs):
	'''For each query, how many segments correspond to the query.'''
	d_fold = {}
	for path in ref_coord_paths:
		for i in range(0, len(path), window_step):
			start, end = i, i+window_size
			if end > len(path):
				end = len(path)
			if end - start < window_size/2:
				continue
			bin = path[start:end]
			orthologs = []
			for gene in bin:
				if rq_ortholog_graph.has_node(gene):
					orthologs += rq_ortholog_graph.neighbors(gene)
			if len(orthologs) < 2:  # ploidy=0
				try:
					d_fold[0] += 1
				except KeyError:
					d_fold[0] = 1
				continue
			qry_clusters = cluster_genes(orthologs, qry_coord_graph, **kargs)
			qry_blocks = list(nx.connected_components(qry_clusters))
			if not qry_blocks:  # not into blocks
				continue
			ref_blocks = map_graph(bin, rq_ortholog_graph, qry_blocks)
			ref_clusters = overlap_blocks(ref_blocks, ref_coord_graph, **kargs)
			ncmpt1 = len(qry_blocks)
#			ncmpt2 = max([len(x)
#						 for x in nx.connected_components(ref_clusters)])
			max_blocks = max(nx.connected_components(ref_clusters), key=lambda x: len(x))
			ncmpt2 = sweep_line(max_blocks, ref_coord_graph)
			ploidy = ncmpt2
			try:
				d_fold[ploidy] += 1
			except KeyError:
				d_fold[ploidy] = 1
	return d_fold

def sweep_line(ref_blocks, ref_coord_graph):
    '''扫描线求最大局部覆盖深度'''
    intervals = []
    for block in ref_blocks:
        idxs = [ref_coord_graph.nodes[g]['index'] for g in block]
        intervals.append((min(idxs), max(idxs)))

    events = []
    for lo, hi in intervals:
        events.append((lo, 1))
        events.append((hi + 1, -1))
    events.sort()

    depth = max_depth = 0
    for _, delta in events:
        depth += delta
        max_depth = max(max_depth, depth)

    return max_depth

def map_graph(bin, rq_ortholog_graph, qry_blocks):
	'''map qry block to ref block'''
	ref_blocks = []
	for block in qry_blocks:
		ref_block = []
		for gene in block:
			ref_block += list(set(rq_ortholog_graph[gene]) & set(bin))
		ref_blocks += [ref_block]
	return ref_blocks


def overlap_blocks(blocks, coord_graph, min_overlap=3, **kargs):
	'''Concatenate blocks that have overlap.'''
	blocks = list(map(tuple, blocks))
	G = nx.Graph()
	for b in blocks:
		G.add_node(b)
	for b1, b2 in itertools.combinations(blocks, 2):
		i1 = [coord_graph.nodes[x]['index'] for x in b1]
		i2 = [coord_graph.nodes[x]['index'] for x in b2]
		min_i1, max_i1 = min(i1), max(i1)
		min_i2, max_i2 = min(i2), max(i2)
		if min(max_i1, max_i2) - max(min_i1, min_i2) + 1 >= min_overlap:  # overlap
			G.add_edge(b1, b2)
	return G


def cluster_genes(genes, coord_graph, max_distance=25, **kargs):
	'''Cluster genes into blocks based on their coordinates.'''
	d_bin = {}
	for gene in genes:
		try:
			chrom = coord_graph.nodes[gene]['chrom']
		except KeyError:
			continue
		try:
			d_bin[chrom] += [gene]
		except KeyError:
			d_bin[chrom] = [gene]
	G = nx.Graph()
	for chrom, genes in list(d_bin.items()):
		genes = sorted(genes, key=lambda x: coord_graph.nodes[x]['index'])
		for i, gene in enumerate(genes[1:]):
			n1, n2 = genes[i], gene
			i1, i2 = coord_graph.nodes[n1]['index'], coord_graph.nodes[n2]['index']
			if i2 - i1 < max_distance:
				G.add_edge(n1, n2)
	return G


if __name__ == '__main__':
	main()
