#coding:utf-8
import sys, os
import re
import copy
import random
import numpy as np
import networkx as nx
from collections import Counter, OrderedDict
import itertools
from Bio import SeqIO, Phylo
from lazy_property import LazyWritableProperty as lazyproperty

from .OrthoFinder import catAln, format_id_for_iqtree, OrthoMCLGroup, OrthoMCLGroupRecord, OrthoFinder, parse_species, exists_and_size_gt_zero
from .small_tools import mkdirs, flatten, test_s, test_f, parse_kargs, rmdirs, lazy_decode
from .RunCmdsMP import run_cmd, run_job, logger
from .small_tools import open_file as open
#from creat_ctl import sort_version

def get_chrom(chrom):
	try:
		return chrom.split('|', 1)[1]
	except IndexError:
		return chrom

class Gene():
	def __init__(self, info):
		self.info = info
		(id, chr, start, end, strand) = self.info
		self.id = self.gene = id
		self.chr = self.chrom = chr
		self.start = start
		self.end = end
		self.strand = strand
		self.coord = (chr, start, end)
		try: self.species, self.raw_gene = id.split('|', 1)
		except ValueError: 
			logger.warn('Species id not found in `{}`, ignoring'.fomart(id))
			self.species = self.raw_gene = id
	def __str__(self):
		return self.id
	def __hash__(self):
		return hash(self.id)
	def __eq__(self, other):
		return self.id == other.id
	def to_bed(self):
		line = [self.chr, self.start, self.end, self.id, 0, self.strand]
		return line
	@property
	def ichr(self):
		return int(re.compile(r'[^\d\s]+(\d+)').match(self.chr).groups()[0])

class KaKs():
	def __init__(self, info, fdtv=False, yn00=False, wgdi=False, method='NG86', **kargs):
		self.info = info
		if yn00:
			self.parse_yn00(method=method)
		elif wgdi:
			self.parse_wgdi(method=method)
		elif fdtv:
			self.parse_4dtv()
		else:
			self.parse_ks()
		self.parse_pair()	
	def parse_4dtv(self):
		(Sequence, fD_Sites, Identical_Sites, TS_Sites, TV_Sites, fDS, fDTS, fDTV, Corrected_4DTV) = self.info
		self.sequence = Sequence
		try: self.ks = float(Corrected_4DTV)
		except ValueError: self.ks = 0
	def parse_yn00(self, method='NG86'):
		(Sequence, dS_YN00, dN_YN00, dS_NG86, dN_NG86, 
			dS_LWL85, dN_LWL85, dS_LWL85m, dN_LWL85m, dS_LPB93, dN_LPB93) = self.info
		self.sequence = Sequence
		d = {'YN00': dS_YN00, 'NG86':dS_NG86, 'LWL85':dS_LWL85, 'LWL85m':dS_LWL85m, 
			 'LPB93':dS_LPB93}
		ks = d[method.upper()]
		try: self.ks = float(ks)
		except ValueError: 
			self.ks = 0
			return
		if self.ks < 0 or np.isnan(self.ks):
			self.ks = 0
	def parse_wgdi(self, method='NG86'):
		try: (g1, g2, dN_NG86, dS_NG86, dN_YN00, dS_YN00) = self.info
		except ValueError: 
			try: self.pair = (g1, g2) = tuple(self.info)
			except  ValueError:
				self.pair = (None, None)
			self.ks = None
			return
		self.pair = (g1, g2)
		d = {'YN00': dS_YN00, 'NG86':dS_NG86}
		ks = d[method.upper()]
		try: self.ks = float(ks)
		except ValueError:
			self.ks = 0
		if self.ks < 0 or np.isnan(self.ks):	# -1
			self.ks = None
	def parse_ks(self):
		(Sequence, Method, Ka, Ks, Ka_Ks, P_Value, Length, 
			S_Sites, N_Sites, Fold_Sites, Substitutions, 
			S_Substitutions, N_Substitutions, Fold_S_Substitutions,
			Fold_N_Substitutions, Divergence_Time, Substitution_Rate_Ratio,
			GC, ML_Score, AICc, Akaike_Weight, Model) = self.info
		self.sequence = Sequence
		self.method = Method
#		try: self.ka = float(Ka)
#		except ValueError: self.ka =None
		try: self.ks = float(Ks)
		except ValueError: self.ks = 0	# too few substitution to calculate, NA
#		try: self.kaks = float(Ka_Ks)
#		except ValueError: self.kaks = None
	def parse_pair(self):
		if not hasattr(self, 'pair'):
			self.pair = re.compile(r'(\S+)\-([A-Z][a-z]*[_\-]\S+\|\S+)').match(self.sequence).groups() #tuple(Sequence.split('-'))
		self.species = self.gene2species(self.pair)
	def gene2species(self, gene_pair):
		sp1, sp2 = gene_pair if gene_pair==(None, None) else [x.split('|')[0] for x in gene_pair]
		return SpeciesPair(sp1, sp2)
	def write(self, fout):
		print('\t'.join(self.info), file=fout)
class KaKsParser:
	def __init__(self, kaks, **kargs):
		self.kaks = kaks
		self.kargs = kargs
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.kaks):
			line = lazy_decode(line)
			temp = line.rstrip().split()
#			print >> sys.stderr, temp
			#print(temp)
			if temp[0] in {'Sequence', 'id1'}:
				if temp[1] == 'dS-YN00':
					self.kargs['yn00'] = True
				elif temp[2] == 'ka_NG86':
					self.kargs['wgdi'] = True
				elif temp[1] == '4D_Sites':
					self.kargs['fdtv'] = True
				continue
#			print >> sys.stderr, self.kargs
			kaks = KaKs(temp, **self.kargs)
			yield kaks
	def to_dict(self):
		d = {}
		for kaks in self:
			ks = kaks.ks
			d[kaks.pair] = ks
			d[tuple(reversed(kaks.pair))] = ks
		return d

def collinearity_ratio(collinearity, chrmap, outMat, min_N=20):
	from creat_ctl import get_good_chrs, Chrs
	d_chrs = {rc.chr:rc.geneN for rc in Chrs(chrmap)}
	good_chrs = get_good_chrs(chrmap, min_genes=200)
	good_chrs = set(good_chrs)
	chrs = set([])
	d_count = {}
	for rc in Collinearity(collinearity, ):
		if rc.N < min_N:
			continue
		chr1, chr2 = rc.chr1, rc.chr2
#		if chr1 == chr2:
#			continue
		if set([chr1, chr2]) - good_chrs:
			continue
		key = tuple(sorted([chr1, chr2]))
		try: d_count[key] += rc.N
		except KeyError: d_count[key] = rc.N
		chrs = chrs | set([chr1, chr2 ])
#	print >> sys.stderr, d_count
	chrs = sorted(chrs)
	print('\t'.join(['']+chrs), file=outMat)
	for chr1 in chrs:
		n1 = d_chrs[chr1]
		line = [chr1]
		for chr2 in chrs:
			n2 = d_chrs[chr2]
			key = tuple(sorted([chr1, chr2]))
			cn = d_count.get(key, 0)
			
			ratio = cn*2.0 / (n1+n2)
#			print  key, cn, n1,n2, ratio
			line += [str(ratio)]
			
		print('\t'.join(line), file=outMat)
class XCollinearity:
	def __init__(self, collinearities, orthologs=None, homo_class=None, **kargs):
		self.orthologs = orthologs
		self.homo_class = homo_class
		self.kargs = kargs
		self.collinearities = self._parse_list(collinearities)
	def __iter__(self):
		return self._parse()
	def _parse_list(self, _collinearities):
		collinearities = []
		if isinstance(_collinearities, str):
			return [_collinearities]
		for collinearity in _collinearities:
			if lazy_decode(open(collinearity).read(1)) == '#': # or self.kargs.get('homology'):
				collinearities += [collinearity]
			else:
				files, _unknown = [], []
				i = 0
				for line in open(collinearity):
					_file = line.strip().split()[0]
					if _file:
						i += 1
					if test_s(_file):
						files += [_file]
					else:
						_unknown += [_file]
				if len(files) == i:
					collinearities += files
				elif len(files) == 0:
					collinearities += [collinearity]
				else:
					logger.warn('Files empty or not exists: {}'.format(_unknown))
					collinearities += files
#				collinearities += [line.strip().split()[0] for line in open(collinearity)]
		return collinearities
	def _parse(self):
		if self.orthologs is not None:
			ortholog_pairs = set(XOrthology(self.orthologs, **self.kargs))
			logger.info('\t{} homologous pairs'.format(len(ortholog_pairs)))
		logger.info('parsing {} collinearity files: {}...'.format(
			len(self.collinearities), self.collinearities[:3]))
		nblock, ngene = 0, 0
		for collinearity in self.collinearities:
	#		logger.info('parsing {}'.format(collinearity))
			for rc in Collinearity(collinearity, **self.kargs):
				nblock += 1
				ngene += rc.N
				if self.orthologs is not None:
					pairs = {Pair(*x) for x in rc.pairs}
					intersect = pairs & ortholog_pairs
					ratio = 1.0*len(intersect) / len(pairs)
					rc.on = len(intersect) # syntenic orthologs
					rc.oi = ratio
					if self.homo_class is not None:
						rc.intersect = intersect
						rc.substract = pairs - ortholog_pairs
				if self.orthologs is not None:
					rc.ton = len(ortholog_pairs)	# all syntenic orthologs
					rc.ortholog_pairs = ortholog_pairs
				yield rc
		logger.info('\t{} collinearity blocks, {} collinearity genes'.format(nblock, ngene))
class XOrthology:
	def __init__(self, orthologs, **kargs):
		self.orthologs = self._parse_list(orthologs)
		self.kargs = kargs
	def __iter__(self):
		return self._parse()
	def _parse_list(self, _orthologs):
#		orthologs = []
		if isinstance(_orthologs, str):
			return [_orthologs]
		else:
			return _orthologs
	def _parse(self):
		for ortholog in self.orthologs:
			logger.info('parsing {}...'.format(ortholog))
			if os.path.isdir(ortholog):	# orthofinder
				for pair in OrthoFinder(ortholog).get_homologs(**self.kargs):
					yield Pair(*pair)
			else:	# orthomcl or similar format
				for rc in Pairs(ortholog, parser=Pair):
#					print(rc, file=sys.stderr)
					yield rc
def identify_orthologous_blocks(collinearities=None, orthologs=None, fout=sys.stdout, 
		gff=None, kaks=None, source=None, min_n=0, min_dist=None, #paralog=False, both=True 
		min_ratio=0.5, max_ratio=1, species=None, homo_class=None, out_stats=None, test_diff=False):
	if species is not None:
		species = parse_species(species)
	if homo_class is not None:
		out_class = open(homo_class, 'w')
	if out_stats is not None:
		out_stats = open(out_stats, 'w')
	if test_diff:
		d_ks = {}
	pre_nb, pre_ng, post_nb, post_ng = 0, 0,0,0
	post_no = 0
	total_oi = 0
	d_sp_count = OrderedDict()
	logger.info('filtering collinearity...')
	rn, rd, ro = 0,0,0
	for rc in XCollinearity(collinearities, orthologs=orthologs, 
				gff=gff, kaks=kaks, source=source, sps=species, homo_class=homo_class):
		pre_nb += 1
		pre_ng += rc.N
		sp_pair = rc.species
		if sp_pair not in d_sp_count:
			d_sp_count[sp_pair]= Count()
		d_sp_count[sp_pair].pre_nb += 1
		d_sp_count[sp_pair].pre_ng += rc.N
		if rc.N < min_n:
			rn += 1
			continue
		if min_dist is not None and rc.is_tandem(min_dist):
			rd += 1
			continue
		if not (min_ratio < rc.oi <= max_ratio):
			ro += 1
			continue
		post_nb += 1
		post_ng += rc.N
		post_no += rc.on # syntenic orthologs
		total_oi += rc.oi * rc.N
		d_sp_count[sp_pair].post_nb += 1
		d_sp_count[sp_pair].post_ng += rc.N
		d_sp_count[sp_pair].post_no += rc.on

		rc.write(fout)
		if homo_class:
			substract = pairs - ortholog_pairs
			for pairs, cls in zip((rc.intersect, rc.substract), ('ortholog', 'non-ortholog')):
				for pair in pairs:
					line = list(pair) + [cls]
					print('\t'.join(line), file=out_class)
	for pair in rc.ortholog_pairs:
		sp_pair = pair.species
		if sp_pair not in d_sp_count:
			continue
		d_sp_count[sp_pair].pre_no += 1				
	logger.info('Synteny: Pre-filter: {} blocks, {} pairs; Post-filter: {} ({:.1%}) blocks, {} ({:.1%}) pairs.'.format(
		pre_nb, pre_ng, post_nb, 1.0*post_nb/pre_nb, post_ng, 1.0*post_ng/pre_ng))
	logger.info('Orthology: Pre-filter: {} pairs; Post-filter: {} ({:.1%}) pairs.'.format(
		rc.ton, post_no, 1.0*post_no/rc.ton))
	logger.info('Post-filter mean OrthoIndex: {:.2f}'.format(total_oi/post_ng))
	if homo_class is not None:
		out_class.close()
	if out_stats is  None:
		return
	logger.info('Output stats..')
	line = ['Species1', 'Species2', 'Pre-filter number of orthologous gene pairs', 'Post-filter number of orthologous gene pairs',
		'Pre-filter number of syntenic blocks', 'Post-filter number of syntenic blocks', 
		'Pre-filter number of syntenic gene pairs', 'Post-filter number of syntenic gene pairs',]
	print('\t'.join(line), file=out_stats)
	for sp_pair, self in d_sp_count.items():
		line = [sp_pair[0], sp_pair[1], self.pre_no, self.post_no, self.pre_nb, self.post_nb, self.pre_ng, self.post_ng]
		print('\t'.join(map(str, line)), file=out_stats)
	out_stats.close()

class Count:
	def __init__(self):
		self.pre_nb, self.pre_ng, self.post_nb, self.post_ng = 0,0,0,0
		self.pre_no, self.post_no = 0,0
class Collinearity():
	'''
	blocks = Collinearity(blockfile)
	for rc in blocks:
		genes1,genes2 = rc.genes1, rc.genes2
	'''
	def __init__(self, collinearity=None, gff=None, chrmap=None, kaks=None, homology=False, source=None, **ks_args):
		self.collinearity = collinearity
		self.gff = gff
		self.chrmap = chrmap
		self.kaks = kaks
		self.d_kaks = self.parse_kaks(**ks_args)
		#if self.gff is not None:
		self.d_gene = self.parse_gff()
		#if self.chrmap is not None:
		self.d_chr = self.map_chr()
		self.homology = homology
		self.source = source
	def __iter__(self):
		return self.parse()
	def __str__(self):
		return self.block
	def __repr__(self):
		return '<MCSCANX Collinearity parser>'
	def write(self, f, ):
		if self.has_head:
			f.write(self.header)
		f.write(self.block)
	def parse(self):
		if not self.homology:
			lines = []
			head = []
			self.has_head = 1
			for line in open(self.collinearity):
				line = lazy_decode(line)
				if re.compile(r'#+ Alignment').match(line):	# mcscanx or wgdi
					if self.source is None and re.compile(r'# Alignment').match(line):
						self.source = 'wgdi'
						self.has_head = 0
					self.header = ''.join(head)
					if lines:
						self.parse_lines(lines)
						yield self
						lines = []
					lines.append(line)
				elif re.compile(r'###$').match(line):	# jcvi
					self.source = 'jcvi'
					self.has_head = 0
					if lines:
						self.parse_lines(lines)
						yield self
						lines = []
					lines.append(line)
				elif line.startswith('#'):
					head.append(line)
				else:
					lines.append(line)
			if lines:
				self.parse_lines(lines)
				yield self
		else: # homology
			for line in open(self.collinearity):
				line = lazy_decode(line)
				self.parse_homology_line(line)
				yield self
	
	def parse_lines(self, lines):
		self.block = ''.join(lines)
		genes1, genes2 = [], []
		for i, line in enumerate(lines):
			if i == 0 and self.source == 'jcvi':
				pass
			elif i == 0:
				pattern = r'#+ Alignment (\d+): score=(\S+) \S+value=(\S+) N=(\S+) (\S+)&(\S+) (plus|minus)'
				try: self.Alignment, self.score, self.e_value, self.N, self.chr1, self.chr2, self.orient = \
								re.compile(pattern).match(line).groups()
				except AttributeError:
					print('unparsed LINE: {}'.format(line), file=sys.stderr)
					raise AttributeError()
				self.chrs = (self.chr1, self.chr2)
				self.sp1 = self.short_sp1 = self.chr1[:2]
				self.sp2 = self.short_sp2 = self.chr2[:2]
				self.score = float(self.score)
				self.e_value = float(self.e_value)
				self.length = self.N = int(self.N)
				self.strand = {'plus': '+', 'minus': '-'}[self.orient]
				self.id = self.Alignment
			else:
				if self.source == 'jcvi':
					gene1, gene2, score = line.strip().split()
				elif self.source == 'wgdi':
#					print >> sys.stderr, line.strip().split('\t')
					gene1, idx1, gene2, idx2, strand = line.strip().split()
					
				else: # mcscan
					pattern = r'.*?\d+.*?\d+:\s+(\S+)\s+(\S+)\s+\d+'
					#print line
					try: gene1, gene2 = re.compile(pattern).match(line).groups()
					except AttributeError:
						print('unparsed LINE: {}'.format(line), file=sys.stderr)
						continue
					tmp = line.strip().split()
					if len(tmp) > 5:
						self.ks = float(tmp[-1])
				genes1.append(gene1)
				genes2.append(gene2)
		if self.source == 'jcvi':
			self.N = len(genes1)
			self.Alignment, self.score, self.e_value = 0,0,0
		self.parse_species(gene1, gene2)
		self.parse_genes(genes1,genes2)
	@property
	def mean_score(self):
		try:return self.score / self.N
		except: return 1
	@property
	def info(self):
		return [self.id, self.species1, self.species2, self.chr1, self.chr2, 
			self.istart1,self.iend1,self.istart2,self.iend2, 
			self.N, self.median_ks, self.mean_ks]
	def is_sp_pair(self, sp1, sp2):
		if (sp1, sp2) == (self.short_sp1, self.short_sp2):
			return (sp1, sp2)
		if (sp2, sp1) == (self.short_sp1, self.short_sp2):
			return (sp2, sp1)
		if (sp1, sp2) == self.species.pair:
			return (sp1, sp2)
		if (sp2, sp1) == self.species.pair:
			return (sp2, sp1)
		return False
	def parse_homology_line(self, line):
		temp = line.strip().split()
		gene1, gene2 = temp[:2]
		genes1 = [gene1]
		genes2 = [gene2]
		self.N = 1
		self.parse_species(gene1, gene2)
		self.parse_genes(genes1,genes2)
	def get_species(self):
		species = set([])
		for rc in self:
			species = species | {rc.species1, rc.species2}
		return species
	@property
	def gene_pairs(self):
		return [tuple(map(self.gene2geneid, pair)) for pair in self.pairs]
	@property
	def gene_genes(self):
		return [list(map(self.gene2geneid, genes)) for genes in self.genes]
	def gene2geneid(self, gene):
		return gene.split('|')[1]
	def parse_genes(self, genes1,genes2):
		self.pairs = list(zip(genes1,genes2))
		self.ks = []
		for pair in self.pairs:
			try: ks = self.d_kaks[pair].ks
			except KeyError: ks = None
			self.ks.append(ks)
		self.genes = [genes1,genes2]
		try:	# Gene obj
			self.genes1 = [self.d_gene[x] for x in genes1]
			self.genes2 = [self.d_gene[x] for x in genes2]
		except KeyError:	# string obj
			self.genes1, self.genes2 = genes1, genes2
		self.segment1, self.segment2 = Segment(self.genes1), Segment(self.genes2)
		self.head, self.tail = (genes1[0], genes2[0]), (genes1[-1], genes2[-1])
		self.head1, self.head2 = self.head
		self.tail1,	self.tail2 = self.tail
		try:
			chr10, start10, end10 = self.d_gene[self.head1].coord
			chr11, start11, end11 = self.d_gene[self.tail1].coord
			self.chr1 = chr10
			self.start1 = min(start10, end10, start11, end11)
			self.end1 = max(start10, end10, start11, end11)
			self.length1 = self.end1 - self.start1 + 1
			idx10 = self.d_gene[self.head1].index
			idx11 = self.d_gene[self.tail1].index
			self.istart1 = min(idx10, idx11)
			self.iend1   = max(idx10, idx11)
			try:	# raw chr id from `chr.list`
				self.chr1 = self.d_chr[chr10]
			except KeyError:
				pass
		except KeyError:
			self.start1, self.end1, self.length1 = None, None, None
			self.istart1, self.iend1 = None, None
		try:
			chr20, start20, end20 = self.d_gene[self.head2].coord
			chr21, start21, end21 = self.d_gene[self.tail2].coord
			self.chr2 = chr20
			self.start2 = min(start20, end20, start21, end21)
			self.end2 = max(start20, end20, start21, end21)
			self.length2 = self.end2 - self.start2 + 1
			idx20 = self.d_gene[self.head2].index
			idx21 = self.d_gene[self.tail2].index
			self.istart2 = min(idx20, idx21)
			self.iend2   = max(idx20, idx21)
			try:
				self.chr2 = self.d_chr[chr20]
			except KeyError:
				pass
		except KeyError:
			self.start2, self.end2, self.length2 = None, None, None
			self.istart2, self.iend2 = None, None
	@property
	def good_ks(self):
		return [ks for ks in self.ks if ks is not None]
	@property
	def mean_ks(self):
		return np.mean(self.good_ks)
	@property
	def median_ks(self):
		return np.median(self.good_ks)
	def parse_species(self, gene1, gene2):
		self.species1 = gene1.split('|')[0]
		self.species2 = gene2.split('|')[0]
		self.species = SpeciesPair(self.species1, self.species2)
	def is_tandem(self, max_dist=10):
		if self.species1 != self.species2 or self.chr1 != self.chr2:
			return False
		dist = max(self.istart1, self.istart2) - min(self.iend1, self.iend2)
		if dist < max_dist:
			return True
		return False
	def parse_gff(self):
		d = {}
		if self.gff is None:
			return d
		genes = set([])
		d_chr = {}
		d_length = {}
		for line in open(self.gff):
			line = lazy_decode(line)
			temp = line.rstrip().split('\t')
			chr, gene, start, end = temp[:4]
			if gene in genes:	# remove repeat
				continue
			genes.add(gene)
			try: strand = temp[4]
			except IndexError: strand = None
			start, end = list(map(int, [start, end]))
#			d[gene] = Gene((gene, chr, start, end, strand))
			g = Gene((gene, chr, start, end, strand))
			try: d_chr[chr] += [g]
			except KeyError: d_chr[chr] = [g]
			try: d_length[chr] = max(d_length[chr], end)
			except KeyError: d_length[chr] = end
		d_chrom = {}
		d_ngenes = {}
		for chr, genes in list(d_chr.items()):
			genes = sorted(genes, key=lambda x:x.start)
			d_chrom[chr] = genes
			d_ngenes[chr] = len(genes)
			for i, gene in enumerate(genes):
				gene.index = i
				d[gene.id] = gene
#			print >> sys.stderr, chr, gene.id, i
		self.chr_length = d_length
		self.d_chrom = d_chrom
		self.chr_ngenes = d_ngenes
		return d
	def map_chr(self):
		d = {}
		if self.chrmap is None:
			return d
		for line in open(self.chrmap):
			temp = line.rstrip().split()
			chrid, chr, sp = temp[:3]
			d[chrid] = chr
		return d
	def parse_kaks(self, **kargs):
		d = {}
		if self.kaks is None:
			return d
		for kaks in KaKsParser(self.kaks, **kargs):
#		for line in open(self.kaks):
#			temp = line.rstrip().split()
#			if temp[0] == 'Sequence':
#				if temp[1] == 'dS-YN00':
#					kargs['yn00'] = True
#				elif temp[1] == '4D_Sites':
#					kargs['fdtv'] = True
#				continue
#			if line.startswith('#') or not temp:
#				continue
#			kaks = KaKs(temp, **kargs)
#			ks = kaks	#.ks
			ks = kaks
			d[kaks.pair] = ks
			d[tuple(reversed(kaks.pair))] = ks
		return d

def anchors2bed(collinearity, gff, chrmap, left_anchors, right_anchors, outbed=sys.stdout):
	left_anchors = left_anchors.split(',')	# ordered
	right_anchors = right_anchors.split(',')
	
	left_gs, right_gs = set([]), set([])
	for block in Collinearity(collinearity, gff, chrmap):
		genes1, genes2 = block.genes
		sp1, sp2 = block.species
		gs1, gs2 = block.genes1, block.genes2	# with coord
		chr1, chr2 = block.chr1, block.chr2
		if set(genes1) & set(left_anchors) and set(genes1) & set(right_anchors):
			pass
		elif set(genes2) & set(left_anchors) and set(genes2) & set(right_anchors):
			genes1, genes2 = genes2, genes1
			sp1, sp2 = sp2, sp1
			gs1, gs2 = gs2, gs1
			chr1, chr2 = chr2, chr1
		else:
			continue
		d_map = {}
		for g1, g2 in zip(gs1, gs2):
			d_map[g1.id] = (g1, g2)
		for anchor in left_anchors:
			if anchor in d_map:	# longest
				left_g1, left_g2 = d_map[anchor]
				break
		for anchor in reversed(right_anchors):
			if anchor in d_map:
				right_g1, right_g2 = d_map[anchor]
				break
		left_g2, right_g2 = sorted([left_g2, right_g2], key=lambda x:x.start)
		g2_chr, g2_start, g2_end = left_g2.chr, left_g2.start, right_g2.end
		g2_range = '{}-{}'.format(left_g2.raw_gene, right_g2.raw_gene)
		g1_range = '{}-{}'.format(left_g1.raw_gene, right_g1.raw_gene)
		id = '{}:{}:{}'.format(sp2, g2_range, g1_range)
		line = [chr2.split('|')[-1], g2_start-1, g2_end, id, sp2, ]
		line = list(map(str, line))
		print('\t'.join(line), file=outbed)
		left_gs.add(left_g1)
		right_gs.add(right_g1)
		anchor_sp = sp1
		anchor_chr = chr1
	left_g1, left_g2 = min(left_gs, key=lambda x:x.start), max(right_gs, key=lambda x:x.end)
	g1_chr, g1_start, g1_end = left_g1.chr, left_g1.start, right_g1.end
	id = '{}:{}'.format(anchor_sp, g1_range)
	line = [anchor_chr.split('|')[-1], g1_start-1, g1_end, id, anchor_sp, ]
	line = list(map(str, line))
	print('\t'.join(line), file=outbed)
	

	
class Gff:
	def __init__(self, gff):
		self.gff = gff
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.gff):
			yield GffLine(line)
			#line = line.rstrip().split()
			#yield Gene(line)
	def get_sps(self, sps, fout):
		sps = set(sps)
		for line in self:
			if line.species in sps:
				line.write(fout)
	def get_genes(self):
		d = {}
		if self.gff is None:
			return d
		for line in self:
			d[line.gene] = line
		return d
	def get_indexed_genes(self):
		d_chrom = OrderedDict()
		for line in self:
			try: d_chrom[line.chrom] += [line]
			except KeyError: d_chrom[line.chrom] = [line]
		d_genes = OrderedDict()
		d_length = OrderedDict()
		d_length2 = OrderedDict()
		species = set([])
		for chrom, lines in list(d_chrom.items()):
			lines = sorted(lines, key=lambda x:x.start)
			for i, line in enumerate(lines):
				line.index = i
				d_genes[line.gene] = line
		#	d_chrom[chrom] = lines
			d_length[chrom] = line.end
			d_length2[(line.species, chrom)] = line.end, len(lines)
			species.add(line.species)
		self.d_length = d_length
		self.d_length2 = d_length2
		#self.d_chrom = d_chrom
		self.species = species
		return d_genes
	def get_index(self):
		d_genes = self.get_indexed_genes()
		d_index = {}
		for gene, line in d_genes.items():
			d_index[(line.chrom, line.index+1)] = line
		return d_index
	def to_wgdi(self, chrLst='chr.list', pep='pep.faa', cds='cds.fa', indir='.', outdir='wgdi', species=None, split=True):
		from creat_ctl import get_good_chrs
		self.gff = os.path.join(indir, self.gff)
		chrLst = os.path.join(indir, chrLst)
		pep = os.path.join(indir, pep)
		cds = os.path.join(indir, cds)

		good_chrs = get_good_chrs(chrLst, min_genes=200)
		d_pep = {rc.id:rc for rc in SeqIO.parse(pep, 'fasta')}
		d_cds = {rc.id:rc for rc in SeqIO.parse(cds, 'fasta')}
		mkdirs(outdir)
		d_handle = {}
		d_genes = self.get_indexed_genes()
		if species is None:
			species = self.species
		else:
			species = parse_species(species)
#		all_gff = '{}/all.gff'.format(outdir, )
		#open files
		for sp in species:
			if split:
				prefix = '{}/{}'.format(outdir, sp)
			else:
				prefix = '{}/{}'.format(outdir, 'all')
			gff = prefix + '.gff'
			lens = prefix + '.lens'
			cds = prefix + '.cds'
			pep = prefix + '.pep'
			if not split and d_handle:
				d_handle[sp] = list(d_handle.values())[0]
			else: 
				d_handle[sp] = open(gff, 'w'), open(lens, 'w'), None,None #open(cds, 'w'), open(pep, 'w')
		# gff
		for line in list(d_genes.values()):
		#	print >> sys.stderr, line
			sp = line.species
			if sp not in set(species):
				continue
			gff,_,cds, pep = d_handle[sp]
#			for d, h in zip((d_pep,d_cds), (pep, cds)):
#				try:
#					SeqIO.write(d[line.gene], h, 'fasta')
#				except KeyError: pass
			chrom = line.chrom
#			chrom = chrom[2:]
			line = [chrom, line.gene, line.start, line.end, line.strand, line.index+1, line.gene]
			print('\t'.join(map(str, line)), file=gff)
#			line[0] = line.chrom
#			print >> all_gff, '\t'.join(map(str, line))
		# lens
		for (sp, chrom), (g_len, bp_len) in list(self.d_length2.items()):
			if sp not in set(species):
				continue
			if chrom not in set(good_chrs):
				continue
			_, lens, _,_ = d_handle[sp]
#			chrom = chrom[2:]
			line = (chrom, g_len, bp_len)
			print('\t'.join(map(str, line)), file=lens)
		# close files
		for sp in species:
			for hd in d_handle[sp]:
				try: hd.close()
				except: pass
			
	def fetch(self, g1,g2):
		'''g1 and g2 is in the same chromosome'''
		d_chrom = {}
		for line in self:
			if line.id == g1:
				target_chrom = line.chrom
				g1_start = line.start
			if line.id == g2:
				g2_start = line.start
			try: d_chrom[line.chrom] += [line]
			except KeyError: d_chrom[line.chrom] = [line]
		assert g1_start < g2_start
		lines = d_chrom[target_chrom]
		lines = sorted(lines, key=lambda x:x.start)
		reach = 0
		for i, line in enumerate(lines):
			line.index = i
			if line.id == g1:
				reach = 1
			if reach:
				yield line
			if line.id == g2:
				reach = 0

	def to_chroms(self, species=None):
		d_chrom = OrderedDict()
		for line in self:
			if species is not None and line.species != species: # 目标物种
				continue
			try: d_chrom[line.chrom] += [line]
			except KeyError: d_chrom[line.chrom] = [line]
		chroms = []
		for chrom, lines in list(d_chrom.items()):
			lines = sorted(lines, key=lambda x:x.start)
			for i, line in enumerate(lines):
				line.index = i
			name = chrom
			chrom = Chromosome(lines)
			chrom.species = line.species
			chrom.name = name
			chroms += [chrom]
		return Chromosomes(chroms)
	def to_graph(self):
		G = GffGraph() #nx.DiGraph()
		d_chrom = OrderedDict()
		for line in self:
			try: d_chrom[line.chrom] += [line]
			except KeyError: d_chrom[line.chrom] = [line]
		for chrom, lines in list(d_chrom.items()):
			lines = sorted(lines, key=lambda x:x.start)
			for i, line in enumerate(lines):
				line.index = i
			path = lines #[line.id for line in lines]
			G.add_path(path)
		return G
class SyntenyGraph(nx.Graph):
	def __init__(self, *args, **kargs):
		super().__init__(*args, **kargs)
	def score_node(self, node):
		return sum([1/attr['weight'] for _, attr in self[node].items()])
	def score_path(self, path):
		return sum([self.score_node(node) for node in path])

class Path:
	def __init__(self, path):
		self.path = path
	def __iter__(self):
		return iter(self.path)
	def __repr__(self):
		return str(self.path)
	def __len__(self):
		return len(self.path)
	def __getitem__(self, index):
		if isinstance(index, int):
			return self.path[index]
		else:
			return self.__class__(self.path[index])

class GffGraph(nx.DiGraph):
	def __init__(self, *args, **kargs):
		super().__init__(*args, **kargs)
	def remove_internals(self, internals):
		for node in internals:
			if node not in self:
				continue
			predecessors = list(self.predecessors(node))
			successors = list(self.successors(node))
			for n1, n2 in itertools.product(predecessors, successors):
				self.add_edge(n1, n2)
			self.remove_node(node)
	def add_path(self, path):
		self.add_edges_from([path[i:i+2] for i in range(len(path)-1)])
	@property
	def starts(self):
		for node, pred in self.pred.items():
			if not pred:
				yield node
	def iter_chrom(self, node):	# linear
		return self.fetch_chrom(node)
	@property
	def chroms(self):
		for start in self.starts:
			yield self.iter_chrom(start)
	def to_wgdi(self, prefix):
		fgff = open(prefix+'.gff', 'w')
		flen = open(prefix+'.lens', 'w')
		for chrom in self.chroms:
			for i,node in enumerate(chrom):
				if i == 0:
					_chrom = node.chrom
					_spec  = node.species
				node.start, node.end = i*100+1, (i+1)*100
				node.chrom = _chrom
				node.gene = _spec + '|' + node.gene.split('|', 1)[-1]
				node.index = i+1
				node.to_wgdi(fgff)
			line = [node.chrom, node.end, node.index]
			flen.write('\t'.join(map(str, line))+'\n')
		fgff.close()
		flen.close()
	def to_idmap(self):
		fidmap = open('id_map.txt', 'w')
		for chrom in self.chroms:
			for i,node in enumerate(chrom):
				if i == 0:
					_chrom = node.chrom
				node.start, node.end = i*100+1, (i+1)*100
				node.chrom = _chrom
				node.index = i+1
				line = [node.id, node.id, node.id, node.id, node.chrom, node.start, node.end, node.strand, None]
				fidmap.write('\t'.join(map(str, line))+'\n')
		fidmap.close()
	def fetch_chrom(self, start, end=None, reverse=False):	# linear
		node = start
		yield node
		while True:
			suc = self.predecessors(node) if reverse else self.successors(node)
			suc = list(suc)
			if not suc:
				break
			node = suc[0]
			yield node
			if end and node == end:
				break
	def lazy_fetch_chrom(self, start, end, **kargs):	# linear
		chrom = list(self.fetch_chrom(start, end, reverse=True, **kargs))
		if chrom[-1] == end:
			return Path(chrom)
		chrom = list(self.fetch_chrom(start, end, reverse=False, **kargs))
		if chrom[-1] == end:
			return Path(chrom)
	def index(self):
		for start in self.starts:
			for i, node in enumerate(self.iter_chrom(start)):
				node.index = i
	def insert_path(self, n1, n2, path):
		self.remove_edge(n1, n2)
		self.add_path([n1] + list(path) + [n2])
	def to_gfa(self, fout):
		for node in self.nodes:
			line = ['S', node, '*']
			print('\t'.join(map(str, line)), file=fout)
		for node1, node2 in self.edges:
			line = ['L', node1, '+', node2, '+', '0M']
			print('\t'.join(map(str, line)), file=fout)

class GffLine:
	def __init__(self, line):
		self.line = lazy_decode(line)
		self._parse()
	def _parse(self):
		temp = self.line.rstrip().split('\t')
		chr, gene, start, end = temp[:4]
		try: strand = temp[4]
		except IndexError: strand = None
		try: start, end = list(map(int, [start, end]))
		except ValueError as e:
			# bed
			gene, start, end = start, end, gene
			try: strand = temp[5]
			except IndexError: strand = None
			try: start, end = list(map(int, [start, end]))
			except ValueError as e:
				print('Error in line:', temp, file=sys.stderr)
				raise ValueError(e)
		g = Gene((gene, chr, start, end, strand))
		self.chrom, self.gene, self.start, self.end, self.strand = \
			chr, gene, start, end, strand
		self.Gene = g
		self.id = gene
		try: self.species, self.raw_gene = gene.split('|', 1)
		except ValueError: pass
	def __hash__(self):
		return hash(self.id)
	def __str__(self):
		return self.id
	def __repr__(self):
		return self.id
	def __lt__(self, other):
		return self.id < other.id
	def __eq__(self, other):
		return self.id == other.id
	def write(self, fout):
		fout.write(self.line)
	def to_wgdi(self, fout):
		line = [self.chrom, self.gene, self.start, self.end, self.strand, self.index, self.id]
		self.line = '\t'.join(map(str, line)) + '\n'
		self.write(fout)

def get_gff(gff, species, fout):
	sps = {line.strip().split()[0] for line in open(species)}
	Gff(gff).get_sps(sps, fout)
def get_ks(ksfile, pairfile, outks, outpair, source='wgdi'): # for wgdi	
	pairs = {CommonPair(*line.strip().split()) for line in pairfile}
	i = 0
	got_pairs = set([])
	for line in open(ksfile):
		i += 1
		if i == 1:
			outks.write(line)
			continue
		temp = line.strip().split()
		pair = CommonPair(*temp)
		if pair in pairs and pair not in got_pairs and len(temp) == 6:
			outks.write(line)
			got_pairs.add(pair)
	for pair in pairs - got_pairs:
		pair.write(outpair)

def slim_tandem(tandem, pairs, outPairs):
	slim_genes = Tandem(tandem).slims()
	for pair in Pairs(pairs):
		if set(pair.pair) & slim_genes:
			continue
		pair.write(outPairs)
	
def split_pair(line, sep=None, parser=None):
	pair = tuple(line.rstrip().split(sep, 1))
	return parser(*pair)
	
class CommonPair(object):
	def __init__(self, *pair):
		self.pair = pair[:2]
	def __iter__(self):
		return iter(self.pair)
	def __getitem__(self, index):
		return self.pair[index]
	def __str__(self):
		return '{}-{}'.format(*self.pair)
	def __format__(self):
		return str(self)
	@property
	def key(self):
		return tuple(sorted(self.pair))
	def __eq__(self, other):
		try: return self.key == other.key
		except AttributeError:
			other = SpeciesPair(*other)
			return self.key == other.key
	def __lt__(self, other):
		return self.key < other.key
	def __hash__(self):
		return hash(self.key)
	def write(self, fout):
		self.line = '{}\t{}\n'.format(*self.pair)
		fout.write(self.line)
class SpeciesPair(CommonPair):
	def __init__(self, *pair):
		super(SpeciesPair, self).__init__(*pair)
class Pair(CommonPair): # gene pair
	def __init__(self, *pair):
		super(Pair, self).__init__(*pair)
		self.gene1, self.gene2 = self.pair
	@property
	def species1(self):
		return self.gene1.split('|')[0]
	@property
	def species2(self):
		return self.gene2.split('|')[0]
	@property
	def species(self):
		return SpeciesPair(self.species1, self.species2)
		
class Pairs(object):
	def __init__(self, pairs, sep=None, parser=Pair):
		self.pairs = pairs
		self.sep = sep # line seperator
		self.parser = parser
	def __iter__(self):
		return self._parse()
	def _parse(self):
		for line in open(self.pairs):
			line = lazy_decode(line)
			yield split_pair(line, self.sep, self.parser)
	def graph(self):
		G = nx.Graph()
		for pair in self:
			G.add_edge(*pair.pair)
		return G
	def subgraphs(self):
		G = self.graph()
		for cmpt in nx.connected_components(G):
			yield G.subgraph(cmpt)
	def slims(self):
		genes = set([])
		for sg in self.subgraphs():
			max_node, _ = max(list(sg.degree().items()), key=lambda x:x[1])
			genes = genes | (set(sg.nodes()) - set([max_node]))
		return genes
class Tandem(Pairs):
	def __init__(self, pairs, sep=',', parser=Pair):
		super(Tandem, self).__init__(pairs, sep, parser)
class SpeciesPairs(Pairs):
	def __init__(self, pairs, sep=None, parser=SpeciesPair):
		super(SpeciesPairs, self).__init__(pairs, sep, parser)

def block_length(collinearity, sp_pairs=None):
	if sp_pairs is None:
		prefix = collinearity
	else:
		prefix = sp_pairs
		
	if sp_pairs is not None: # parse species pair file
		sp_pairs = set(SpeciesPairs(sp_pairs))
	d_genes = {}
	for rc in Collinearity(collinearity):
		spp = rc.species
		if not sp_pairs is None and not spp in sp_pairs:
			continue
		try: d_genes[spp] += [rc.N]
		except KeyError: d_genes[spp] = [rc.N]
	#print [spp], sp_pairs	
	#for sp in sp_pairs:
	#	print sp
	#print d_genes
	prefix += '.block_length'
	datafile = prefix + '.density.data'
	outfig = prefix + '.density.pdf'
	with open(datafile, 'w') as f:
		print('{}\t{}'.format('pair', 'value'), file=f)
		for spp, values in list(d_genes.items()):
			print(spp, values)
			for value in values:
				for v in range(value):
					print('{}\t{}'.format(str(spp), value), file=f)
	rsrc = prefix + '.density.r'
	#xlabel, ylabel = 'Block length (gene number)', 'Percent of genes'
	xlabel, ylabel = 'Block length (gene number)', 'Cumulative number of genes'
	with open(rsrc, 'w') as f:
		print('''datafile = '{datafile}'
data = read.table(datafile, head=T)
library(ggplot2)
#p <- ggplot(data, aes(x=value, color=pair)) + geom_line(stat="density", size=1.15) + xlab('{xlabel}') + ylab('{ylabel}')  + scale_colour_hue(l=45)
p <- ggplot(data, aes(x=value, fill=pair)) + geom_histogram() + xlab('{xlabel}') + ylab('{ylabel}')
ggsave('{outfig}', p, width=12, height=7)
'''.format(datafile=datafile, outfig=outfig, xlabel=xlabel, ylabel=ylabel, ), file=f)
	cmd = 'Rscript {}'.format(rsrc)
	os.system(cmd)

		
class Segment:
	def __init__(self, genes):
		self.genes = genes
	def __iter__(self):
		return iter(self.genes)
	def __len__(self):
		return len(self.genes)
	def __hash__(self):
		return hash(self.key)
	def __str__(self):
		return '{}:{}-{}({})'.format(self.chrom, self.start, self.end, self.strand)
	def __getitem__(self, index):
		if isinstance(index, int):
			return self.genes[index]
		else:
			return self.__class__(self.genes[index])

	@property
	def key(self):
		return tuple(map(str, self.genes))
	@property
	def head(self):
		return self.genes[0]
	@property
	def tail(self):
		return self.genes[-1]
	@property
	def chrom(self):
		return self.head.chr
	@property
	def indices(self):
		return [gene.index for gene in self]
	@property
	def start(self):
		return min(self.indices)
	@property
	def end(self):
		return max(self.indices)
	@property
	def span(self):
		return self.end - self.start + 1
	@property
	def strand(self):
		if self.indices[0] > self.indices[-1]:
			return '-'
		return '+'
	def reverse(self):
		self.genes = self.genes[::-1]
	def distance(self, other):
		if not self.chrom == other.chrom:
			return None
		return other.start - self.end
	def min_distance(self, other):
		if not self.chrom == other.chrom:
			return None
		return max(self.start, other.start) - min(self.end, other.end)
	def overlap(self, other):
		if not self.chrom == other.chrom:
			return False
		return max(0, min(self.end, other.end) - max(self.start, other.start))
	def contains(self, other):
		if not self.chrom == other.chrom:
			return False
		if other.start >= self.start and other.end <=self.end:
			return True
		return False

class Chromosome(Segment):
	def __init__(self, genes):
		self.genes = genes
	@lazyproperty
	def name(self):
		return self.chrom

class Chromosomes:
	def __init__(self, chroms):
		self.names = [chr.name for chr in chroms]
		self.chroms = chroms
	def __iter__(self):
		return iter(self.chroms)
	def __len__(self):
		return sum([len(chr) for chr in self.chroms])
	def sort(self):
		from creat_ctl import sort_version
		d = dict(list(zip(self.names, self.chroms)))
		sorted_names = sort_version(self.names)
		self.names = sorted_names
		self.chroms = [d[name] for name in self.names]
		
		
def cluster_pairs(collinearity, logs='b'):
	import networkx as nx
	G = cluster_graph(collinearity, logs=logs)
	for cmpt in nx.connected_components(G):
		yield cmpt
def cluster_subgraphs(collinearity, logs='b', **kargs):
	G = cluster_graph(collinearity, logs=logs, **kargs)
	for cmpt in nx.connected_components(G):
		yield G.subgraph(cmpt)
def cluster_graph(collinearity, logs='b', **kargs): 	# logs: b: both, o: orthologs
	import networkx as nx
	G = nx.Graph()
	for rc in Collinearity(collinearity, **kargs):
		sp1,sp2 = rc.species
		if logs == 'o' and sp1 == sp2:
			continue
		G.add_edges_from(rc.pairs)
	return G
def cluster_add_outgroup(collinearities, orthogroup, outgroup, fout=sys.stdout, min_ratio=0):
	outgroup = set(parse_group(outgroup))
	logger.info('outgroup: {}'.format(outgroup))
	G = nx.Graph()
	for rc in XCollinearity(collinearities):
		sp1,sp2 = rc.species
		if sp1 in outgroup and sp2 in outgroup:
			continue
		if not (sp1 in outgroup or sp2 in outgroup):
			continue
		G.add_edges_from(rc.pairs)
	logger.info('{} nodes in outgroup graph'.format(len(G)))
	nog, ng = 0,0
	for og in OrthoMCLGroup(orthogroup):
		outgrp_genes = []
		for g in og.genes:
			if g in G:
				outgrp_genes += G.neighbors(g)
		nxg = len(og.genes)
		add_genes = []
		for g, count in list(Counter(outgrp_genes).items()):
			ratio = 1.0*count/nxg
			if ratio > min_ratio:
				add_genes += [g]
				ng += 1
		if add_genes:
			og.genes += sorted(add_genes)
			nog += 1
		og.write(fout)
	logger.info('add {} outgroup genes for {} orthogroups'.format(ng, nog))
def parse_group(groups):
	xgroup = []
	if isinstance(groups, str):
		groups = [groups]
	if groups is not None:
		for group in groups:
			if exists_and_size_gt_zero(group):
				xgroup += [line.strip().split()[0] for line in open(group)]
			else:
				xgroup += [group]
	return xgroup

def cluster_by_mcl(collinearities, orthologs=None, inflation=2, outgroup=None, ingroup=None, outpre='cluster'):
	ingroup = set(parse_group(ingroup))
	outgroup = set(parse_group(outgroup))
	logger.info('outgroup: {}'.format(outgroup))
	logger.info('ingroup: {}'.format(ingroup))

	network = '{}.network'.format(outpre)
	fout = open(network, 'w')
	np = 0
	i,j,k = 0,0,0
	for rc in XCollinearity(collinearities, orthologs=orthologs):
		sp1,sp2 = rc.species
		if sp1 == sp2:	# exclude paralogs
			i += 1
			continue
		if sp1 in outgroup or sp2 in outgroup:  # exclude outgoup
			j += 1
			continue
		if ingroup and not (sp1 in ingroup and sp2 in ingroup):	# only include ingroup
			k += 1
			continue
		for g1, g2 in rc.pairs:
			np += 1
			line = [g1, g2]
			if orthologs:
				line += [str(rc.oi)]
			fout.write('{}\n'.format('\t'.join(line)))
	fout.close()
	logger.info('excluded: {} paralogs, {} in outgroup, {} not in ingroup'.format(i,j,k))
	cluster = '{}.mcl'.format(outpre)
	cmd = '''mcl {input} --abc -I {inflation} -o - -te {ncpu} | \
awk 'NF>1' | awk '{{split($0,a,"\\t");sl=asort(a);for (i=1;i<=sl;i++){{printf("%s ", a[i])}}; printf "\\n"}}' | \
awk '{{print "SOG"NR": "$0}}' > {output}'''.format(
				inflation=inflation, input=network, output=cluster, ncpu=10)
	run_cmd(cmd, log=True, )
	nc = len([1 for line in open(cluster)])
	logger.info('{} syntenic gene pairs reslut in {} SOGs'.format(np, nc))

	
def test_closest(collinearity, kaks, spsd, min_size=0):
	ColinearGroups(collinearity, spsd, kaks=kaks, min_size=min_size).get_min_ks()
def cg_trees(collinearity, spsd, seqfile, gff, tmpdir='./tmp'):
	ColinearGroups(collinearity, spsd, seqfile=seqfile, gff=gff, tmpdir=tmpdir).chrom_trees() #get_trees()
def anchor_trees(collinearity, spsd, seqfile, gff, tmpdir='./tmp'):
	ColinearGroups(collinearity, spsd, seqfile=seqfile, gff=gff, min_size=5, tmpdir=tmpdir).anchor_trees() #get_trees()
def gene_trees(collinearity, spsd, seqfile, orthologs, tmpdir='./tmp'):
	ColinearGroups(collinearity, spsd, seqfile=seqfile, orthologs=orthologs, tmpdir=tmpdir).get_trees()
def to_phylonet(collinearity, spsd, seqfile, outprefix, tmpdir='./phylonet_tmp'):
	ColinearGroups(collinearity, spsd, seqfile=seqfile, tmpdir=tmpdir).to_phylonet(outprefix)
def to_ark(collinearity, spsd, gff, max_missing=0.2):
	ColinearGroups(collinearity, spsd, gff=gff).to_ark(max_missing=max_missing)
def gene_retention(collinearity, spsd, gff):
	ColinearGroups(collinearity, spsd, gff=gff).gene_retention()

GenetreesTitle = ['OG', 'genes', 'genetree', 'min_bootstrap', 'topology_species',
				'chromosomes', 'topology_chromosomes']

class ColinearGroups:
	def __init__(self, collinearity=None, spsd=None, 
				kaks=None, seqfile=None, gff=None, 
				min_size=0, tmpdir='./tmp', 
				orthologs=None, 	# 直系同源关系。共线性的备选，无基因组或染色体时使用
				noparalog=True,	# no paralogs
				nosamechr=False,	# no same chromosome
				):
		self.collinearity = collinearity
		self.kaks = kaks
		self.seqfile = seqfile
		self.gff = gff
		self.min_size = min_size
		self.tmpdir = tmpdir
		self.orthologs = orthologs
		self.noparalog = noparalog
		self.nosamechr = nosamechr
		sp_dict = parse_spsd(spsd)
		self.sp_dict = sp_dict #Counter(sp_dict)
		self.spsd = spsd
		self.max_ploidy = max(list(sp_dict.values())+[1])
		self.prefix = spsd
		#print >>sys.stderr, self.sp_dict
	@property
	def groups(self):
		G = nx.Graph()
		for rc in Collinearity(self.collinearity):
			if self.noparalog and len(set(rc.species)) == 1: # discard paralog
				continue
			if self.nosamechr and rc.chr1 == rc.chr2:
				continue
			for pair in rc.pairs:
				G.add_edge(*pair)
		i = 0
		for cmpt in nx.connected_components(G):
			i += 1
			ogid = 'SOG{}'.format(i)
			yield OrthoMCLGroupRecord(genes=cmpt, ogid=ogid)
	def to_synet(self, fout=sys.stdout):
		d_profile = dict([(sp, []) for sp in list(self.sp_dict.keys())])
		i = 0
		for group in self.groups:
			i += 1
			counter = group.counter
			for sp in list(d_profile.keys()):
				value = '1' if sp in counter else '0'
				d_profile[sp] += [value]
		desc = 'ntaxa={};ncluster={}'.format(len(d_profile), i)
		for sp, values in list(d_profile.items()):
			print('>{} {}\n{}'.format(sp, desc, ''.join(values)), file=fout)
					
	def infomap(self):
		mkdirs(self.tmpdir)
		d_id = {}
		i = 0
		graphfile = '{}/infomap.graph'.format(self.tmpdir)
		f = open(graphfile, 'w')
		for rc in Collinearity(self.collinearity):
			if len(set(rc.species)) == 1: # discard paralog
				continue
			for g1, g2 in rc.pairs:
				for g in [g1, g2]:
					if not g in d_id:
						i += 1
						d_id[g] = i
				i1, i2 = d_id[g1], d_id[g2]
				print('{} {}'.format(i1, i2), file=f)
				print('{} {}'.format(i2, i1), file=f)
		f.close()
		cmd = 'infomap {} {} --clu -N 10  -2'.format(graphfile, self.tmpdir)
		run_cmd(cmd)
					
	@property
	def raw_graph(self):
		G = nx.Graph()
		sp_pairs = set([])
		for rc in Collinearity(self.collinearity, kaks=self.kaks):
			if rc.N < self.min_size:	# min length
				continue
			sp_pairs.add(rc.species)
			for pair in rc.pairs:
				G.add_edge(*pair)
		if self.orthologs is not None:
			for pair in Pairs(self.orthologs):
				if pair.species in sp_pairs:
					continue
				G.add_edge(*pair.pair)
		return G
	def gene_retention(self, winsize=100, winstep=None, min_genes=0.02):
		if winstep is None:
			winstep = winsize/10
		self.root = self.get_root()
		target_sps = sorted(set(self.sp_dict)-set([self.root]), key=lambda x: list(self.sp_dict.keys()).index(x))
		d_sp = OrderedDict([(sp, []) for sp in target_sps])
#		d_retention = copy.deepcopy(d_sp)
		sp_comb = [(sp1, sp2) for sp1, sp2 in itertools.combinations(target_sps, 2)]
#		d_diff = OrderedDict([(spc, []) for spc in sp_comb])
#		d_loss = copy.deepcopy(d_sp)
		# out
		out_rete = self.prefix + '.retention'
		out_diff = self.prefix + '.diff'
		out_loss = self.prefix + '.loss'
		f_rete = open(out_rete, 'w')
		line = ['ichr', 'chr', 'win', 'sp', 'retention']
		print('\t'.join(line), file=f_rete)
		f_diff = open(out_diff, 'w')
		line = ['ichr', 'chr', 'win', 'spc', 'diff']
		print('\t'.join(line), file=f_diff)
		f_loss = open(out_loss, 'w')
		line = ['ichr', 'chr', 'sp', 'loss']
		print('\t'.join(line), file=f_loss)
		
		gff = Gff(self.gff)
		chroms = gff.to_chroms(species=self.root)
		chroms.sort()
		graph = self.raw_graph
		ichr = 0
		for chrom in chroms:
			if 1.0* len(chrom)/len(chroms) < min_genes:	# too short chrom
				continue
			ichr += 1
			for i in range(0, len(chrom), winstep):
				window = chrom[i: i+winsize]
				d_win = copy.deepcopy(d_sp)
				size = len(window)
				if size < winsize/2:
					continue
				d_win = self.count_window(window, graph, d_win)
				for sp, counts in list(d_win.items()):
					retention = [v for v in counts if v>0]
					rate = 1e2*len(retention) / size
#					d_retention[sp] += [rate]
					line = [ichr, chrom.name, i, sp, rate]
					line = list(map(str, line))
					print('\t'.join(line), file=f_rete)
				for sp1, sp2 in sp_comb:
					counts1, counts2 = d_win[sp1], d_win[sp2]
					try: diff = self.count_diff(counts1, counts2)
					except ZeroDivisionError: continue
					line = [ichr, chrom.name, i, sp1+'-'+sp2, diff]
					line = list(map(str, line))
					print('\t'.join(line), file=f_diff)
			d_win = copy.deepcopy(d_sp)
			d_win = self.count_window(chrom, graph, d_win)
			for sp, counts in list(d_win.items()):
				for loss in self.count_loss(counts):
					line = [ichr, chrom.name, sp, loss]
					line = list(map(str, line))
					print('\t'.join(line), file=f_loss)
		f_rete.close()
		f_loss.close()
		f_diff.close()
	def count_window(self, window, graph, d_win):
		for gene in window:
			for sp in list(d_win.keys()):
				d_win[sp] += [0]
			if not gene.id in graph:
				continue
			for syn_gene in graph[gene.id]:
				sp = syn_gene.split('|')[0]
				if sp not in d_win:
					continue
				d_win[sp][-1] += 1
		return d_win
	def count_diff(self, counts1, counts2):
		retent, diff, loss = 0, 0,0
		for v1, v2 in zip(counts1, counts2):
			if v1 == v2 == 0:
				loss += 1
			elif v1 ==0 or v2 == 0:
				diff += 1
			else:
				retent += 1
		return 1e2*diff/(diff+retent+loss)
		
	def count_loss(self, counts):
		last_v, last_i = '', 0
		for i,v in enumerate(counts):
			if last_v == 0 and v > 0:
				yield i - last_i
			if v == 0 and last_v > 0:
				last_i = i
			last_v = v
		
	@property
	def graph(self):
		G = nx.Graph()
		d_ks = {}
		sp_pairs = set([])
		for rc in Collinearity(self.collinearity, kaks=self.kaks):
			if rc.N < self.min_size:	# min length
				continue
			if set(rc.species) - set(self.sp_dict):	# both be in sp_dict
				continue
			if len(set(rc.species)) == 1: # discard paralog
				continue
			sp_pairs.add(rc.species)
			for pair, ks in zip(rc.pairs, rc.ks):
				G.add_edge(*pair)
				key = tuple(sorted(pair))
				d_ks[key] = ks
		#print sp_pairs
		self.d_ks = d_ks
		if self.orthologs is not None:	# 直系同源关系，无共线性信息
			for pair in Pairs(self.orthologs):
				if set(pair.species) - set(self.sp_dict):	# 只保留目标物种
					continue
				if pair.species in sp_pairs:	# 只导入无共线性信息的
					continue
				G.add_edge(*pair.pair)
		return G
	def filter_blocks(self):
		for rc in Collinearity(self.collinearity):
			if rc.N < self.min_size:	# min length
				continue
			if self.noparalog and len(set(rc.species)) == 1:
				continue
			if self.nosamechr and rc.chr1 == rc.chr2:
				continue
			if self.sp_dict and set(rc.species) - set(self.sp_dict):
				continue
			yield rc
	def to_graph(self):
		G = SyntenyGraph()
		for rc in self.filter_blocks():
			for pair in rc.pairs:
				G.add_edge(*pair, weight=1/rc.score)
		return G
	@property
	def chr_graph(self):
		G = nx.Graph()
		for rc in Collinearity(self.collinearity):
			if rc.N < self.min_size:	# 至少20
				continue
			if set(rc.species) - set(self.sp_dict):	# both be in sp_dict
				continue
			chr1 = (rc.species1, rc.chr1)
			chr2 = (rc.species2, rc.chr2)
			G.add_edge(chr1, chr2)
		return G
	def chr_subgraphs(self, min_tile=0.2, min_count=15):
		from creat_ctl import sort_version
		G = nx.Graph()
		for sg in self.subgraphs(same_number=False, same_degree=False, max_missing=0):
			chroms = [(gene2species(gene), self.d_gff[gene].chrom) for gene in sg.nodes()]
			for chr1, chr2 in itertools.combinations(chroms, 2):
				try: 
					G[chr1][chr2]['count'] += 1	# 对染色体组合进行计数
				except KeyError:
					G.add_edge(chr1, chr2, count=1)
		#for n1,n2 in G.edges():
		#	print n1,n2, 
		#	print G[n1][n2]
		counts = [G[n1][n2]['count'] for n1,n2 in G.edges()]
		#min_count = np.percentile(counts, min_tile)
		print('min_count of cluster', min_count)
		for cmpt in nx.connected_components(G):
			cmpt = sorted(cmpt)
			sps = [sp for sp,chr in cmpt]
			sps_count = Counter(sps)
			less = False
			print(sps_count)
			for sp, count in list(self.sp_dict.items()):
				if sps_count.get(sp, 0) < count:
					less = True
					break
			if less:	# 每个物种的染色体数目不得小于倍性
				continue
			#print cmpt
			
			d_count = {}
			groups = []
			less = False
			for sp, group in itertools.groupby(cmpt, key=lambda x:x[0]):
				#print list(group)
				# 物种内按倍性组合
				combs = itertools.combinations(group, self.sp_dict[sp])
				flt_combs = []
				for comb in combs: # 过滤掉无连接的
					#print comb
					counts = []
					for chr1, chr2 in itertools.combinations(comb, 2):
						try: count = G[chr1][chr2]['count']
						except KeyError: count = 0
						#print chr1, chr2, count
						if count < min_count:
							break
						counts += [count]
					else:
						print(comb, counts)
						flt_combs += [comb]
				#print flt_combs
				d_count[sp] = len(flt_combs)
				groups += [flt_combs]
				if not flt_combs:
					less = True
					break
			#print groups
			#for group in groups:
			#	for g1 in group:
			#		print g1
			print(d_count)
			if less: # 如果某个物种缺失，则无法形成染色体组合
				continue
			#products = itertools.product(*groups)
			#print len(list(products))
			i = 0
			for group in itertools.product(*groups):
				comb = list(flatten(group))	# 染色体组合
				counts = []
				less = False
				for chr1, chr2 in itertools.combinations(comb, 2):
					try: count = G[chr1][chr2]['count']
					except KeyError: count = 0
					if count < min_count:
						less = True
						break
					counts += [count]
				else:
					print(comb, counts)
				if less: # 一对染色体数量不足，则弃去这组合
					continue
				i += 1
				#print group
				chroms = [chr for sp, chr in comb]
				yield sort_version(chroms)
			#	yield group
			print(i)
	def anchor_trees(self):
		#self.chr_subgraphs()
		# chromosome tree
		max_trees = 10
		i = 0
		j = 0
		cmd_list = []
		treefiles = []
		treefiles2 = []
		d_gene_count = {}
		d_gene_count2 = {}
		#print >> sys.stdout, len(d_chromfiles), 'chromosome groups'
		# 按染色体串联建树，允许部分基因丢失
		for chroms in self.chr_subgraphs():
			j += 1
			if len(chroms) > len(set(chroms)):
				continue
			i += 1
			if i > max_trees:
				continue
			alnfiles = self.chrom_tree(chroms)
			ngene = len(alnfiles)
			print(len(self.d_chroms), list(self.d_chroms.items())[:10], file=sys.stderr)
			prefix = 'CHR_' + '-'.join(chroms) + '_'+str(ngene) + '_' + str(len(alnfiles))
			cmds = self.concat_tree(alnfiles, prefix, idmap=self.d_chroms, astral=True)
			treefile = self.iqtree_treefile
			treefiles += [treefile]
			d_gene_count[treefile] = len(alnfiles)
			# astral
			treefile = self.astral_treefile
			treefiles2 += [treefile]
			d_gene_count2[treefile] = len(alnfiles)
			
			cmd_list += [cmds]
			#print >> sys.stderr, 'dot', chroms
#			cmd_list += self.dot_plot(chroms)
			#print >> sys.stderr, prefix, len(alnfiles)
		print(j, 'chromosome groups', file=sys.stdout)
		cmd_file = '{}/chrom-cmds.list'.format(self.tmpdir)
		if cmd_list:
			run_job(cmd_file, cmd_list=cmd_list, tc_tasks=100)
		print(i, 'chromosome groups', file=sys.stdout)
		print('# iqtree', file=sys.stdout)
		self.print_topology(treefiles, d_gene_count=d_gene_count)
		print('# astral', file=sys.stdout)
		self.print_topology(treefiles2, d_gene_count=d_gene_count2)
		
		# clean
		self.clean(self.tmpdir)
		
	def subgraphs(self, same_number=True, same_degree=False, max_missing=0.2):
		'''same_number和max_missing互斥：当same_number为真时，无missing；
		当same_number为假时，由max_missing控制物种缺失率；
		max_missing=0不允许缺失物种'''
		self.count = []
		G = self.graph
		for cmpt in nx.connected_components(G):
			sg = G.subgraph(cmpt)
			if self.orthologs is None:
				try:
					chroms = [self.d_gff[gene].chrom for gene in sg.nodes()]
				except KeyError:
					chroms = [d_gff[gene].chrom for gene in sg.nodes() if gene in self.d_gff]

				if not len(chroms) == len(set(chroms)):	# 分列于不同的染色体
					continue
			sp_count = d_count = Counter(genes2species(sg.nodes()))
		#	print d_count
			if len(d_count) == len(self.sp_dict):	# same species set
				self.count += [tuple(sorted(d_count.items()))]
		#	if same_number and not len(sg.nodes()) == sum(self.sp_dict.values()):# 基因总数符合
		#		continue
			if same_number and not d_count == self.sp_dict:  # 基因数量和倍性吻合
				continue
			d_degree = sg.degree()
			if same_degree and not min(d_degree.values()) == len(sg.nodes())-1:	 # 互连
				continue
			if not same_number:	# 限制物种缺失率，不允许过多缺失。
				target_sps = [sp for sp,count in list(sp_count.items()) if 0<count<=self.sp_dict[sp]]
				missing = 1 - 1.0*len(target_sps) / len(self.sp_dict)
				if missing > max_missing:
					continue
				sps = [gene2species(gene) for gene in cmpt]
				target_sps = set(target_sps)
				genes = [gene for sp, gene in zip(sps, cmpt) if sp in target_sps]
				sg = G.subgraph(genes)
			yield sg
	def to_ark(self, fout=sys.stdout, outfmt='grimm', min_genes=200, max_missing=0.2):
		from creat_ctl import is_chr0, sort_version
		logger.info('loading collinear graph')
		d_idmap = {}
		i = 0
		mapfile = '{}.groups'.format(self.spsd)
		fmapout = open(mapfile, 'w')
		for sg in self.subgraphs(same_number=False, same_degree=False, max_missing=max_missing): # 0.8 before 2020-6-29
			i += 1
			genes = sg.nodes()
			for gene in genes:
				d_idmap[gene] = i
			group = OrthoMCLGroupRecord(ogid=i, genes=sorted(genes))
			group.write(fmapout)
		fmapout.close()

		logger.info('loading gff')
		d_chroms = {}
		for line in Gff(self.gff):
			sp, chrom = line.species, line.chrom
			if not sp in self.sp_dict:
				continue
			gene = line.Gene
			try: d_chroms[sp][chrom] += [gene]
			except KeyError: 
				try: d_chroms[sp][chrom] = [gene]
				except KeyError: d_chroms[sp] = {chrom: [gene]}
		
		logger.info('output markers')
		for sp in self.sp_dict:
			print('>{}'.format(sp), file=fout)
			d_chrs = d_chroms[sp]
			chroms = list(d_chrs.keys())
			chroms = sort_version(chroms)
			print('>{}'.format(sp), file=sys.stderr)
			total = 0
			for chrom in chroms:
				if is_chr0(chrom):
					continue
				genes = d_chrs[chrom]
				if len(genes) < min_genes:
					continue
				genes = sorted(genes, key=lambda x:x.start)
				markers = []
				for gene in genes:
					if not gene.id in d_idmap:
						continue
					marker = str(d_idmap[gene.id])
					if gene.strand == '-':
						marker = '-' + marker
					markers += [marker]
				print(' '.join(markers) + ' $', file=fout)
				print(chrom, len(markers), file=sys.stderr)
				total += len(markers)
			print('total', total, file=sys.stderr)
	def to_phylonet(self, outprefix, min_ratio=0.9):
		'''使用多拷贝基因用于phylonet'''
		mkdirs(self.tmpdir)
		self.d_seqs = d_seqs = seq2dict(self.seqfile)
		self.root = root_sp = self.get_root()
		G = self.graph
		d_idmap = {}
		d_idmap2 = {}
		treefiles = []
		cmd_list = []
		i,j = 0,0
		for genes in nx.connected_components(G):
			sps = [gene2species(gene) for gene in genes]
			sp_count = Counter(sps)
			target_sps = [sp for sp,count in list(sp_count.items()) if 0<count<=self.sp_dict[sp]]
			if 1.0*len(target_sps) / len(self.sp_dict) < min_ratio:
				continue
			
			target_sps = set(target_sps)
			if not self.root in target_sps:
				continue
			i += 1
			#print >>sys.stderr, i, sp_count
			
			target_genes = [gene for sp, gene in zip(sps, genes) if sp in target_sps]
			
			og = 'OG_{}'.format(i)
			outSeq = '{}/{}.fa'.format(self.tmpdir, og)
			root = None
			d_num = {sp:0 for sp in target_sps}
			fout = open(outSeq, 'w')
			for gene in sorted(target_genes):
				#j += 1
				rc = d_seqs[gene]
				sp = gene2species(rc.id)
				j = d_num[sp] + 1
				d_num[sp] = j
				sid = '{}.{}'.format(sp, j)
				gid = format_id_for_iqtree(rc.id)
				d_idmap2[gid] = sid
				try: d_idmap[sp] += [sid]
				except KeyError: d_idmap[sp] = [sid]
				rc.id = gid
				SeqIO.write(rc, fout, 'fasta')
				if sp == root_sp:
					root = rc.id
			fout.close()
			
			cmds = []
			alnSeq = outSeq + '.aln'
			alnTrim = alnSeq + '.trimal'
			iqtreefile = alnTrim + '.treefile'
			treefile = rooted_treefile = iqtreefile
			treefiles += [treefile]
			if not os.path.exists(iqtreefile):
				cmd = 'mafft --auto {} > {} 2> /dev/null'.format(outSeq, alnSeq)
				cmds += [cmd]
				cmd = 'trimal -automated1 -in {} -out {} &> /dev/null'.format(alnSeq, alnTrim)
				cmds += [cmd]
				opts = ''
				if not root is None:
					opts = '-o {}'.format(root)
				cmd = 'iqtree -redo -s {} -nt AUTO -bb 1000 {} -mset JTT &> /dev/null'.format(alnTrim, opts)
				cmds += [cmd]
			cmds = ' && '.join(cmds)
			cmd_list += [cmds]
		run_job(cmd_list=cmd_list, tc_tasks=100)
		genetrees = '{}.genetrees'.format(outprefix)
		self.cat_genetrees(treefiles, genetrees, idmap=d_idmap2, plain=False, format_confidence='%d')
		taxamap = '{}.taxamap'.format(outprefix)
		with open(taxamap, 'w') as fout:
			print(self.to_taxa_map(d_idmap), file=fout)
		
		self.clean(self.tmpdir)
	def to_taxa_map(self, d_idmap):	# PHYLONET Taxa Map
		map = []
		for sp, indvs in list(d_idmap.items()):
			indvs = sorted(set(indvs))
			indvs = ','.join(indvs)
			map += ['{}:{}'.format(sp, indvs)]
		return '<{}>'.format(';'.join(map))
	@lazyproperty
	def d_gff(self):
		return Gff(self.gff).get_genes()

			
	
	def get_trees(self):	# gene trees
		'''完全符合倍性比的基因树'''
		from creat_ctl import sort_version
		if not os.path.exists(self.tmpdir):
			os.mkdir(self.tmpdir)
		self.d_gff = d_gff = Gff(self.gff).get_genes()
		#print >> sys.stderr, d_gff.items()[:100]
		self.d_seqs = d_seqs = seq2dict(self.seqfile)
		self.root = root_sp = self.get_root()
		d_species = {}
		cmd_list = []
		treefiles = []
		iqtreefiles = []
		i = 0
		chrom_lists = []
		d_chroms = {}
		d_alnfiles = {}
		gene_groups = []
		chrom_groups = []
		ogs = []
		for sg in self.subgraphs():
			genes = sg.nodes()
			try:	
				chroms = [d_gff[gene].chrom for gene in genes]
				if not len(chroms) == len(set(chroms)):	# 一条染色体一个基因
					continue
			except KeyError:
				chroms = [d_gff[gene].chrom for gene in genes if gene in d_gff]
			i += 1
			og = 'OG_{}'.format(i)
			ogs += [og]
			gene_groups += [genes]
			outSeq = '{}/{}.fa'.format(self.tmpdir, og)
			chroms = tuple(sort_version(chroms))
			chrom_lists += [chroms]
			root = None
			fout = open(outSeq, 'w')
			for gene in genes:
				rc = d_seqs[gene]
				sp = gene2species(rc.id)
				try:
					chrom = d_gff[gene].chrom
				except KeyError:
					chrom = None
				chrom_id = '{}-{}'.format(sp, chrom)
				rc.id = format_id_for_iqtree(gene)		# 被改变了
				d_species[rc.id] = sp
				d_chroms[rc.id] = chrom_id
				d_species[chrom_id] = sp
				d_species[chrom] = sp
				SeqIO.write(rc, fout, 'fasta')
				if sp == root_sp:
					root = rc.id
			fout.close()
#			if root is not None:
#				chrom_groups += [tuple(sort_version(set(chroms)-{d_chroms[root]}))]
			cmds = []
			alnSeq = outSeq + '.aln'
			alnTrim = alnSeq + '.trimal'
			iqtreefile = alnTrim + '.treefile'
			treefile = rooted_treefile = alnTrim + '.tre'
			treefiles += [treefile]
			iqtreefiles += [iqtreefile]
			d_alnfiles[alnTrim] = chroms
			if not os.path.exists(iqtreefile):
				cmd = '. ~/.bashrc; mafft --auto {} > {} 2> /dev/null'.format(outSeq, alnSeq)
				cmds += [cmd]
				cmd = 'trimal -automated1 -in {} -out {} &> /dev/null'.format(alnSeq, alnTrim)
				cmds += [cmd]
				opts = ''
				if not root is None:
					opts = '-o {}'.format(root)
				cmd = 'iqtree -redo -s {} -nt AUTO -bb 1000 {} -mset JTT &> /dev/null'.format(alnTrim, opts)
				cmds += [cmd]
			if not test_s(rooted_treefile):
				if root is None:
					cmd = 'nw_reroot {} '.format(iqtreefile)
				else:
					cmd = 'nw_reroot {intre} {root} | nw_prune - {root} '.format(
						intre=iqtreefile, root=root,)
				cmd += ' | nw_topology -I - | nw_order - | nw_order - -c d | nw_order - > {}'.format(rooted_treefile)
			else:
				cmd = ''
			#if not os.path.exists(iqtreefile):
			cmds += [cmd]
			cmds = ' && '.join(cmds)
			cmds += '\nrm '+outSeq
			cmd_list += [cmds]
		if cmd_list:
			cmd_file = '{}/cmds.list'.format(self.tmpdir)
			run_job(cmd_file, cmd_list=cmd_list, tc_tasks=100)
		#print >>sys.stderr, Counter(chrom_lists)
		#print >>sys.stderr, Counter(self.count)
		d_count = Counter(self.count)
		self.print_self(d_count)
		self.d_species = d_species		# gene id / chrom_id -> sp
		self.d_chroms = d_chroms		# gene id -> chrom_id
		self.d_alnfiles = d_alnfiles	# alnfile -> chrom list
		print(i, 'groups', file=sys.stdout)
		self.print_topology(treefiles)
		# clean
		f = open('genetrees.list', 'w')
		line = GenetreesTitle
		print('\t'.join(line), file=f)
		j = 0
		for og, genes, chroms, treefile, iqtreefile in zip(ogs, gene_groups, chrom_lists, treefiles, iqtreefiles):
			#print chroms
			genes = ','.join(genes)
			chroms = ','.join(chroms) if chroms else ''
			if test_s(iqtreefile):
				tree = [line.strip() for line in open(iqtreefile)][0]
				min_bs = self.get_min_bootstrap(iqtreefile)
				try: topl = self.get_topology(treefile)
				except ValueError as e: 
					logger.warn('{}: {}'.format(treefile, e))
					continue
				topl_chr = self.get_topology(treefile, idmap=self.d_chroms)
			else:
				continue
			j += 1
			line = [og, genes, tree, min_bs, topl, chroms, topl_chr]
			line = list(map(str, line))
			print('\t'.join(line), file=f)
		f.close()
		print(j, 'groups with treefile', 1e2*j/i, file=sys.stdout)
		self.clean(self.tmpdir)
		return treefiles
	def print_self(self, d_count):
		'''统计基因比'''
		my_counts = []
		for key, count in list(d_count.items()):
			sps = [sp for sp in self.sp_dict]
			dcounts = dict([(sp, _count) for sp, _count in key])
			counts = [dcounts[sp] for sp in sps]
			my_counts += [[counts, count]]
		print(sps, file=sys.stdout)
		for counts, count in sorted(my_counts):
			print(counts, count, file=sys.stdout)
	def chrom_tree(self, target_chroms, min_ratio=0.3, min_seqs=4):
		'''按染色体串联（允许丢失）'''
		prefix = '-'.join(target_chroms)
		d_gff = self.d_gff
		try: d_seqs = self.d_seqs
		except AttributeError: self.d_seqs = d_seqs = seq2dict(self.seqfile)
		self.root = root_sp = self.get_root()
		
		d_species = {}
		cmd_list = []
		i = 0
		j = 0
		d_chroms = {}
		d_alnfiles = {}
		for sg in self.subgraphs(same_number=False):
			chroms = [d_gff[gene].chrom for gene in sg.nodes()]
			if not set(chroms) & set(target_chroms):
				continue
			i += 1
			
			og = '{}_OG_{}'.format(prefix, i)
			d_count = Counter(chroms)
			diff = set(chroms) - set(target_chroms)
			ratio = 1e0 * len(chroms) / len(target_chroms)
			
			if diff:
				continue
			if not max(d_count.values())<2:
				continue
			j += 1
			if len(chroms) < min_seqs:
				continue
			if not ratio>min_ratio:	# 子集
				continue
			#print >>sys.stderr, chroms
			outSeq = '{}/gene-{}.fa'.format(self.tmpdir, og)
			root = None
			fout = open(outSeq, 'w')
			for gene in sg.nodes():
				rc = d_seqs[gene]
				sp = gene2species(gene) #gene2species(rc.id)
				chrom = d_gff[gene].chrom
				chrom_id = '{}-{}'.format(sp, chrom)
				rc.id = format_id_for_iqtree(gene)
				d_species[rc.id] = sp
				d_chroms[rc.id] = chrom_id
				d_species[chrom_id] = sp
				d_species[chrom] = sp
				SeqIO.write(rc, fout, 'fasta')
				if sp == self.root:
					root = rc.id
			fout.close()
			
			cmds = []
			alnSeq = outSeq + '.aln'
			alnTrim = alnSeq + '.trimal'
			d_alnfiles[alnTrim] = chroms
			if True: #not os.path.exists(alnTrim):
				cmd = '. ~/.bashrc; mafft --auto {} > {} 2> /dev/null'.format(outSeq, alnSeq)
				cmds += [cmd]
				cmd = 'trimal -automated1 -in {} -out {} &> /dev/null'.format(alnSeq, alnTrim)
				cmds += [cmd]
				opts = ''
				if not root is None:
					opts = '-o {}'.format(root)
				cmd = 'iqtree -redo -s {} -nt 1  {} -mset JTT &> /dev/null'.format(alnTrim, opts)
				cmds += [cmd]
			#cmds += ['rm '+outSeq]
			cmds = ' && '.join(cmds)
			cmds += '\nrm '+outSeq
			cmd_list += [cmds]
		print(prefix, '{} used / {} available / {} genes'.format(len(d_alnfiles), j, i), file=sys.stdout)
		if cmd_list:
			cmd_file = '{}/gene-{}-aln-cmds.list'.format(self.tmpdir, prefix)
			run_job(cmd_file, cmd_list=cmd_list, tc_tasks=100, by_bin=1)
		
		alnfiles = list(d_alnfiles.keys())
		self.d_species = d_species		# gene id / chrom_id -> sp
		self.d_chroms = d_chroms		# gene id -> chrom_id
		self.d_alnfiles = d_alnfiles	# alnfile -> chrom list
		#print set(d_species.values())
		#print set(d_chroms.values())
		return alnfiles
		

	def chrom_trees(self, min_genes=2):
		'''按染色体串联（各种都做，包括允许基因丢失）'''
		self.treefiles = self.get_trees()	# 基因树，完全符合倍性比
		#self.print_topology(self.treefiles)
		d_chromfiles = {}
		for alnfile, chrom in list(self.d_alnfiles.items()):
			try: d_chromfiles[chrom] += [alnfile]
			except KeyError: d_chromfiles[chrom] = [alnfile]
	#	print >> sys.stderr, set(map(len, d_chromfiles.keys()))
		cmd_list = []
		treefiles = []
		d_gene_count = {}
		treefiles2 = []  # astral
		d_gene_count2 = {}
		d_gene_count3 = {}
		i = 0
		xxchroms = []
		d_concat_alnfiles = {}
		# 按染色体串联建树，只用完全符合倍性比的基因，至少俩基因
		for chroms, alnfiles in sorted(list(d_chromfiles.items()), key=lambda x: -len(x[1])):
			if len(chroms) > len(set(chroms)) or len(alnfiles) < min_genes:
				continue
			xxchroms += [chroms]
			i += 1
			prefix = '-'.join(chroms) + '_' + str(len(alnfiles))
			
			cmds = self.concat_tree(alnfiles, prefix, idmap=self.d_chroms, astral=True)
		#	d_gene_count3[concat_alnfile] = len(alnfiles)
			print(prefix, len(alnfiles), file=sys.stdout) #, alnfiles
			treefile = self.iqtree_treefile
			d_gene_count[treefile] = len(alnfiles)
			treefiles += [treefile]
			# astral
			treefile = self.astral_treefile
			treefiles2 += [treefile]
			d_gene_count2[treefile] = len(alnfiles)

			cmd_list += [cmds]
			
		cmd_file = '{}/merged-cmds.list'.format(self.tmpdir)
		if cmd_list:
			run_job(cmd_file, cmd_list=cmd_list, tc_tasks=50)
		print(sum(d_gene_count.values()), 'groups', '/', i, 'clusters', file=sys.stdout)
		print('# iqtree', file=sys.stdout)
		#print >> sys.stderr, treefiles, d_gene_count
		self.print_topology(treefiles, d_gene_count=d_gene_count)
		print('# astral', file=sys.stdout)
		self.print_topology(treefiles2, d_gene_count=d_gene_count2)

		self.clean(self.tmpdir)

		# chromosome tree
		max_trees = 30
		i = 0
		cmd_list = []
		treefiles = []
		treefiles2 = []
		d_gene_count = {}
		d_gene_count2 = {}
		print(len(d_chromfiles), 'chromosome groups', file=sys.stdout)
		# 按染色体串联建树，允许部分基因丢失
		for chroms, alnfiles in sorted(list(d_chromfiles.items()), key=lambda x: -len(x[1])):
			ngene = len(alnfiles)
			if len(chroms) > len(set(chroms)):
				continue
			i += 1
			if i > max_trees:
				continue
			alnfiles = self.chrom_tree(chroms)
			print(len(self.d_chroms), list(self.d_chroms.items())[:10], file=sys.stderr)
			prefix = 'CHR_' + '-'.join(chroms) + '_'+str(ngene) + '_' + str(len(alnfiles))
			cmds = self.concat_tree(alnfiles, prefix, idmap=self.d_chroms, astral=True)
			treefile = self.iqtree_treefile
			treefiles += [treefile]
			d_gene_count[treefile] = len(alnfiles)
			# astral
			treefile = self.astral_treefile
			treefiles2 += [treefile]
			d_gene_count2[treefile] = len(alnfiles)
			
			cmd_list += [cmds]
			#print >> sys.stderr, 'dot', chroms
			cmd_list += self.dot_plot(chroms)
			#print >> sys.stderr, prefix, len(alnfiles)
		cmd_file = '{}/chrom-cmds.list'.format(self.tmpdir)
		if cmd_list:
			run_job(cmd_file, cmd_list=cmd_list, tc_tasks=50)
		print(i, 'chromosome groups', file=sys.stdout)
		print('# iqtree', file=sys.stdout)
		self.print_topology(treefiles, d_gene_count=d_gene_count)
		print('# astral', file=sys.stdout)
		self.print_topology(treefiles2, d_gene_count=d_gene_count2)
		
		# clean
		self.clean(self.tmpdir)
		return	# 终止
		# phase, 不太成功
		phased_chroms, d_rename = self.phase_trees(xxchroms)
		alnfiles = [d_concat_alnfiles[chroms] for chroms in phased_chroms]
		genes = sum([d_gene_count3[alnfile] for alnfile in alnfiles])
		print('{} groups in {} clusters phased'.format(genes, len(alnfiles)), file=sys.stderr) 
		concat_alnfile = '{}/{}.aln'.format(self.tmpdir, 'phased')
		print(alnfiles, file=sys.stderr)
		with open(concat_alnfile, 'w') as fout:
			catAln(alnfiles, fout, idmap=d_rename)
		root = None
		for rc in SeqIO.parse(concat_alnfile, 'fasta'):
			sp, chrom = rc.id.split('-', 1)
			if sp == self.root:
				root = rc.id
				break
		cmds = []
		iqtreefile = concat_alnfile + '.treefile'
		treefile = rooted_treefile = concat_alnfile + '.tre'
		opts = ''
		if not root is None:
			opts = '-o ' + root
		if not os.path.exists(iqtreefile):
			cmd = '. ~/.bashrc; iqtree -redo -s {} -nt AUTO -bb 1000 {} &> /dev/null'.format(concat_alnfile, opts)
			cmds += [cmd]
		if root is None:
			cmd = 'nw_reroot {} '.format(iqtreefile)
		else:
			cmd = 'nw_reroot {intre} {root}'.format(intre=iqtreefile, root=root,)
		cmd += ' | nw_topology -I - | nw_order - | nw_order - -c d > {}'.format(rooted_treefile)
		if not os.path.exists(iqtreefile):
			cmds += [cmd]
		cmds = ' && '.join(cmds)
		run_cmd(cmds, log=True)

	def dot_plot(self, chroms):
		xchroms = self.groupby_species(chroms)
		cmds = []
		for chroms1, chroms2 in itertools.combinations_with_replacement(xchroms, 2):
			#print >> sys.stderr, chroms1, chroms2
			prefix = 'dotplot.{}-{}'.format('_'.join(chroms1), '_'.join(chroms2))
			ctl = prefix + '.ctl'
#			with open(ctl, 'w') as fout:
#				print >> fout, '1500\n1500\n{}\n{}'.format(','.join(chroms1), ','.join(chroms2))
			cmd = 'python /share/home/nature/src/dot_plotter.py -s pairs.collinearity -g pairs.gff -c {} \
				--kaks kaks.homology.kaks --ks-hist --max-ks 3 -o {} --plot-ploidy'.format(ctl, prefix)
			#run_cmd(cmd)
			cmds += [cmd]
		return cmds
	def clean(self, tmpdir):
		suffixes = [ 'fa', 'aln', #'trimal',  # 'aln'
				'bionj', 'contree', 'ckp.gz', 'iqtree', 'log', 'mldist', 'model.gz', 'splits.nex', 'uniqueseq.phy'
					]
		for suffix in suffixes:
			cmd = 'rm {}/*.{}'.format(tmpdir, suffix)
			run_cmd(cmd)
		prefixes = ['gene', ]
		for prefix in prefixes:
			cmd = 'rm {}/{}*'.format(tmpdir, prefix)
			run_cmd(cmd)

	def concat_tree(self, alnfiles, prefix, idmap=None, astral=False):
		concat_alnfile = '{}/{}.concat'.format(self.tmpdir, prefix)
		with open(concat_alnfile, 'w') as fout:
			catAln(alnfiles, fout, idmap=idmap)
		root = None
		for rc in SeqIO.parse(concat_alnfile, 'fasta'):
			sp, chrom = rc.id.split('-', 1)
			if sp == self.root:
				root = rc.id
				break
		cmds = []
		iqtreefile = concat_alnfile + '.treefile'
		opts = ''
		if not root is None:
			opts = '-o ' + root
		if True: #not os.path.exists(iqtreefile):
			cmd = 'iqtree -redo -s {} -nt AUTO -bb 1000 {} -mset JTT &> /dev/null'.format(concat_alnfile, opts)
			cmds += [cmd]
		self.iqtree_treefile = treefile = rooted_treefile = concat_alnfile + '.tre'
		
		if root is None:
			cmd = 'nw_reroot {} '.format(iqtreefile)
		else:
			cmd = 'nw_reroot {intre} {root} | nw_prune - {root}'.format(
				intre=iqtreefile, root=root,)
		cmd += ' | nw_topology -I - | nw_order - | nw_order - -c d | nw_order - > {}'.format(rooted_treefile)
		if True: #not os.path.exists(iqtreefile):
			cmds += [cmd]
		
		# astral
		if astral:
			iqtreefiles = [alnfile + '.treefile' for alnfile in alnfiles]
			genetrees = '{}/{}.genetrees'.format(self.tmpdir, prefix)
			self.cat_genetrees(iqtreefiles, genetrees, idmap=self.d_chroms, plain=False)
			sptree = genetrees + '.astral'
			cmd = '''mem=50g
	astral_root=~/bin/Astral-MP-5.14.5
	java -Xmx50g -D"java.library.path=$astral_root/lib" -jar $astral_root/astral.*.jar -i {} -o {}'''.format(
				genetrees, sptree)
			cmds += [cmd]
			if root is None:
				cmd = 'nw_reroot {} '.format(sptree)
			else:
				cmd = 'nw_reroot {intre} {root} | nw_prune - {root}'.format(
					intre=sptree, root=root,)
			self.astral_treefile = treefile = rooted_treefile = sptree + '.nwk'
			
			cmd += ' | nw_topology -I - | nw_order - | nw_order - -c d | nw_order - > {}'.format(rooted_treefile)
			cmds += [cmd]
		cmds = ' && '.join(cmds)
		return cmds
	def cat_genetrees(self, treefiles, genetrees, idmap=None, **kargs):
		with open(genetrees, 'w') as fout:
			for iqtreefile in treefiles:
				if not os.path.exists(iqtreefile):
					logger.warn('{} not exists'.format(iqtreefile))
					continue
				newick = self.get_topology(iqtreefile, idmap=idmap, **kargs)
				print(newick, file=fout)
	def phase_trees(self, xxchroms):
		'''[('Cc7', 'Cs3', 'Cs8', 'Ns1', 'Ns2'), ('Cc2', 'Cs3', 'Cs9', 'Ns14', 'Ns18')]'''
		array = []
		for chroms in xxchroms:
			chroms = self.groupby_species(chroms)
			array += [chroms]
		array = np.array(array)
		d_phased = {}
		for i in range(array.shape[1]):
			xchroms = array[:, i]
			phased = self.phase_chroms(xchroms)
			d_phased.update(phased)

		phased_chroms = []
		d_rename = {}	# idmap
		
		for xchroms in array:
			if not all([chroms in d_phased for chroms in xchroms]):	# all phased
				continue
			flatten = tuple()
			phased2 = []
			for chroms in xchroms:
				flatten += chroms
				phased = d_phased[chroms]
				phased2 += [phased]
				for i, chrom in enumerate(phased):
					sp = self.d_species[chrom]
					new_name = '{}-{}'.format(sp, i)
					chrom_id = '{}-{}'.format(sp, chrom)
					self.d_species[new_name] = sp
					d_rename[chrom_id] = new_name
			print(xchroms, '->', phased2, file=sys.stderr)
			phased_chroms += [flatten]
		return phased_chroms, d_rename
	def groupby_species(self, chroms):
		'''('Cc7', 'Cs3', 'Cs8', 'Ns1', 'Ns2')'''
		array = [(chrom, self.d_species[chrom]) for chrom in chroms]
		chroms = []
		for sp, item in itertools.groupby(array, key=lambda x:x[1]):
			chroms.append(tuple([chrom for chrom,sp in item]))
		return chroms	# [['Cc7'], ['Cs3', 'Cs8'], ['Ns1', 'Ns2']]
	def phase_chroms(self, xchroms):
		'''chroms = [('Cs3', 'Cs8'), ('Cs3', 'Cs9')]'''	# sorted
		xchroms = list(map(tuple, xchroms))
		length = len(xchroms)
		d_index = {}
		d_phased = {}
		chroms = xchroms[0]
		ploidy = len(chroms)
		if ploidy == 1:		# no need to phase
			for chroms in xchroms:
				d_phased[chroms] = chroms
			return d_phased
		d_phased[chroms] = chroms  # init
		index = list(range(ploidy))
		for i, chrom in enumerate(chroms):
			d_index[chrom] = i
		xchroms.pop(0)
		while True:
			pops = []
			for j, chroms in enumerate(xchroms):
				has_index = [d_index[chrom] for chrom in chroms if chrom in d_index]
				if len(has_index) > len(set(has_index)):	# conflict
					continue
				elif len(chroms) - len(has_index) == 1:		# phasable
					chrom_wo_index = [chrom for chrom in chroms if chrom not in d_index][0]
					to_index = list(set(index) - set(has_index))[0]
					d_index[chrom_wo_index] = to_index
					array = [(d_index[chrom], chrom) for chrom in chroms]
					phased = [chrom for idx, chrom in sorted(array)]
					d_phased[chroms] = tuple(phased)
					print(chroms, '->', d_phased[chroms], file=sys.stderr)
					pops += [j]
				elif len(chroms) == len(has_index): 	# phased
					pops += [j]
			if len(pops) == 0:
				break
			for idx in reversed(pops):
				xchroms.pop(idx)
		
		print('{} / {} phased for {}..'.format(len(d_phased), length, chroms[0]), file=sys.stderr)
		return d_phased
	def count_topology(self, treefiles, d_gene_count={}):
		d_top_count = {}
		d_top_count2 = {}
		for treefile in treefiles:
			if not os.path.exists(treefile) or os.path.getsize(treefile) == 0:
				logger.warn('{} not exists'.format(treefile)) 
				continue
			try: topology = self.get_topology(treefile)
			except ValueError as e:
				logger.warn('{}: {}'.format(treefile, e))
				continue
			gene_count = d_gene_count.get(treefile, 1)
			try: d_top_count[topology] += gene_count
			except KeyError: d_top_count[topology] = gene_count
			try: d_top_count2[topology] += 1
			except KeyError: d_top_count2[topology] = 1
		return d_top_count, d_top_count2
	def print_topology(self, treefiles, **kargs):
		d_top_count, d_top_count2 = self.count_topology(treefiles, **kargs)
		for top, count in sorted(list(d_top_count.items()), key=lambda x:-x[1]):
			print(top, count, d_top_count2[top], file=sys.stdout)
	def get_min_bootstrap(self, treefile):
		tree = Phylo.read(treefile, 'newick')
		bootraps = [clade.confidence for clade in tree.get_nonterminals() if clade.confidence>=0]
		#print treefile, tree, bootraps
		min_bs = min(bootraps) if bootraps else 0
		return min_bs
	def get_topology(self, treefile, idmap=None, plain=True, **kargs):
		# to_strings(self, confidence_as_branch_length=False, branch_length_only=False, plain=False, plain_newick=True, ladderize=None, max_confidence=1.0, format_confidence='%1.2f', format_branch_length='%1.5f')
		from Bio.Phylo.NewickIO import Writer
		if idmap is None:
			try: idmap = self.d_species
			except AttributeError: idmap = {}
		tree = Phylo.read(treefile, 'newick')
		if idmap:
			for clade in tree.get_terminals():
				if idmap and clade.name not in idmap:
					logger.warn('ID `{}` in {} is not exists in idmap'.format(clade.name, treefile))
				clade.name = idmap.get(clade.name, clade.name)
		newick = list(Writer([tree]).to_strings(plain=plain, **kargs))[0]
		return newick
	def get_root(self):
		for sp, ploidy in list(self.sp_dict.items()):
			if ploidy == 1:
				return sp
			return sp
	def get_min_ks(self):
		self.d_gff = d_gff = Gff(self.gff).get_genes()
		d_matrix = {}
		keys = list(self.sp_dict.keys())
		for sp1, sp2 in itertools.product(keys, keys):
			d_matrix[(sp1, sp2)] = 0
		i = 0
		for sg in self.subgraphs():
			i += 1
			#print >> sys.stderr, '\t'.join(sorted(sg.nodes()))
			for gene in sg.nodes():
				d_sg_ks = {}
				for neighbor in sg.neighbors(gene):
					key = tuple(sorted([gene, neighbor]))
					ks = self.d_ks[key]
					d_sg_ks[(gene, neighbor)] = ks
				min_pair = min(d_sg_ks, key=lambda x: d_sg_ks[x])
		#		print >> sys.stderr, gene, d_sg_ks, min_pair
				sp_pair = tuple(genes2species(min_pair))
				d_matrix[sp_pair] += 1
		print(i, 'groups', file=sys.stderr)
		print(d_matrix)
def orthomcl_to_astral(source='orthomcl', **kargs):
	ToAstral(source=source, **kargs).run()
def orthomcl_stats(source='orthomcl', **kargs):
	ToAstral(source=source, **kargs).stat()

class ToAstral(ColinearGroups):
	def __init__(self, input=None, pep=None, spsd=None, cds=None, tmpdir='tmp', root=None, both=True, suffix=None, 
			ncpu=50, max_taxa_missing=0.5, max_mean_copies=10, max_copies=5, singlecopy=False, onlyaln=False,
			source=None, orthtype='orthologues', fast=True, concat=False, clean=False, overwrite=False, 
			trimal_opts='-automated1', iqtree_opts=''):
		self.input = input
		self.pep = pep
		self.cds = cds
		self.spsd = spsd
		self.root = root.split() if isinstance(root, str) else root
		self.both = both
		self.ncpu = ncpu
		self.tmpdir = tmpdir
		self.max_taxa_missing = max_taxa_missing
		self.max_mean_copies = max_mean_copies
		self.max_copies = max_copies
		if singlecopy:
			self.max_copies = 1
		self.singlecopy = singlecopy
		self.concat = concat
		self.fast = fast
		self.sp_dict = parse_spsd(spsd, skip=False)
		self.suffix = input if suffix is None else suffix 
		self.orthtype = orthtype
		self.source = source
		self.trimal_opts = trimal_opts
		self.iqtree_opts = iqtree_opts
		self.clean = clean
		self.onlyaln = onlyaln
		self.overwrite = overwrite
	def lazy_get_groups(self, orthtype='Orthogroups'):
		species = parse_species(self.spsd) #subset
		if os.path.isdir(self.input): # orthofinder
			source = 'orthofinder' + '-' + orthtype.lower()
			result = OrthoFinder(self.input)
			if species is None:
				species = result.Species
			if orthtype.lower() == 'orthogroups':
				groups = result.get_orthogroups(sps=species)
			elif orthtype.lower() == 'orthologues':
				groups = result.get_orthologs_cluster(sps=species)
			else:
				raise ValueError("Unknown type: {}. MUST in ('orthogroups', 'orthologues')".format(orthtype))
		elif self.source is not None and  self.source.lower() == 'orthomcl':
			source = 'orthomcl' + '-' + orthtype.lower()
			if species is None:
				species = OrthoMCLGroup(self.input).get_species()
			groups = OrthoMCLGroup(self.input, sps=species)
			
		else:
			source = 'mcscanx'
			if species is None:
				species = Collinearity(self.input).get_species()
			result = ColinearGroups(self.input, spsd=self.spsd)
			groups = result.groups
		self.species = species
		self.source = source
		return groups
	def stat(self, nbin=20):
		ng, rg = 0, 0
		d_taxon = {}
		d_gene = {}
		for og in self.lazy_get_groups(orthtype=self.orthtype):
			rg += 1
			got_sp = [(sp, genes) for sp, genes in list(og.spdict.items()) if len(genes) <= self.sp_dict.get(sp, self.max_copies) ]
			taxa_occupancy = 1e2*len(got_sp) / len(self.species)
			d_taxon[og.ogid] = taxa_occupancy
			taxa_missing = 1 - taxa_occupancy/100
			if taxa_missing > self.max_taxa_missing:
				continue
			ng += 1
			for sp, genes in got_sp:
				ngs = len(genes)
				try: d_gene[sp] += [ngs]
				except KeyError: d_gene[sp] = [ngs]
		# 
		xs = 'sc' if self.singlecopy else 'mc'
		self.suffix = '{}.{}'.format(self.suffix, xs)
		logger.info('{} taxa; {} -> {} genes'.format(len(self.species), rg, ng))

		# taxa missing
		step = 100//nbin
		d_bin = {}
			
		for val in d_taxon.values():
			bin = int(val) // step * step
			try: d_bin[bin] += 1
			except KeyError: d_bin[bin] = 1
		to_file = self.suffix + '.taxa_occupancy'
		logger.info('outputing global taxa occupancy to `{}`'.format(to_file))
		f = open(to_file, 'w')
		line = ['taxa_occupancy', 'gene_count']
		f.write('\t'.join(map(str, line)) + '\n')
		for bin, count in sorted(d_bin.items()):
			line = [bin, count]
			f.write('\t'.join(map(str, line)) + '\n')
		f.close()

		# gene occupancy
		self.suffix = '{}.mm{}'.format(self.suffix, self.max_taxa_missing)
		go_file = self.suffix + '.gene_occupancy'
		logger.info('outputing gene occupancy to `{}`'.format(go_file))
		f = open(go_file, 'w')
		tiles = [5,25,75,90, 95]
		line = ['species', 'gene_count', 'gene_occupancy%', 'median_copy', 'mean_copy', 
			] + ['tile{}'.format(tile) for tile in tiles]
		f.write('\t'.join(map(str, line)) + '\n')
		for sp in sorted(self.species):
			xg = d_gene[sp]
			xng = len(xg)
			line = [sp, xng, round(1e2*xng/ng,2), np.median(xg), round(np.mean(xg), 2),
					] + [ np.percentile(xg, tile) for tile in tiles]
			f.write('\t'.join(map(str, line)) + '\n')
		f.close()
		
	def run(self):
		#logger.info('VARS: {}'.format(self.__dict__))
		mafft_template = '. ~/.bashrc; mafft --auto {} > {} 2> /dev/null'
		pal2nal_template = 'pal2nal.pl -output fasta {} {} > {}'
		trimal_template = 'trimal %s -in {} -out {} > /dev/null' % (self.trimal_opts, )
		iqtree_template = 'iqtree2 -redo -s {} %s -nt 1 {} > /dev/null' % (self.iqtree_opts, )
		reroot_template = 'mv {tree} {tree}.bk && nw_reroot -l {tree}.bk {root} | nw_order -c n - > {tree}'
		mkdirs(self.tmpdir)
		d_pep = seq2dict(self.pep)
		d_cds = seq2dict(self.cds) if self.cds else {}
		d_idmap = {}
		pepTreefiles, cdsTreefiles = [], []
		pepAlnfiles, cdsAlnfiles = [], []
		cmd_list = []
		roots = []
		i,j = 0, 0
		for og in self.lazy_get_groups(orthtype=self.orthtype):
			#species = og.species
			#nsp = len(set(species))
			# compatible with single-copy, low-copy, and limited-copy
			got_sp = [(sp, genes) for sp, genes in list(og.spdict.items()) if len(genes) <= self.sp_dict.get(sp, self.max_copies) ]
			taxa_missing = 1 - 1.0*len(got_sp) / len(self.species)
			
			if taxa_missing > self.max_taxa_missing:
				continue
			if og.mean_copies > self.max_mean_copies:
				continue
			random.shuffle(got_sp)

			iters = []
			for (sp, genes) in got_sp:
				for g in genes:
					iters += [(g, sp)]

			if len(iters) < 4:
				continue
			i += 1
			#print [sp for (sp, genes) in got_sp]

			ogid = og.ogid
			pepSeq = '{}/{}.pep'.format(self.tmpdir, ogid)
			cdsSeq = '{}/{}.cds'.format(self.tmpdir, ogid)
			f_pep = open(pepSeq, 'w')
			f_cds = open(cdsSeq, 'w')
			d_root = {}
			for gene, sp in iters:
				try: rc = d_pep[gene]
				except KeyError:
					logger.warn('{} not found in {}; skipped'.format(gene, self.pep))
					continue
				rc.id = format_id_for_iqtree(gene)
				d_idmap[rc.id] = sp
				SeqIO.write(rc, f_pep, 'fasta')
				if self.cds:
					rc = d_cds[gene]
					rc.id = format_id_for_iqtree(gene)
					SeqIO.write(rc, f_cds, 'fasta')
				if self.root and sp in set(self.root):
					d_root[sp] = rc.id
			f_pep.close()
			f_cds.close()
			
			if d_root:
				j += 1
			root = ' '.join(list(d_root.values()))
			pepAln = pepSeq + '.aln'
			cdsAln = cdsSeq + '.aln'
			pepTrim = pepAln + '.trimal'
			cdsTrim = cdsAln + '.trimal'
			pepTreefile = pepTrim + '.treefile'
			cdsTreefile = cdsTrim + '.treefile'
			treefile = cdsTreefile if self.cds and not self.both else pepTreefile
			cmd = '[ ! -s {} ]'.format(treefile) if not self.overwrite else '[ true ]'
			cmds = [cmd]
			cmd = mafft_template.format(pepSeq, pepAln)
			cmds += [cmd]
			iqtree_opts0 = '' #' -o {} '.format(root) if root else ''
			pep = True
			if self.cds:
				iqtree_opts = iqtree_opts0 + ' -mset GTR ' if self.fast else iqtree_opts0 
				cmd = pal2nal_template.format(pepAln, cdsSeq, cdsAln)
				cmds += [cmd]
				cmd = trimal_template.format(cdsAln, cdsTrim)
				cmds += [cmd]
				cmd = iqtree_template.format(cdsTrim, iqtree_opts)
				cmds += [cmd]
				if root:
					cmd = reroot_template.format(tree=cdsTreefile, root=root)
					cmds += [cmd]
				cdsTreefiles += [cdsTreefile]
				cdsAlnfiles += [cdsTrim]
				pep = True if self.both else False
			if pep and not self.onlyaln:
				iqtree_opts = iqtree_opts0 + ' -mset JTT ' if self.fast else iqtree_opts0
				cmd = trimal_template.format(pepAln, pepTrim)
				cmds += [cmd]
				cmd = iqtree_template.format(pepTrim, iqtree_opts)
				cmds += [cmd]
				if root:
					cmd = reroot_template.format(tree=pepTreefile, root=root)
					cmds += [cmd]
				pepTreefiles += [pepTreefile]
				pepAlnfiles += [pepTrim]
			roots += [root]
			cmds = ' && '.join(cmds)
			cmd_list += [cmds]
		logger.info('total {} groups, {} rooted'.format(i,j))
		pepTreefiles = [t for _, t in sorted(zip(roots, pepTreefiles), reverse=1)]	# prefer to rooted
#		with open('/tmp/sort.tree', 'w') as f:
#			print(pepTreefiles, file=f)
		cdsTreefiles = [t for _, t in sorted(zip(roots, cdsTreefiles), reverse=1)]
		if self.suffix is None:
			self.suffix = '{}_to_astral'.format(self.source)
		xs = 'sc' if self.singlecopy else 'mc'
		self.suffix = '{}.{}'.format(self.suffix, xs)

		nbin = 10
		cmd_file = '{}/{}.cmds.list'.format(self.tmpdir, self.suffix)
		run_job(cmd_file, cmd_list=cmd_list, tc_tasks=self.ncpu, by_bin=nbin, fail_exit=False)

		# cat genetrees
		pepGenetrees = '{}.pep.mm{}.genetrees'.format(self.suffix, self.max_taxa_missing)
		cdsGenetrees = '{}.cds.mm{}.genetrees'.format(self.suffix, self.max_taxa_missing)
		if not self.onlyaln:
			for treefiles, genetrees in zip([pepTreefiles, cdsTreefiles], [pepGenetrees, cdsGenetrees]):
				logger.info('combining {} gene trees into `{}`'.format(len(treefiles), genetrees))
				self.cat_genetrees(treefiles, genetrees, idmap=d_idmap, plain=False, format_confidence='%d')

		# concat alignments
		cdsCatAln = '{}.cds.mm{}.concat.aln'.format(self.suffix, self.max_taxa_missing)
		pepCatAln = '{}.pep.mm{}.concat.aln'.format(self.suffix, self.max_taxa_missing)
		if self.singlecopy and self.concat:
			for alnfiles, _catAln in zip([pepAlnfiles, cdsAlnfiles], [pepCatAln, cdsCatAln]):
				logger.info('concatenating {} alignments into `{}`'.format(len(alnfiles), _catAln))
				with open(_catAln, 'w') as outAln:
					catAln(alnfiles, outAln, idmap=d_idmap)
		# clean
		if self.cds and not self.both:
			rmdirs(pepGenetrees, pepCatAln)
		elif not self.cds:
			rmdirs(cdsGenetrees, cdsCatAln)

		if self.clean:
			logger.info('cleaning `{}`'.format(self.tmpdir))
			rmdirs(self.tmpdir)

def parse_spsd(spsd, skip=False):
	d = OrderedDict()
	if spsd is None:
		return d
	for line in open(spsd):
		temp = line.strip().split()
		if not temp:
			continue
		try:
			sp, ploidy = temp[:2]
		except ValueError:
			if skip:
				continue
			sp = temp[0]
			ploidy = 1
		d[sp] = int(ploidy)
	return d
def get_chrs(collinearity):
	d = {}
	for rc in Collinearity(collinearity):
		chr1, chr2 = rc.chrs
		for g1, g2 in rc.pairs:
			d[g1] = chr1
			d[g2] = chr2
	return d
def get_pair(collinearity, minN=0):
	for rc in Collinearity(collinearity):
		if rc.N < minN:
			continue
		for g1, g2 in rc.pairs:
			yield g1, g2
def gene2species(gene, sep="|"):
	return gene.split(sep)[0]
def genes2species(genes, sep="|"):
	return [gene2species(gene, sep) for gene in genes]

def seq2dict(seq):
	from Bio import SeqIO
	return dict([(rc.id, rc)for rc in SeqIO.parse(seq, 'fasta')])

def test():
	collinearity, gff, chrmap = sys.argv[1:4]
	outTab = sys.stdout
	blocks = Collinearity(collinearity, gff, chrmap)
	for rc in blocks: #.parse():
		line = [rc.Alignment, rc.chr1, rc.start1, rc.end1, rc.chr2, rc.start2, rc.end2]
		line = list(map(str, line))
		print('\t'.join(line), file=outTab)
		for gene in rc.genes1:
			print(gene.info)
		
def list_blocks(collinearity, outTsv, gff=None, kaks=None):
	'''以共线性块为单位，输出信息'''
	line = ["Alignment", "species1", "species2", "chr1", "chr2", "start1", "end1", "length1", "start2", "end2", "length2", "strand", "N_gene", "mean_Ks", 'median_Ks', 'score', 'e_value']
	line = ["id", "species1", "species2", "chr1", "chr2", "start1", "end1", "length1", "start2", "end2", "length2", "strand", "length", "ks_average", 'ks_median', 'score', 'e_value']
	print('\t'.join(line), file=outTsv)
	for rc in Collinearity(collinearity,gff=gff,kaks=kaks):
		sp1, sp2 = rc.species
		chr1, chr2 = rc.chrs
		Alignment, score, e_value, N, strand = rc.Alignment, rc.score, rc.e_value, rc.N, rc.strand
		start1, end1, length1 = rc.start1, rc.end1, rc.length1
		start2, end2, length2 = rc.start2, rc.end2, rc.length2
		mean_ks = rc.mean_ks
		median_ks = rc.median_ks
		line = [Alignment, sp1, sp2, chr1, chr2, start1, end1, length1, start2, end2, length2, strand, N, mean_ks, median_ks, score, e_value]
		line = list(map(str, line))
		print('\t'.join(line), file=outTsv)
def gene_class(collinearity, inTsv, outTsv, byAlignment=True):
	'''将共线性块的分类信息传递到基因对'''
	d_info = {}
	for line in open(inTsv):
		temp = line.strip().split('\t')
		Alignment, gClass = temp[0], temp[-1]
		d_info[Alignment] = gClass
	for rc in Collinearity(collinearity):
		Alignment = rc.Alignment
		if Alignment not in d_info:
			continue
		for g1, g2 in rc.pairs:
			line = [g1, g2, d_info[Alignment]]
			print('\t'.join(line), file=outTsv)
def list_pairs(collinearity, outTsv, gff=None, kaks=None, blocks=None):
	'''提取block内基因对的信息'''
	line = ['gene1', 'gene2', 'Ks', "chr1", "start1", "end1", "strand1", "chr2", "start2", "end2", "strand2", "Alignment"] # + ["Alignment", "species1", "species2", "chr1", "chr2", "start1", "end1", "length1", "start2", "end2", "length2", "strand", "N_gene", "mean_Ks", 'median_Ks', 'score', 'e_value']
	print('\t'.join(line), file=outTsv)
	if blocks is not None:
		d_blocks = {}
		for line in open(blocks):
			temp = line.strip().split('\t')
			Alignment = temp[0]
			d_blocks[Alignment] = None
	for rc in Collinearity(collinearity,gff=gff,kaks=kaks):
		sp1, sp2 = rc.species
		chr1, chr2 = rc.chrs
		Alignment, score, e_value, N, strand = rc.Alignment, rc.score, rc.e_value, rc.N, rc.strand
		if blocks is not None and not Alignment in d_blocks:
			continue
		start1, end1, length1 = rc.start1, rc.end1, rc.length1
		start2, end2, length2 = rc.start2, rc.end2, rc.length2
		mean_ks = rc.mean_ks
		median_ks = rc.median_ks
		line0 = [Alignment, sp1, sp2, chr1, chr2, start1, end1, length1, start2, end2, length2, strand, N, mean_ks, median_ks, score, e_value]
		for g1, g2, ks in zip(rc.genes1, rc.genes2, rc.ks):
			line = [g1.id, g2.id, ks, g1.chr, g1.start, g1.end, g1.strand,  g2.chr, g2.start, g2.end, g2.strand, Alignment] #+ line0
			line = list(map(str, line))
			print('\t'.join(line), file=outTsv)
def block_ks(collinearity, kaks, outkaks, min_n=10):
	for rc in Collinearity(collinearity, kaks=kaks):
		if rc.N < min_n:
			continue
		for pair in rc.pairs:
			try:
				info = rc.d_kaks[pair]
			except KeyError:
				continue
			info.write(outkaks)
def bin_ks_by_chrom(collinearity, gff, kaks, sp1, sp2, out=sys.stdout, bin_size=500000):
	lines = []
	for rc in Collinearity(collinearity,gff=gff,kaks=kaks):
		if not rc.is_sp_pair(sp1, sp2):
			continue
		chr1, chr2 = list(map(get_chrom, [rc.chr1, rc.chr2]))
		same_order = rc.is_sp_pair(sp1, sp2) == (sp1, sp2)
		for g1, g2, ks in zip(rc.genes1, rc.genes2, rc.ks):
			g = g1 if same_order else g2
			g.ks = ks
			g.bin = g.start // bin_size
			lines += [g]
	lines = sorted(lines, key=lambda x: (x.chr, x.start))
	bin = 20
	for chrom, genes in itertools.groupby(lines, key=lambda x: x.chr):
		genes = list(genes)
		for i in range(0, len(genes), bin):
			gs = genes[i:i+bin]
#		for BIN, gs in itertools.groupby(genes, key=lambda x: x.bin):
			gs = list(gs)
			if len(gs) < 10:
				continue
			median_ks = np.median([g.ks for g in gs])
#			start = BIN*bin_size
#			end = start + bin_size
#			line = [chrom, start, end, median_ks]
			line = [chrom, gs[0].start, gs[-1].end, median_ks]
			line = list(map(str, line))
			print('\t'.join(line), file=out)

def count_genes(collinearity, sp1, sp2):
	d_count = {}
	for rc in Collinearity(collinearity):
		if not rc.is_sp_pair(sp1, sp2):
			continue
		chr1, chr2 = list(map(get_chrom, [rc.chr1, rc.chr2]))
		if rc.is_sp_pair(sp1, sp2) != (sp1, sp2):
			chr1, chr2 = chr2, chr1
		ngene = rc.N
		try: d_count[chr2] += ngene
		except KeyError: d_count[chr2] = ngene
	for chrom, ngene in sorted(d_count.items()):
		print(chrom, ngene)
def main():
	import sys
	subcmd = sys.argv[1]
	kargs = parse_kargs(sys.argv)
	if subcmd == 'list_blocks':	# 列出所有共线性块
		list_blocks(collinearity=sys.argv[2], outTsv=sys.stdout, gff=sys.argv[3], kaks=sys.argv[4])
	elif subcmd == 'gene_class': # 按共线性块的分类对基因对进行分类
		gene_class(collinearity=sys.argv[2], inTsv=sys.argv[3], outTsv=sys.stdout)
	elif subcmd == 'list_pairs':	# 列出指定共线性块的基因对
		list_pairs(collinearity=sys.argv[2], outTsv=sys.stdout, gff=sys.argv[3], kaks=sys.argv[4], blocks=sys.argv[5])
	elif subcmd == 'get_gff':	# 获取指定物种集的gff
		gff = sys.argv[2]
		species = sys.argv[3]
		fout = sys.stdout
		get_gff(gff, species, fout)
	elif subcmd == 'slim_tandem':	# 串联重复簇只保留一个基因，其他从基因对中剔除，用于共线性分析
		tandem, pairs = sys.argv[2:4]
		outPairs = sys.stdout
		slim_tandem(tandem, pairs, outPairs)
	elif subcmd == 'test_closest':	# 获取Ks最小的物种对
		collinearity, kaks = sys.argv[2:4]
		spsd = sys.argv[4]
		test_closest(collinearity, kaks, spsd)
	elif subcmd == 'cg_trees':	# 按照物种倍性构建基因树/染色体树，以单个基因作为anchor
		collinearity, spsd, seqfile, gff = sys.argv[2:6]
		try: tmpdir = sys.argv[6]
		except IndexError: tmpdir = 'tmp'
		cg_trees(collinearity, spsd, seqfile, gff, tmpdir, **kargs)
	elif subcmd == 'anchor_trees':	# 按照物种倍性构建基因树/染色体树，以染色体作为anchor（因此过滤标准严格）
		collinearity, spsd, seqfile, gff = sys.argv[2:6]
		try: tmpdir = sys.argv[6]
		except IndexError: tmpdir = 'tmp'
		anchor_trees(collinearity, spsd, seqfile, gff, tmpdir)
	elif subcmd == 'gene_trees':	# 按照物种倍性构建基因树
		collinearity, spsd, seqfile = sys.argv[2:5]
		try: orthologs = sys.argv[5]
		except IndexError: orthologs = None
		try: tmpdir = sys.argv[6]
		except IndexError: tmpdir = 'tmp'
		gene_trees(collinearity, spsd, seqfile, orthologs, tmpdir)
	elif subcmd == 'to_phylonet':	# 用于phylonet
		collinearity, spsd, seqfile, outprefix = sys.argv[2:6]
		to_phylonet(collinearity, spsd, seqfile, outprefix)
	elif subcmd == 'block_ks':  # 过滤掉block过短的ks
		collinearity, kaks = sys.argv[2:4]
		outkaks = sys.stdout
		try: min_n = int(sys.argv[4])
		except IndexError: min_n = 10
		block_ks(collinearity, kaks, outkaks, min_n=min_n)
		
	elif subcmd == 'count_genes':
		collinearity, sp1, sp2 = sys.argv[2:5]
		count_genes(collinearity, sp1, sp2)
	elif subcmd == 'to_ark':	# 转换ARK格式
		collinearity, spsd, gff = sys.argv[2:5]
		try: max_missing = float(sys.argv[5])
		except IndexError: max_missing = 0.2
		to_ark(collinearity, spsd, gff, max_missing=max_missing)
	elif subcmd == 'to_synet':	# synet tree
		collinearity, spsd = sys.argv[2:4]
		ColinearGroups(collinearity, spsd).to_synet(fout=sys.stdout)
	elif subcmd == 'cluster':	# cluster by infomap
		collinearity = sys.argv[2]
		ColinearGroups(collinearity).infomap()
	elif subcmd == 'block_length': # block_length distribution
		collinearity, sp_pairs = sys.argv[2:4]
		block_length(collinearity, sp_pairs)
	elif subcmd == 'gene_retention':	# 基因保留与丢失
		collinearity, spsd, gff = sys.argv[2:5]
		gene_retention(collinearity, spsd, gff)
	elif subcmd == 'anchors2bed':	# 提取一段block的坐标
		collinearity, gff, chrmap, left_anchors, right_anchors = sys.argv[2:7]
		anchors2bed(collinearity, gff, chrmap, left_anchors, right_anchors, outbed=sys.stdout)
	elif subcmd == 'bin_ks':
		collinearity, gff, kaks, sp1, sp2 = sys.argv[2:7]
		bin_ks_by_chrom(collinearity, gff, kaks, sp1, sp2)
	elif subcmd == 'to_astral':
		input, pep = sys.argv[2:4]
		ToAstral(input, pep, **kargs).run()
	elif subcmd == 'to_wgdi':
		try: gff, chrLst, pep, cds = sys.argv[2:6]
		except: gff, chrLst, pep, cds = 'all_species_gene.gff', 'chr.list', 'pep.faa', 'cds.fa'
		Gff(gff).to_wgdi(chrLst, pep, cds, **kargs)
	elif subcmd == 'cr':
		collinearity, chrmap = sys.argv[2:4]
		collinearity_ratio(collinearity, chrmap, outMat=sys.stdout, **kargs)
	elif subcmd == 'ortho_block':
		collinearity = sys.argv[2:-1]
		OFdir = sys.argv[-1]
#		collinearity, OFdir = sys.argv[2:4]
		fout = sys.stdout
		identify_orthologous_blocks(collinearity, [OFdir], fout, **kargs)
	elif subcmd == 'get_ks':
		ksfile = sys.argv[2]
		get_ks(ksfile, pairfile=sys.stdin, outks=sys.stdout, outpair=sys.stderr, **kargs)
	else:
		raise ValueError('Unknown sub command: {}'.format(subcmd))
if __name__ == '__main__':
	main()
