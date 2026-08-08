# coding: utf-8
"""
paralog: output paralogous gene pairs produced at each branch from HOGs,
         and index collinearity blocks by branch-paralog content.

Internal nodes: sibling child HOGs (same parent) -> genes of the same species
  across different child HOGs are paralogs on that branch.
Leaf species: genes of the same species within the same HOG -> paralogs.
"""
from collections import defaultdict, Counter
from .hog import HOG
from .mcscan import XCollinearity
from .RunCmdsMP import logger


# ---------------------------------------------------------------------------
#  paralog output (used directly or by indexer)
# ---------------------------------------------------------------------------

class Paralog:
	def __init__(self, ogfile=None, orthfiles=None, sptreefile=None,
				 nodes=None, species=None, prefix='paralog', **kargs):
		self.ogfile = ogfile
		self.orthfiles = orthfiles
		self.sptreefile = sptreefile
		self.nodes = nodes
		self.species = species
		self.prefix = prefix
		self.kargs = kargs

	def run(self):
		"""Load HOGs, write paralog TSV, return {branch: frozenset(pairs)}, root_name.

		With hog_tsv: load HOGs from an existing HOGs.tsv instead of rebuilding.
		With inparalog: read an existing --no-index paralog.tsv (first 3 cols).
		"""
		hog_tsv = self.kargs.get('hog_tsv')
		inparalog = self.kargs.get('inparalog')
		fpath = self.prefix + '.paralog.tsv'

		if inparalog:
			branch_pairs, count = self._load_pairs(inparalog)
			root = self._root_name()
			logger.info('Loaded {} paralog pairs from {}'.format(count, inparalog))
			return {b: frozenset(ps) for b, ps in branch_pairs.items()}, root

		hog = HOG(ogfile=self.ogfile, orthfiles=self.orthfiles,
				  sptreefile=self.sptreefile, **self.kargs)
		if hog_tsv:
			hog.from_tsv(hog_tsv)
		else:
			write_tsv = bool(self.kargs.get('write_hog'))
			if write_tsv:
				hog.outtsv = self.prefix + '.HOGs.tsv'
			hog.pipe(write_tsv=write_tsv)
			if write_tsv:
				logger.info('HOGs written to {}'.format(hog.outtsv))
		logger.info('Loaded {} HOGs'.format(len(hog.all_hogs)))

		branch_pairs, count = self._write_and_group(hog, fpath)
		logger.info('Output {} paralog pairs to {}'.format(count, fpath))
		root = hog.tree.name if hog.tree else self._root_name()
		return {b: frozenset(ps) for b, ps in branch_pairs.items()}, root

	def _root_name(self):
		"""Root node name from species tree, or 'Root' if no tree."""
		if not self.sptreefile:
			return 'Root'
		from .tree import number_nodes
		return number_nodes(self.sptreefile).name

	def _load_pairs(self, fpath):
		"""Read paralog pairs from an existing paralog.tsv (first 3 columns).

		Applies the same --nodes/--species filters as HOG rebuilding.
		Return ({branch: set of canonical (g1,g2)}, count).
		"""
		node_set = set(self.nodes) if self.nodes else None
		sp_set = set(self.species) if self.species else None
		branch_pairs = defaultdict(set)
		count = 0
		with open(fpath) as f:
			for line in f:
				if line.startswith('#'):
					continue
				parts = line.rstrip().split('	')
				if len(parts) < 3:
					continue
				g1, g2, node = parts[0], parts[1], parts[2]
				if node_set and node not in node_set:
					continue
				if sp_set:
					sp = g1.split('|')[0] if '|' in g1 else g1.rsplit('_', 1)[0]
					if sp not in sp_set:
						continue
				pair = (g1, g2) if g1 < g2 else (g2, g1)
				branch_pairs[node].add(pair)
				count += 1
		return branch_pairs, count

	def _write_and_group(self, hog, fpath):
		"""Write paralog TSV, return {branch: set of canonical (g1,g2)}, count."""
		branch_pairs = defaultdict(set)
		with open(fpath, 'w') as fout:
			fout.write('#gene1\tgene2\tnode\tspecies\tHOG_id\n')
			count = 0
			for g1, g2, node_id, sp, hog_id in hog.iter_branch_paralogs(
					self.nodes, self.species):
				fout.write('{}\t{}\t{}\t{}\t{}\n'.format(
					g1, g2, node_id, sp, hog_id))
				pair = (g1, g2) if g1 < g2 else (g2, g1)
				branch_pairs[node_id].add(pair)
				count += 1
		return branch_pairs, count


# ---------------------------------------------------------------------------
#  branch-paralog indexer (default mode)
# ---------------------------------------------------------------------------

class ParalogIndexer:
	"""Assign collinearity blocks to branches based on paralog content."""

	def __init__(self, ogfile, orthfiles, sptreefile,
				 self_synteny=None, min_n=0, gff=None, min_dist=None,
				 pi_cutoff=0.05, nodes=None, species=None,
				 prefix='paralog_index', **hog_kargs):
		self.ogfile = ogfile
		self.orthfiles = orthfiles
		self.sptreefile = sptreefile
		self.self_synteny = self_synteny or orthfiles
		if not self.self_synteny:
			import sys
			logger.error('Paralog indexing needs synteny blocks: '
						 'provide -ss or -s')
			sys.exit(1)
		self.min_n = min_n
		self.gff = gff
		self.min_dist = min_dist
		self.pi_cutoff = pi_cutoff
		self.nodes = nodes
		self.species = species
		self.prefix = prefix
		self.hog_kargs = hog_kargs

		# lazy
		self._branch_pairs = None
		self._root_branch = None

	def _load_branch_pairs(self):
		"""Load HOGs and paralog pairs via Paralog."""
		if self._branch_pairs is not None:
			return
		paralog = Paralog(ogfile=self.ogfile, orthfiles=self.orthfiles,
						  sptreefile=self.sptreefile, prefix=self.prefix,
						  nodes=self.nodes, species=self.species,
						  **self.hog_kargs)
		self._branch_pairs, self._root_branch = paralog.run()
		# count paralog pairs per (branch, species)
		self._branch_sp_counts = defaultdict(int)
		for branch, pairs in self._branch_pairs.items():
			for g1, g2 in pairs:
				sp = g1.split('|')[0]
				self._branch_sp_counts[(branch, sp)] += 1

	def _compute_pi(self, block_pairs, branch_pairs):
		"""Return (BPI, n_intersect) for branch paralogy index."""
		if not block_pairs:
			return 0.0, 0
		common = sum(1 for p in block_pairs if p in branch_pairs)
		return common / len(block_pairs), common

	def _canonical_pair(self, g1, g2):
		return (g1, g2) if g1 < g2 else (g2, g1)

	def assign(self):
		"""Assign blocks to branches, return {branch: [(block, sp1, sp2, N, pi), ...]}."""
		self._load_branch_pairs()
		assigned = defaultdict(list)
		threshold = self.pi_cutoff
		root = self._root_branch
		# branch column order: internal nodes first, then leaves, each in
		# tree traversal order; leaves follow tree topology clustering.
		# Without a species tree, fall back to sorted branch order.
		if self.sptreefile:
			from .tree import number_nodes
			tree = number_nodes(self.sptreefile)
			internal = [n.name for n in tree.traverse() if not n.is_leaf()]
			leaves = tree.get_leaf_names()[::-1]  # reverse: outgroups first
			tree_order = internal + leaves
			branches = [b for b in tree_order if b in self._branch_pairs]
			branches += sorted(set(self._branch_pairs) - set(branches))
		else:
			branches = sorted(self._branch_pairs)
		self._pi_branches = branches
		self._pi_rows = []  # (block_id, N, pi_vector) for heatmap

		for rc in XCollinearity(self.self_synteny, gff=self.gff):
			if rc.N < self.min_n:
				continue
			if rc.species1 != rc.species2:
				continue
			if self.min_dist and rc.is_tandem(self.min_dist):
				continue
			if self.species and rc.species1 not in self.species:
				continue

			block_pairs = [self._canonical_pair(g1, g2) for g1, g2 in rc.pairs]

			best_branch = root
			best_pi = 0.0
			best_nparalog = 0
			pi_vector = []
			for branch in branches:
				pi, n_paralog = self._compute_pi(block_pairs,
												 self._branch_pairs[branch])
				pi_vector.append(pi)
				if pi > best_pi:
					best_pi = pi
					best_nparalog = n_paralog
					best_branch = branch

			if best_pi < threshold:
				best_branch = root

			# store data immediately — rc is a shared mutable object
			assigned[best_branch].append(
				(rc.block, rc.species1, rc.N, best_nparalog, best_pi))
			self._pi_rows.append((rc.id, rc.N, best_branch, tuple(pi_vector)))
		self._pi_branches = branches
		logger.info('Assigned blocks to {} branches'.format(len(assigned)))
		return assigned

	def write_heatmap(self, assigned):
		"""Write block x branch PI matrix heatmap.

		Rows (blocks) are sorted by their PI vector (descending) so blocks
		with similar branch-PI profiles end up adjacent.  With cluster=True
		rows are ordered by hierarchical clustering instead.
		"""
		fpath = self.prefix + '.heatmap.tsv'
		branches = self._pi_branches
		branch_idx = {b: i for i, b in enumerate(branches)}
		if getattr(self, 'heatmap_cluster', False):
			# hierarchical clustering of rows
			import numpy as np
			from scipy.cluster.hierarchy import linkage, leaves_list
			M = np.array([r[3] for r in self._pi_rows], dtype=float)
			if M.size:
				Z = linkage(M, method='ward')
				order = leaves_list(Z)
				rows = [self._pi_rows[i] for i in order]
			else:
				rows = []
		else:
			# sort by assigned-branch column index, then PI vector (desc);
			# root-assigned blocks (no paralog signal) go last
			rows = sorted(self._pi_rows,
						  key=lambda r: (branch_idx.get(r[2], len(branches)),
										 tuple(-v for v in r[3])))

		with open(fpath, 'w') as fout:
			fout.write('#block	N	' + '	'.join(branches) + '\n')
			for bid, N, bch, vec in rows:
				fout.write('{}	{}	{}\n'.format(
					bid, N, '	'.join('{:.4f}'.format(v) for v in vec)))
		logger.info('Heatmap matrix written to {} ({} blocks x {} branches)'.format(
			fpath, len(rows), len(branches)))

		# figure
		import numpy as np
		import matplotlib
		matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		M = np.array([r[3] for r in rows], dtype=float)
		if M.size == 0:
			logger.warning('No blocks for heatmap')
			return

		# row heights: uniform, or scaled by gene count when --scale given
		Ns = np.array([r[1] for r in rows], dtype=float)
		method = getattr(self, 'heatmap_scale', None)
		if method:
			if method == 'linear':
				h = Ns
			elif method == 'log2':
				h = np.log2(1 + Ns)
			elif method == 'log10':
				h = np.log10(1 + Ns)
			elif method == 'sqrt':
				h = np.sqrt(Ns)
			else:  # 'log' = natural log
				h = np.log1p(Ns)
			h = h / h.max() * (1 - 1e-3) + 1e-3  # tiny floor, no artificial base
		else:
			h = np.ones(len(rows))
		# cumulative y positions (bottom of each row)
		y = np.concatenate([[0.0], np.cumsum(h)])

		fig, ax = plt.subplots(figsize=(max(6, 0.4 * len(branches)),
										min(7, max(3, 0.02 * len(rows)))))
		cmap = plt.get_cmap('YlOrRd')
		# draw each row as a broken barh segment per branch column
		for i in range(len(rows)):
			col_frac = M[i]
			segs = []
			for j, v in enumerate(col_frac):
				segs.append((j - 0.5, 1.0))  # (xstart, width)
			colors = [cmap(v) for v in col_frac]
			ax.broken_barh(segs, (y[i], h[i]),
						   facecolors=colors, edgecolors='none',
						   linewidths=0)
		ax.set_xlim(-0.5, len(branches) - 0.5)
		ax.set_ylim(0, y[-1])
		ax.invert_yaxis()  # first row on top, consistent with imshow default
		ax.set_xticks(range(len(branches)))
		ax.set_xticklabels(branches, rotation=90, fontsize=6)
		ax.xaxis.tick_top()  # ticks on top
		ax.tick_params(axis='x', which='both', top=True, bottom=False,
					   labeltop=True, labelbottom=False)
		ax.set_yticks([])
		ax.set_xlabel('Branch')
		ax.set_ylabel('Block')
		import matplotlib as mpl
		fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, 1),
										   cmap=cmap),
					 ax=ax, label='BPI', shrink=0.6)
		fig.tight_layout()
		fig.savefig(self.prefix + '.heatmap.pdf')
		fig.savefig(self.prefix + '.heatmap.png', dpi=150)
		plt.close(fig)
		logger.info('Heatmap written to {}.heatmap.pdf/png'.format(self.prefix))

	def write_stats(self, assigned):
		"""Write per-branch per-species statistics TSV, return stats dict."""
		fpath = self.prefix + '.stats.tsv'
		stats = defaultdict(lambda: [0, 0, 0, 0.0])
		for branch, items in assigned.items():
			for block, sp, N, n_paralog, pi in items:
				if sp is None:
					continue
				s = stats[(branch, sp)]
				s[0] += 1
				s[1] += N
				s[2] += n_paralog
				s[3] += pi

		with open(fpath, 'w') as fout:
			fout.write('#branch	species	paralog_pairs	syntenic_blocks	'
					   'syntenic_gene_pairs	syntenic_paralog_pairs	'
					   'mean_PI	weighted_PI\n')
			for (branch, sp), (blocks, gp, pp, sum_pi) in sorted(stats.items()):
				mean_pi = sum_pi / blocks if blocks else 0.0
				n_paralogs = self._branch_sp_counts.get((branch, sp), 0)
				wpi = pp / gp if gp else 0.0
				fout.write('{}	{}	{}	{}	{}	{}	{:.4f}	{:.4f}\n'.format(
					branch, sp, n_paralogs, blocks, gp, pp, mean_pi, wpi))
		logger.info('Stats written to {}'.format(fpath))
		return stats

	def plot_tree(self, assigned, stats):
		"""Render species tree with pies: [syntenic, non-syntenic] paralog pairs.

		Pie size maps to total paralog_pairs per branch; pie fractions show
		the proportion of paralog pairs found in syntenic blocks.
		"""
		import os
		if 'QT_QPA_PLATFORM' not in os.environ:
			os.environ['QT_QPA_PLATFORM'] = 'offscreen'
		from ete3 import TreeStyle, PieChartFace, TextFace, NodeStyle, faces as etefaces
		from .tree import number_nodes

		# aggregate per branch over species
		branch_syn = defaultdict(int)   # branch -> syntenic paralog pairs
		branch_all = defaultdict(int)   # branch -> total paralog pairs
		for (branch, sp), (blocks, gp, pp, sum_pi) in stats.items():
			branch_syn[branch] += pp
			branch_all[branch] += self._branch_sp_counts.get((branch, sp), 0)

		tree = number_nodes(self.sptreefile)
		# number of descendant species per node — normalize pie size
		node_nspecies = {}
		for node in tree.traverse():
			node_nspecies[node.name] = len(node.get_leaf_names())

		ns = NodeStyle()
		ns['size'] = 0
		ns['hz_line_width'] = 1
		ns['vt_line_width'] = 1

		# scale pie size by paralog pairs per descendant species
		branch_per_sp = {b: 1.0 * n / max(node_nspecies.get(b, 1), 1)
						 for b, n in branch_all.items()}
		max_n = max(branch_per_sp.values()) if branch_per_sp else 1
		scale = 60.0 / max_n  # max pie 60px

		def layout(node):
			node.set_style(ns)
			nid = node.name
			if nid in branch_all:
				all_n = branch_all[nid]
				if all_n == 0:
					return
				syn_n = branch_syn.get(nid, 0)
				non_n = all_n - syn_n
				pcts = [100.0 * syn_n / all_n, 100.0 * non_n / all_n]
				size = max(8, int(scale * branch_per_sp[nid]))  # pie diameter
				pie = PieChartFace(pcts, size, size,
								   colors=['#377eb8', '#d9d9d9'],
								   line_color=None)
				pie.opacity = 0.9
				if node.is_leaf():
					etefaces.add_face_to_node(pie, node, column=0)
				else:
					etefaces.add_face_to_node(pie, node, column=0, position='branch-right')
			if node.is_leaf():
				name_face = TextFace(node.name, fsize=10)
				etefaces.add_face_to_node(name_face, node, column=1)
			else:
				name_face = TextFace(nid, fsize=7, fgcolor='#888888')
				etefaces.add_face_to_node(name_face, node, column=0, position='branch-top')

		ts = TreeStyle()
		ts.layout_fn = layout
		ts.show_leaf_name = False
		ts.scale = 200
		ts.branch_vertical_margin = 10

		out_pdf = self.prefix + '.tree.pdf'
		out_png = self.prefix + '.tree.png'
		tree.render(out_png, tree_style=ts)
		tree.render(out_pdf, tree_style=ts)
		logger.info('Tree plot written to {}'.format(out_pdf))

	def write_blocks(self, assigned):
		"""Write assigned blocks per branch."""
		for branch, items in assigned.items():
			fpath = '{}.{}.blocks'.format(self.prefix, branch)
			with open(fpath, 'w') as fout:
				for block, sp, N, n_paralog, pi in items:
					fout.write(block)
		logger.info('Block files written for {} branches'.format(len(assigned)))




# ---------------------------------------------------------------------------
#  CLI entry
# ---------------------------------------------------------------------------

def xmain(**kargs):
	no_index = kargs.pop('no_index', False)
	hog_tsv = kargs.get('hog_tsv')
	inparalog = kargs.get('inparalog')
	if not no_index and not hog_tsv and not inparalog:
		if not kargs.get('ogfile') or not kargs.get('orthfiles'):
			import sys
			logger.error('Either -og/-s or --hog/--paralog must be provided')
			sys.exit(1)
	if no_index:
		kargs.pop('output', None)
		Paralog(**kargs).run()
	else:
		tree_plot = kargs.pop('tree_plot', False)
		heatmap = kargs.pop('heatmap', False)
		heatmap_cluster = kargs.pop('heatmap_cluster', False)
		heatmap_scale = kargs.pop('heatmap_scale', None)
		indexer = ParalogIndexer(**kargs)
		indexer.heatmap_cluster = heatmap_cluster
		indexer.heatmap_scale = heatmap_scale
		assigned = indexer.assign()
		stats = indexer.write_stats(assigned)
		indexer.write_blocks(assigned)
		if tree_plot:
			if kargs.get('sptreefile'):
				indexer.plot_tree(assigned, stats)
			else:
				logger.warning('--tree-plot requires -t; skipped')
		if heatmap:
			indexer.write_heatmap(assigned)
