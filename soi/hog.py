import sys
import re
import numpy as np
import networkx as nx
from collections import defaultdict, Counter
from itertools import combinations

from .OrthoFinder import OrthoMCLGroup
from .mcscan import ColinearGroups
from .tree import number_nodes
from .RunCmdsMP import logger

def xmain(**kargs):
	HOG(**kargs).pipe()
	
class HOG:
	def __init__(self, ogfile=None, orthfiles=None, sptreefile=None, outpre = "HOGs",
		  paralog=False, max_copies=5, out_stats=False, bar_plot=False, tree_plot=False, 
			min_child_species=1, cross_speciation=False, drop_no_cross=False, **kargs):
		self.ogfile = ogfile
		self.orthfiles = orthfiles
		self.sptreefile = sptreefile
		self.outtsv = outpre + '.tsv'
		self.noparalog = not paralog
		self.max_copies = max_copies
		self.min_child_species = min_child_species
		self.cross_speciation = cross_speciation
		self.drop_no_cross = drop_no_cross
		self.out_stats = outpre + '.stats.tsv' if out_stats else None
		self.bar_plot = outpre + '.bar' if bar_plot else None
		self.tree_plot = outpre + '.tree' if tree_plot else None
	def from_tsv(self, fpath):
		"""Load HOGs from a previously written HOGs.tsv (skip rebuilding).

		Fills self.tree/species/all_hogs/children in the same shape as pipe().
		Without a species tree, self.tree stays None and species are inferred
		from gene prefixes.
		"""
		if self.sptreefile:
			from .tree import number_nodes
			self.tree = number_nodes(self.sptreefile)
			self.species = self.tree.get_leaf_names()
		else:
			self.tree = None
			self.species = []
		self.all_hogs = {}
		self.all_genes = []
		with open(fpath) as f:
			header = f.readline()
			for line in f:
				parts = line.rstrip().split('	')
				if len(parts) < 5:
					continue
				hog_id, og_id, node_id, parent, gene_str = parts[:5]
				genes = gene_str.split()
				self.all_genes.extend(genes)
				species = sorted(set(g.split('|')[0] for g in genes))
				self.all_hogs[hog_id] = HOGrecord(
					hog_id=hog_id, og_id=og_id, node_id=node_id,
					genes=genes, species=species,
					parent=None if parent in ('None', 'Root') else parent,
					children=[])
		# rebuild children links from parent column
		for hog in self.all_hogs.values():
			if hog.parent and hog.parent in self.all_hogs:
				self.all_hogs[hog.parent].children.append(hog.hog_id)
		logger.info('Loaded {} HOGs from {}'.format(len(self.all_hogs), fpath))
		return self

	def pipe(self, write_tsv=True):
		logger.info(f'Reading and Numbering species tree from {self.sptreefile}')
		if not self.sptreefile:
			raise ValueError('rebuild needs -t species tree (HOG rebuilding requires a species tree)')
		self.tree = sptree = number_nodes(self.sptreefile)
		self.species = sptree.get_leaf_names()

		# precompute child species sets for cross-speciation check
		node_child_species = {}
		if self.cross_speciation or self.drop_no_cross:
			for node in sptree.traverse():
				if not node.is_leaf():
					node_child_species[node.name] = [set(c.get_leaf_names()) for c in node.get_children()]
		
		logger.info(f'Loading orthologs from {self.orthfiles}')
		graph = ColinearGroups(self.orthfiles, spsd=self.species, 
					noparalog=self.noparalog).graph
		
		logger.info(f'Splitting HOGs from {self.ogfile}')
		prefix = 'hog'
		self.all_hogs = {}
		self.all_genes = []
		n_drop_no_cross = 0
		n_no_split = 0
		for og in OrthoMCLGroup(self.ogfile):
			og_genes = og.genes
			self.all_genes  += og_genes
			og_species = set(og.species)
			og_id = og.ogid
			og_spdict = og.spdict
			subgraph = graph.subgraph(og_genes).copy()

			node_to_hogs = defaultdict(list)

			for node in sptree.traverse(strategy="postorder"):
				if self.noparalog and node.is_leaf():
					continue
				node_id = node.name
				node_species = set(node.get_leaf_names())
				subset_species = node_species & og_species
				if not subset_species:
					continue

				subset_genes = []
				for sp in subset_species:
					subset_genes.extend(og_spdict[sp])

				if not subset_genes:
					continue

				# cross-speciation: genes must span all child branches of this node
				do_split = True
				if not node.is_leaf() and (self.cross_speciation or self.drop_no_cross):
					child_sps = node_child_species[node_id]
					crosses = all(subset_species & cs for cs in child_sps)
					if not crosses:
						if self.drop_no_cross:
							n_drop_no_cross += 1
							continue
						do_split = False
						n_no_split += 1

				if do_split:
					hog_subgraph = subgraph.subgraph(subset_genes)
					connected_components = list(nx.connected_components(hog_subgraph))
				else:
					connected_components = [set(subset_genes)]
				connected_components.sort(key=lambda comp: min(comp))
				for idx, cc_genes in enumerate(connected_components):
					# species actually present in this connected component
					cc_species = [sp for sp in og_spdict
								  if sp in subset_species and set(og_spdict[sp]) & cc_genes]
					# filter orphan HOGs: internal node with genes from < min_child_species species
					if not node.is_leaf() and self.min_child_species > 1:
						if len(cc_species) < self.min_child_species:
							continue
					hog_id = f"{og_id}.{node_id}.{prefix}{idx}"

					hog = HOGrecord(
						hog_id=hog_id,
						og_id=og_id,
						node_id=node_id,
						genes=list(cc_genes),
						species=cc_species,
						parent=None,
						children=[]
					)
					self.all_hogs[hog_id] = hog
					node_to_hogs[node_id].append(hog)

			for node in sptree.traverse(strategy="postorder"):
				node_id = node.name
				current_hogs = node_to_hogs.get(node_id, [])

				if not current_hogs:
					continue

				parent_node = node.up
				if not parent_node:
					continue

				p_node_id = parent_node.name
				parent_hogs = node_to_hogs.get(p_node_id, [])

				for hog in current_hogs:
					hog_genes = set(hog["genes"])
					for phog in parent_hogs:
						phog_genes = set(phog["genes"])
						if hog_genes.issubset(phog_genes):
							hog["parent"] = phog["hog_id"]
							phog["children"].append(hog["hog_id"])
							break

			# drop leaf HOGs with no parent (genes absent from parent-node
			# HOGs, e.g. filtered out by min_child_species at the parent level)
			n_dropped = 0
			leaf_ids = {n.name for n in sptree.traverse() if n.is_leaf()}
			for node_id, hogs in list(node_to_hogs.items()):
				if node_id not in leaf_ids:
					continue
				for hog in list(hogs):
					if not hog["parent"]:
						del self.all_hogs[hog["hog_id"]]
						node_to_hogs[node_id].remove(hog)
						n_dropped += 1
			if n_dropped:
				logger.debug(f"Dropped {n_dropped} orphan leaf HOGs without parent "
							 f"(SOG {og_id})")
		logger.info(f"Processed {len(self.all_hogs)} HOGs")
		if n_drop_no_cross:
			logger.info(f"Dropped {n_drop_no_cross} gene sets that do not cross all child branches (--drop-no-cross)")
		if n_no_split:
			logger.info(f"Kept {n_no_split} gene sets as single HOGs — do not cross all child branches (--cross-speciation)")
		logger.info("All HOGs with hierarchy built successfully!")
		
		if write_tsv:
			with open(self.outtsv, 'w', encoding='utf-8') as f:
				self.write_all_hogs_in_one_file(f)
			logger.info(f"HOGs written to {self.outtsv}")

		# compute and output copy-number statistics (only when needed)
		need_stats = self.out_stats or self.bar_plot or self.tree_plot
		if need_stats:
			leaf_data, internal_data, node_names = self.compute_copy_stats()
			if self.out_stats and (leaf_data or internal_data):
				self.write_stats_table(leaf_data, internal_data, node_names)
			if self.bar_plot:
				self.plot_stats(leaf_data, internal_data, node_names)
			if self.tree_plot:
				self.plot_tree(leaf_data, internal_data, node_names)
		return self.all_hogs
		
	def write_all_hogs_in_one_file(self, fout=sys.stdout):
		header = ["HOG", "OG", "Tree_Node", "Parent", "Genes"]
		fout.write("\t".join(header) + "\n")

		for hog in self.all_hogs.values():
			hog_id = hog["hog_id"]
			og_id = hog["og_id"]
			node_id = hog["node_id"]
			if hog["parent"] is not None:
				parent_id = hog["parent"]
			elif node_id == self.tree.name:
				parent_id = "Root"
			else:
				parent_id = "None"
			gene_str = " ".join(sorted(hog["genes"]))

			row = [hog_id, og_id, node_id, parent_id, gene_str]
			fout.write("\t".join(row) + "\n")

	# ------------------------------------------------------------------
	#  shared branch-level grouping
	# ------------------------------------------------------------------

	def _group_hogs_by_branch(self):
		"""Group HOGs by branch for downstream aggregation.

		Internal nodes: for each node N, group child HOGs by their parent
		  HOG (the HOG at N's parent).  Each parent->children group
		  represents a duplication event on the branch leading into N.

		Leaf nodes: for each species, group its genes by HOG.
		  Multiple genes of the same species in the same HOG indicate
		  within-species paralogy.

		Returns:
		  internal_groups: {node_id: {parent_hog_id: [child_hog_records]}}
		  leaf_groups:     {species: {hog_id: [genes_of_that_species]}}
		"""
		all_nodes = list(self.tree.traverse(strategy="postorder")) if self.tree else []
		leaf_set = {n.name for n in all_nodes if n.is_leaf()}
		internal_set = {n.name for n in all_nodes if not n.is_leaf()}

		# Without a species tree: infer leaf species from gene prefixes,
		# internal nodes from node_ids not matching any species.
		if not self.tree:
			sp_from_genes = set()
			for hog in self.all_hogs.values():
				for sp in hog["species"]:
					sp_from_genes.add(sp)
			leaf_set = sp_from_genes
			internal_set = set()
			for hog in self.all_hogs.values():
				nid = hog["node_id"]
				if nid and nid not in leaf_set and nid not in ("Root", "None"):
					internal_set.add(nid)

		# --- Internal nodes ---
		internal_groups = defaultdict(lambda: defaultdict(list))
		for hog in self.all_hogs.values():
			nid = hog["node_id"]
			if nid not in internal_set:
				continue
			pid = hog["parent"]
			if not pid or pid == "Root":
				if not getattr(self, '_enable_root_paralog', False):
					continue
				# Root-level sibling HOGs: group by og_id as the "parent" key
				pid = hog["og_id"]
			parent_hog = self.all_hogs.get(pid)
			if not parent_hog:
				continue
			internal_groups[nid][pid].append(hog)

		# --- Leaf species ---
		# Only use HOGs at the immediate parent node of each species.
		# Genes appearing only at higher nodes (e.g. ginseng in N6)
		# are not counted — those duplications belong to deeper branches.
		leaf_parent = {}
		if self.tree:
			for leaf in all_nodes:
				if leaf.is_leaf() and leaf.up:
					leaf_parent[leaf.name] = leaf.up.name
		else:
			# Without a tree: the direct parent node of each species is the
			# deepest internal node whose HOGs contain that species' genes.
			sp_to_nodes = {}
			for hog in self.all_hogs.values():
				nid = hog["node_id"]
				if nid in internal_set:
					for sp in hog["species"]:
						sp_to_nodes.setdefault(sp, set()).add(nid)
			# deepest = the node farthest from root; without tree topology
			# we cannot rank depth — use the node appearing in the most
			# species' gene sets (closest to leaves) as a heuristic
			node_sp_count = {}
			for nids in sp_to_nodes.values():
				for nid in nids:
					node_sp_count[nid] = node_sp_count.get(nid, 0) + 1
			for sp, nids in list(sp_to_nodes.items()):
				if nids:
					leaf_parent[sp] = max(nids,
										  key=lambda n: node_sp_count.get(n, 0))

		leaf_groups = defaultdict(lambda: defaultdict(list))
		for hog in self.all_hogs.values():
			for sp in hog["species"]:
				if sp not in leaf_parent:
					continue
				if hog["node_id"] != leaf_parent[sp]:
					continue
				sp_genes = [g for g in hog["genes"] if g.startswith(sp)]
				if sp_genes:
					leaf_groups[sp][hog["hog_id"]].extend(sp_genes)

		return internal_groups, leaf_groups

	# ------------------------------------------------------------------
	#  copy-number statistics (uses _group_hogs_by_branch)
	# ------------------------------------------------------------------

	def compute_copy_stats(self):
		"""Compute per-node copy-number distribution.

		Internal nodes: child count per parent HOG → distribution.
		Leaf species:   gene copies per HOG → distribution.

		Returns (leaf_data, internal_data, node_names).
		"""
		internal_groups, leaf_groups = self._group_hogs_by_branch()
		all_nodes = list(self.tree.traverse(strategy="postorder"))
		leaf_set = {n.name for n in all_nodes if n.is_leaf()}
		max_c = self.max_copies

		# Internal
		internal_data = []
		internal_names = []
		for node in all_nodes:
			nid = node.name
			if nid in leaf_set:
				continue
			pcounts = internal_groups.get(nid, {})
			if not pcounts:
				continue
			dist = Counter()
			for children in pcounts.values():
				c = min(len(children), max_c)
				dist[c] += 1
			arr = np.array(sorted(dist.items()), dtype=int)
			if len(arr):
				internal_data.append(arr)
				internal_names.append(nid)

		# Leaf
		leaf_data = []
		leaf_names = []
		for node in all_nodes:
			if not node.is_leaf():
				continue
			sp = node.name
			sp_hog_genes = leaf_groups.get(sp, {})
			if not sp_hog_genes:
				continue
			dist = Counter()
			for genes in sp_hog_genes.values():
				c = min(len(genes), max_c)
				dist[c] += 1
			arr = np.array(sorted(dist.items()), dtype=int)
			if len(arr):
				leaf_data.append(arr)
				leaf_names.append(sp)

		return leaf_data, internal_data, leaf_names + internal_names

	# ------------------------------------------------------------------
	#  branch paralog iterator
	# ------------------------------------------------------------------

	def iter_branch_paralogs(self, nodes=None, species=None):
		"""Yield paralogous gene pairs produced at each branch.

		Internal nodes: sibling child HOGs (same parent) → gene pairs
		  across different child HOGs for the same species.

		Leaf species: genes of the same species within the same HOG
		  → pairwise paralogs.

		Yields: (gene1, gene2, node_id, species, hog_id)
		"""
		internal_groups, leaf_groups = self._group_hogs_by_branch()
		node_set = set(nodes) if nodes else None
		sp_set = set(species) if species else None

		# Helper: group genes by species for a list of HOG records
		_hog_sp_genes = {}
		def _get_sp_genes(hog_rec):
			if hog_rec.hog_id not in _hog_sp_genes:
				sg = defaultdict(list)
				for g in hog_rec["genes"]:
					sp = g.split('|')[0] if '|' in g else g.rsplit('_', 1)[0]
					sg[sp].append(g)
				_hog_sp_genes[hog_rec.hog_id] = dict(sg)
			return _hog_sp_genes[hog_rec.hog_id]

		# --- Internal nodes: cross-child-HOG paralogs ---
		for nid, parent_hogs in internal_groups.items():
			if node_set and nid not in node_set:
				continue
			for parent_hog_id, children in parent_hogs.items():
				if len(children) < 2:
					continue
				for ci in range(len(children)):
					sp_genes_i = _get_sp_genes(children[ci])
					for cj in range(ci + 1, len(children)):
						sp_genes_j = _get_sp_genes(children[cj])
						common = set(sp_genes_i) & set(sp_genes_j)
						if sp_set:
							common &= sp_set
						for sp in common:
							for g1 in sp_genes_i[sp]:
								for g2 in sp_genes_j[sp]:
									yield (g1, g2, nid, sp, parent_hog_id)

		# --- Leaf species: within-HOG paralogs ---
		for sp, hog_genes in leaf_groups.items():
			if sp_set and sp not in sp_set:
				continue
			for hog_id, genes in hog_genes.items():
				if len(genes) < 2:
					continue
				for g1, g2 in combinations(sorted(set(genes)), 2):
					yield (g1, g2, sp, sp, hog_id)

	def write_stats_table(self, leaf_data, internal_data, node_names):
		"""Write TSV like save_depth_table: rows = nodes, cols = copy numbers."""
		max_c = self.max_copies
		all_data = leaf_data + internal_data
		header = ["Node"] + [str(i) for i in range(1, max_c)] + [f"{max_c}+", "Multi%"]
		rows = [header]
		for name, arr in zip(node_names, all_data):
			counts = {i: 0 for i in range(1, max_c + 1)}
			for copies, cnt in arr:
				c = min(int(copies), max_c)
				counts[c] += int(cnt)
			total = sum(counts.values())
			multi = sum(v for k, v in counts.items() if k > 1)
			ratio = f"{100.0 * multi / total:.1f}" if total else "0"
			row = [name] + [str(counts[i]) for i in range(1, max_c + 1)] + [ratio]
			rows.append(row)
		text = "\n".join(["\t".join(r) for r in rows]) + "\n"
		with open(self.out_stats, 'w', encoding='utf-8') as f:
			f.write(text)
		logger.info(f"Copy-number stats written to {self.out_stats}")

	def plot_stats(self, leaf_data, internal_data, node_names):
		"""Bar chart: one panel per node showing copy-number distribution."""
		from .ploidy_plotter import plot_bars
		all_data = leaf_data + internal_data
		if not all_data:
			logger.warning("No data for plot")
			return
		n = len(all_data)
		nrow = int(np.ceil(np.sqrt(n)))
		ncol = int(np.ceil(n / nrow))
		outfigs = [f"{self.bar_plot}.pdf", f"{self.bar_plot}.png"]
		plot_bars(all_data, node_names, outfigs=outfigs,
				nrow=nrow, ncol=ncol, max_ploidy=self.max_copies,
				xlabel='Copy number', ylabel='HOG count')
		logger.info(f"Bar plot written to {self.bar_plot}.pdf, {self.bar_plot}.png")

	def plot_tree(self, leaf_data, internal_data, node_names):
		"""Render species tree with pie charts at nodes showing copy-number distribution."""
		import os
		if 'QT_QPA_PLATFORM' not in os.environ:
			os.environ['QT_QPA_PLATFORM'] = 'offscreen'
		from ete3 import TreeStyle, PieChartFace, TextFace, NodeStyle, faces as etefaces
		pie_colors = ['#d9d9d9', '#377eb8', '#4daf4a', '#ff7f00', '#e41a1c',
					  '#984ea3', '#a65628', '#f781bf', '#999999', '#66c2a5']
		max_c = self.max_copies
		all_data = leaf_data + internal_data

		node_pcts = {}
		for name, arr in zip(node_names, all_data):
			counts = {i: 0 for i in range(1, max_c + 1)}
			for copies, cnt in arr:
				c = min(int(copies), max_c)
				counts[c] += int(cnt)
			total = sum(counts.values())
			if total == 0:
				continue
			pcts = [100.0 * counts[i] / total for i in range(1, max_c + 1)]
			node_pcts[name] = pcts

		# node style: hide dots, keep lines
		ns = NodeStyle()
		ns['size'] = 0
		ns['hz_line_width'] = 1
		ns['vt_line_width'] = 1

		def layout(node):
			node.set_style(ns)
			nid = node.name
			if nid in node_pcts:
				pcts = node_pcts[nid]
				non_zero = [(p, c) for p, c in zip(pcts, pie_colors[:max_c]) if p > 0.5]
				if non_zero:
					filtered_pcts = [p for p, c in non_zero]
					filtered_colors = [c for p, c in non_zero]
					pie = PieChartFace(filtered_pcts, 20, 20,
									   colors=filtered_colors, line_color=None)
					pie.opacity = 0.9
					pie.margin_left = 4
					pie.margin_right = 4
					if node.is_leaf():
						etefaces.add_face_to_node(pie, node, column=0) #, position='aligned')
					else:
						etefaces.add_face_to_node(pie, node, column=0, position='branch-right')
			if node.is_leaf():
				name_face = TextFace(node.name, fsize=10)
				etefaces.add_face_to_node(name_face, node, column=1) #, position='aligned')
			else:
				name_face = TextFace(nid, fsize=7, fgcolor='#888888')
				etefaces.add_face_to_node(name_face, node, column=0, position='branch-top')

		tree = self.tree
		ts = TreeStyle()
		ts.layout_fn = layout
		ts.show_leaf_name = False
		ts.scale = 200
		ts.branch_vertical_margin = 10

		out_pdf = self.tree_plot + '.pdf'
		out_png = self.tree_plot + '.png'
		tree.render(out_png, tree_style=ts)
		tree.render(out_pdf, tree_style=ts)
		logger.info(f"Tree plot written to {out_pdf}, {out_png}")


class HOGrecord:
    __slots__ = ["hog_id", "og_id", "node_id", "genes", "species", "parent", "children"]
    
    def __init__(self, hog_id, og_id, node_id, genes, species, parent=None, children=None):
        self.hog_id = hog_id
        self.og_id = og_id
        self.node_id = node_id
        self.genes = genes
        self.species = species
        self.parent = parent
        self.children = children if children is not None else []

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __hash__(self):
        return hash(self.hog_id)

    def __eq__(self, other):
        return isinstance(other, HOGrecord) and self.hog_id == other.hog_id

    def __repr__(self):
        return f"{self.hog_id}"



def main():
	HOG(ogfile=sys.argv[1], orthfiles=sys.argv[2], sptreefile=sys.argv[3]).pipe()
if __name__ == '__main__':
	main()
