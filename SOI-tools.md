
#### Macro-synteny phylogeny ####
Before run the pipeline:
1. install the lasest verion of [SOI](https://github.com/zhangrengang/orthoindex#installation) by `pip3 install .`, 
and have a check: `which soi-syn`;
2. complete the [example pipeline](https://github.com/zhangrengang/evolution_example) to get the orthologous synteny file `collinearity.ortho`;
3. prepare file `species.config` to set the expected subgenome numbers for the targeted species (TAB seperated):
```
Vitis_vinifera  1
Aralia_elata    2
Centella_asiatica       2
```
The first species will be set as the outgroup by default.
Then run the pipeline:
```
cd phylogenomics
soi-syn anchor_trees collinearity.ortho species.config ../pep.faa ../all_species_gene.gff output_dir
```
After the pipeline completed, you can find tree files in the `output_dir`:
```
OG*treefile		# gene tree file for each anchor gene
chr*treefile	# macro-synteny tree file by concatenating the anchor genes from the same chromosome set
CHR*treefile	# macro-synteny tree file, but allowing gene missing
```
For example, `CHR_Ae11-Ae15-Ca2-Ca7-Vv1_143_283.concat.treefile`, `Ae11-Ae15-Ca2-Ca7-Vv1` is the chromosome set,
`143` is the number of anchor genes, and `283` is the number of all syntenic genes allowing missing.

#### Allele identification ####
For diploid and polyploid assembly, we may aim to identify allele genes. As allele genes are not always syntenic (e.g., non-syntenic due to samll-scale inversion), we combine synteny and orthology to identify a full set of allele genes.

After `soi filter` has identified orthologous synteny in the [example pipeline](https://github.com/zhangrengang/evolution_example), we run
```
soi-syn retrieve_allele collinearity.ortho ../OrthoFinder/OrthoFinder/Results_*/ ../all_species_gene.gff sps="Panax_ginseng Panax_notoginseng" min_block=10 win_size=10 > allele.txt
```
`sps="Panax_ginseng Panax_notoginseng"` sets one or more target species or pseudo-species (space seperated); 
species labels need to be different for different subgenomes (label gene ID like `SP1|sgA.g1`, `SP2|sgB.g1`);

`min_block=10` sets minimum length of syntenic blocks; increasing the value will result in more reliable alleles;

`win_size=10` sets upstream and downstream 10 genes to retrieve orthologs using syntenic genes as anchors; increasing the value will retrieve more alleles.

The output file is like:
```
                # primary alleles               # secondary alleles             # sources of primary gene pairs
chrom   idx     Panax_ginseng   Panax_notoginseng       Panax_ginseng   Panax_notoginseng       Panax_ginseng-Panax_notoginseng
1.0     2552.0  GWHGBEIL002170.1        PN011184        -       -       orthology:15422
1.0     2552.0  GWHGBEIL002171.1        PN011183        -       -       orthology:26434
1.0     2552.5  GWHGBEIL002169.1        PN011185        -       -       synteny-1:51
1.0     2557.0  GWHGBEIL002172.1        PN000315        -       -       synteny-1:40865
1.0     2558.0  GWHGBEIL002173.1        PN000316        -       -       synteny-1:None
1.0     2560.0  GWHGBEIL002175.1        PN021828        -       -       synteny-1:None
1.0     2561.0  GWHGBEIL002176.1        PN021827        -       -       synteny-1:None
1.0     2562.0  GWHGBEIL002177.1        PN021826        -       -       synteny-1:35052
1.0     2564.0  GWHGBEIL002178.1        PN021824        -       -       synteny-1:37814
1.0     2569.5  GWHGBEIL002183.1        PN026247        -       -       synteny-1:2117
1.0     2570.5  GWHGBEIL002184.1        PN026246        -       -       synteny-1:42347
1.0     2571.5  GWHGBEIL002179.1        PN026953        -       -       orthology:40400
....
```
The 3rd and 4th columns indicate allelic gene pairs. 
The last column indicates the source of alleles; for example, `orthology:15422` means ortholog #15422, 
`synteny-1:51` means it is from syntenic block #1 and ortholog #51; `synteny-1:None` means that it is not pre-inferred as an ortholog.
The column number will extend for polyploids.

#### Orthology format conversion ####
The output formats of different orthology inference tools are different. It maybe better to convert them into a unified format (gene pairs). Here is an example:
```
soi-syn homologs OrthoFinder/OrthoFinder/Results_*/ > homologs.txt
```
The output file is in a simplest pair format (`gene1<TAB>gene2`), containing all orthologs and inparalogs inferred by OrthoFinder. 
Outputs from some other orthology inference tools are also supported. See details for [these formats](README.md#orthology-format).
