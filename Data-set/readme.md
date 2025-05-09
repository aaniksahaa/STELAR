### Datasets for "STELAR: A statistically consistent coalescent-based species trees estimation method by maximizing triplet consistency".
___________________________________________________________________________________________________________________________________________
### Note:
-  (astral, stelar, mpest, superTriplests species trees are represented with `.astral.5.6.1` , `.stelar`, `.mpest` and `.SuperTriplets` extensions respectively)`
- All the dataset analyzed here, have already being created in other studies as well. We have just `rooted, trimmed` the dataset for the purpouse of our study Source of the original 
dataset is [here](https://sites.google.com/eng.ucsd.edu/datasets/home?authuser=0). 
___________________________________________________________________________________________________________________________________________


**11-taxon datasets**: 10 replicates	 (four model conditions: 100gene-100bp, 100gene-1000bp, 1000gene-100bp, 1000gene-1000bp)

File contents:
  1. True species tree: 11-taxa/true-species-tree.strip
  2. estimated genes strongILS: 11-taxa/estimated_genes_strongILS/<X>/Rep<Y>_XgenesLS - where Y is the replicate number (1...10), and Y is the number of genes.
  3. simulated genes strongILS: 11-taxa/simulated_genes_strongILS/<X>/Rep<Y>_XgenesLS - where Y is the replicate number (1...10), and Y is the number of genes.
____________________________________________________________________________________________________________________________________________
**15-taxon datasets**: 10 replicates (four model conditions: 100gene-100bp, 100gene-1000bp, 1000gene-100bp, 1000gene-1000bp)

File contents:
  1. True species trees: 15-taxon/true-species.tre
  2. True gene trees: 15-taxon/true-genetrees/R<X>/<Y>/<Y>.tre - where X is the replicate number (1...10), and Y is the gene id (1...1000).
  3. MPEST species trees: 15-taxon/model-condition/estimated-genetrees/R<X>.mpest - where X is the replicate number (1...10), and Y is the gene id (1...1000).
  4. Astral species trees: 15-taxon/model-condition/estimated-genetrees/R<X>.astral - where X is the replicate number (1...10), and Y is the gene id (1...1000).
  5. Stelar species trees: 15-taxon/model-condition/estimated-genetrees/R<X>.stelar - where X is the replicate number (1...10), and Y is the gene id (1...1000).

____________________________________________________________________________________________________________________________________________

**Mammalian (37 taxa)**: 20 replicates. coalescence level = [noscale, scale2d, scale2u], number of genes = [25, 50, 100, 200, 400, 500], base pair = [500, 1000, 1500]
File contents:
  1. True species tree: 37-taxa/true.tre.strip
  2. Gene trees: 37-taxa/X.y.Z/R<K>/Best.1 - where X is the coalescence level, Y is the number of genes, Z is the base pair, K is the replicate number (1...20).
  3. Astral inferred species tree: 37-taxa/X.y.Z/R<K>/astral.1 - where X is the coalescence level, Y is the number of genes, Z is the base pair, K is the replicate number (1...20).
  4. Stelar inferred species tree: 37-taxa/X.y.Z/R<K>/stelar.1 - where X is the coalescence level, Y is the number of genes, Z is the base pair, K is the replicate number (1...20).
4. MPEST inferred species tree: 37-taxa/X.y.Z/R<K>/Best.1.tre - where X is the coalescence level, Y is the number of genes, Z is the base pair, K is the replicate number (1...20).
____________________________________________________________________________________________________________________________________________

Pleasei feel free to email bayzid AT cse.buet.ac.bd in case you have any questions.
