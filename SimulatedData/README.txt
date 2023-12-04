## Simulated data for benchmarking the performance of different multi-task and single-task learning algorithms.

# Linear cell lineage:
celltype_order.txt

- hsc corresponds to Cell type 1,
- cmp corresponds to Cell type 2,
- gmp corresponds to Cell type 3.

hsc (C1) -> cmp (C2) -> gmp (C3)

# Dataset 1: consisting of 65 genes and 2000 cells for three cell types. 
hsc.table
cmp.table
gmp.table

# Dataset 2: downsampling dataset 1 to 1000 cells
hsc_n1000.table
cmp_n1000.table  
gmp_n1000.table  

# Dataset 3: downsampling dataset 1 to 200 cells
hsc_n200.table
cmp_n200.table 
gmp_n200.table  

# Simulated networks for each cell type:
Simulated_GRNs/
- hsc_network.tsv
- cmp_network.tsv  
- gmp_network.tsv  

# Input files for scMTNI:
## Genes for each cell type:
- cmp_allGenes.txt  
- gmp_allGenes.txt  
- hsc_allGenes.txt
## Regulators for each cell type
- cmp_allregulators.txt  
- gmp_allregulators.txt  
- hsc_allregulators.txt
## The orthogroup related input files
- testdata_ogids.txt
## Config file:
- config_condor.txt

