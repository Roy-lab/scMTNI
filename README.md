# single-cell Multi-Task learning Network Inference (scMTNI)


We have developed single-cell Multi-Task learning Network Inference (scMTNI), a multi-task learning framework for joint inference of cell type-specific gene regulatory networks that leverages the cell lineage structure and scRNA-seq and scATAC-seq mea- surements to enable robust inference of cell type-specific gene regulatory networks. scMTNI takes as input a cell lineage tree, cell type-specific scRNA-seq data and optional cell type-specific prior networks that can be derived from bulk or single-cell ATAC-seq datasets. 

The scMTNI model has the following benefits: 
- 1 uses multi-task learning allowing the learning procedure to be informed by the shared infor- mation across cell types, 
- 2 incorporates the lineage structure to influence the extent of sharing between the learned networks, 
- 3 incorporates prior information, such as motif-based prior network derived from scATAC-seq data, thereby integrating scRNA-seq and scATAC-seq data to infer gene regulatory network dynamics across cell lineages.

Preprint: https://biorxiv.org/cgi/content/short/2022.07.25.501350v1

![alt text](figures_doc/scMTNI_v2.png)

## Step 1. Install
The code is compiled and tested for Linux environment. 
GSL (GNU Scientific Library) is used to handle matrix-related and vector-related operations.
It requires GCC version of gcc-6.3.1 and GNU extension with std=gnu++14 setting. The typical install time on a "normal" desktop computer is a few minutes.

1) git clone https://github.com/Roy-lab/scMTNI.git 
2) cd scMTNI/Code/ 
3) make


## Step 2. Prepare input files
The data for demo is in ExampleData/. The demo data contains 100 regulators and 300 genes.

## 2.1 integrating scRNA-seq and scATAC-seq using LIGER
Apply LIGER to integrate the scRNA-seq and scATAC-seq datasets, check LIGER (https://github.com/welch-lab/liger) for details. 
Input example files for scATAC-seq and scRNA-seq: ExampleData/LIGER/scATACseq.txt, ExampleData/LIGER/scRNAseq.txt

```
Rscript --vanilla Scripts/Integration/LIGER_scRNAseq_scATAC.R
```
The output files are in ExampleData/LIGER/. The liger cluster assginment is in ExampleData/LIGER/ligerclusters.txt

## 2.2 generating the prior network using scATAC-seq data and motifs
Check https://github.com/Roy-lab/scMTNI/blob/master/Scripts/genPriorNetwork/readme.md for details. 
Due to limitation of file size in Github, bam files are currently not provided in ExampleData/. For demo, please directly use the output prior networks ExampleData/cluster*_network.txt 
```
bash Scripts/genPriorNetwork/genPriorNetwork_scMTNI.sh
```
The example output files are ExampleData/cluster*_network.txt 

## 2.3 Prepare all input files and config file for scMTNI
### First prepare filelist.txt
The first column is the cell name, the second column is the location and filename of the expression data for each cell type. The example file ExampleData/filelist.txt:
```cluster8	ExampleData/cluster8.table
cluster3	ExampleData/cluster3.table
cluster2	ExampleData/cluster2.table
cluster1	ExampleData/cluster1.table
cluster6	ExampleData/cluster6.table
cluster9	ExampleData/cluster9.table
cluster10	ExampleData/cluster10.table
cluster7	ExampleData/cluster7.table
```

### Then prepare all the other input files based on ExampleData/filelist.txt and regulators list ExampleData/regulators.txt
Prepare input files with prior network:
```
indir=ExampleData/
filelist=${indir}/filelist.txt
regfile=${indir}/regulators.txt
python Scripts/PreparescMTNIinputfiles.py --filelist $filelist --regfile $regfile --indir $indir --outdir Results --splitgene 50 --motifs 1
```

Prepare input files without prior network:
```
python Scripts/PreparescMTNIinputfiles.py --filelist $filelist --regfile $regfile --indir $indir --outdir Results --splitgene 50 --motifs 0
```

### Prepare cell lineage tree:
The cell lineage tree file should have 5 columns describing the tree:
- 1. Child cell
- 2. Parent cell
- 3. Branch-specific gain rate (The probability that an edge is gained in a child given that the edge is absent in the predecessor cell)
- 4. Branch-specific loss rate (The probability that an edge is lost in a child given that the edge is present in the predecessor cell)

The example file for cell lineage tree ExampleData/celltype_tree_ancestor.txt

```cluster3	cluster8	0.2	0.2
cluster2	cluster3	0.2	0.2
cluster1	cluster2	0.2	0.2
cluster6	cluster2	0.2	0.2
cluster9	cluster6	0.2	0.2
cluster10	cluster6	0.2	0.2
cluster7	cluster10	0.2	0.2
```

## Step 3. Run
The input data for demo is in ExampleData/. The expected output is in Results/. The estimuated run time for the demo is around 7 minute.
The output network for each cell type is Results/cluster*/fold0/var_mb_pw_k50.txt

### Example uasge of scMTNI with prior network
```
Code/scMTNI -f ExampleData/testdata_config.txt -x50 -v1 -l ExampleData/TFs_OGs.txt -n ExampleData/AllGenes.txt -d ExampleData/celltype_tree_ancestor.txt -m ExampleData/testdata_ogids.txt -s ExampleData/celltype_order.txt -p 0.2 -c yes -b -0.9 -q 2 
```
The above example will run scMTNI using all regulators and targets. 


Since scMTNI learns regulators on a per-target basis, the algorithm can easily be parallelized by running the algorithm for each target gene (or sets of genes) separately. For example, to run scMTNI using 10 genes, we can replace the -n parameter with a file that contains only 10 genes as in ExampleData/AllGenes0.txt:

```
Code/scMTNI -f ExampleData/testdata_config.txt -x50 -v1 -l ExampleData/TFs_OGs.txt -n ExampleData/AllGenes0.txt -d ExampleData/celltype_tree_ancestor.txt -m ExampleData/testdata_ogids.txt -s ExampleData/celltype_order.txt -p 0.2 -c yes -b -0.9 -q 2 
```

### Example uasge of scMTNI without prior network
```
Code/scMTNI -f ExampleData/testdata_config_noprior.txt -x50 -v1 -l ExampleData/TFs_OGs.txt -n ExampleData/AllGenes.txt -d ExampleData/celltype_tree_ancestor.txt -m ExampleData/testdata_ogids.txt -s ExampleData/celltype_order.txt -p 0.2 -c yes -b -0.9 -q 0
```

### Parameter Explanations
f : config file with six columns, rows for each cell. Each cell's row should have the following species-specific entries:
- 1. Cell Name
- 2. Location of expression data with file name (cell.table)
- 3. Location to place outputs
- 4. List of regulators to be used
- 5. List of target genes to be used
- 6. List of motifs to be used. This file should have three tab-separated columns, listing the regulator, target, and motif score

x : Maximum # of regulators to be used for a given target.

v : Used for cross-validation. Can be left at 1.

p : default 0.5. The probability that an edge is present in the root cell.

l : List of the orthogroups (id #s) to be considered as regulators. Note: a regulator must also be present in the species-specific list of regulators given in the species-specific config file (parameter f)

n : List of the orthogroups (id #s) to be considered as targets. Note: a target must also be present in the species-specific list of targets given in the species-specific config file (parameter f)

d : The cell lineage tree to be used. This file should have 5 columns describing the tree:
- 1. Child cell
- 2. Parent cell
- 3. Branch-specific gain rate (The probability that an edge is gained in a child given that the edge is absent in the predecessor cell)
- 4. Branch-specific loss rate (The probability that an edge is lost in a child given that the edge is present in the predecessor cell)

m : A file describing the gene relationships. The first column of this file is of the format OGID{NUMBER}_{DUP}. Each NUMBER represents an orthogroup. For orthogroups with duplications, DUP is the duplication count/id. If there are no duplications in the dataset being used, DUP will always be 1.

s : A list of the cells present in the gene file (parameter m), in the order they exist in the gene file



## Step 4. Evaluation
### 4.0 Generate consensus network for subsample results

```
# Input:
# 1 [number of subsamples]
# 2 [number of regulators]
# 3 [cell file order file path]
# 4 [scMTNI results directory path]
# 5 [number of ogids]

nseed=10
maxReg=50
cellfile=ExampleData/celltype_order.txt
indir=Results_subsample/
nogid=50

bash Scripts/Postprocessing/Makeconsensusnetwork.sh $nseed $maxReg $cellfile $indir $nogid 

```

### 4.1 Compute AUPR:
```
bash Scripts/Evaluation/aupr_wrapper_list_intersection.sh
```
 
### 4.2 Compute F-score for top k edges compared to gold standard datasets

```
cellfile=ExampleData/celltype_order.txt
python Scripts/Evaluation/fscore_filterPred.py  --inferred $predicted_net --gold ${GSfile} --regulators $regulators --targets $targets --outdir $outpath
```

## Step 5. Network dynamics analysis
### 5.1 edge-based k-means clustering analysis:
Apply k-means clustering on edge*cell confidence matrix to find subnetworks with different patterns of conservation.

```
cd Scripts/Network_Analysis/
matlab -nodisplay -nosplash -nodesktop -r "StablityKmeansClustering; quit()"
```

### 5.2 topic model-based dynamic network analysis:
Apply LDA topic models to examine subnetwork level rewiring 

Prepare input matrix for subsample results
```
# Input:
# 1 [cell file order file path]
# 2 [scMTNI results directory path]
# 3 [edge filter threshold; i.e. 'full', '\_cf0.8'] 
Rscript Scripts/Network_Analysis/prepareNetmatrix.R ExampleData/celltype_order.txt Results_subsample/ \_cf0.8
```

Apply LDA
```
cd Scripts/Network_Analysis/
matlab -nodisplay -nosplash -nodesktop -r "LDA_analysis; quit()"
```

Obtain giant component of inferred networks
```
bash Scripts/Network_Analysis/getGiantComponents_H72noprior.sh
```

Obtain top regulators per component
```
bash Scripts/Network_Analysis/getTopRegPerComponent.sh
```

Make networks per topic: 
```
Rscript Scripts/Network_Analysis/makeAllGraphs.R
```




