# 04OCT2022 SGM
# wrapper for getting consensus across all cell types
# Output used as input for clustering

extractTopconsensusedges=Scripts/Postprocessing/extractTopconsensusedges.R

cf=cf
n=0.8
outputpath=Results_subsample/analysis
#mkdir ${outputpath}
Rscript --vanilla $extractTopconsensusedges $cf $n ${outputpath}
