# 04OCT2022 SGM
# wrapper for getting consensus across all cell types
# Output used as input for clustering

extractTopconsensusedges=/mnt/dv/wid/projects7/Roy-singlecell2/ashtonlab/results/H72_de/scmtni/merlin/scripts/extractTopconsensusedges.R

cf=cf
n=0.8
outputpath=/mnt/dv/wid/projects7/Roy-singlecell2/ashtonlab/results/H72_de/scmtni/merlin/subsample/analysis
#mkdir ${outputpath}
Rscript --vanilla $extractTopconsensusedges $cf $n ${outputpath}

