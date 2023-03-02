
nseed=$1 #50
maxReg=$2 #50
cell=$3 # celltype_order.txt
indir=$4 # /subsample/
nogid=$5

#########################################################
# Step0: Merge separate networks from different gene sets
#########################################################

notdone=0
for r in `seq 0 1 $((nseed-1))` 
#for r in {0..1}
do
    outdir=$indir/randseed${r}/$cell/fold0
    mkdir -p $outdir
    rm $outdir/var_mb_pw_k${maxReg}.txt
    echo $outdir/var_mb_pw_k${maxReg}.txt
    for og in `seq 0 1 $((nogid-1))`  #{0..239}
	#for og in 0
    do
        path0=$indir/randseed${r}/ogid${og}/
        path=${path0}/Results/$cell/fold0
        echo $path/var_mb_pw_k${maxReg}.txt 
        if [[ -s $path/var_mb_pw_k${maxReg}.txt ]];then
           cat $path/var_mb_pw_k${maxReg}.txt >> $outdir/var_mb_pw_k${maxReg}.txt
        else
           echo $path/var_mb_pw_k${maxReg}.txt "empty"
           notdone=1
        fi
    done
    echo "Output: " $outdir
    echo " "
done


###############################
# Step1: Make consensus network
###############################

# Make a list of network files:
estEdge=../../Scripts/Postprocessing/estimateedgeconf/estimateEdgeConf
mkdir -p $indir/analysis/
cd $indir/analysis/
echo $cell
mkdir -p $cell
rm $cell/network_files.txt

for rseed in `seq 0 1 $((nseed-1))`
#for rseed in {0..1}
	do
		ls ../randseed${rseed}/$cell/fold0/var_mb_pw_k${maxReg}.txt >> $cell/network_files.txt
    done

${estEdge} $cell/network_files.txt 0 $cell/net_ alledges
sed "s/_${cell}//g;s/_${cell}//g" ${cell}/net_alledge.txt > ${cell}/consensus_edges.txt



