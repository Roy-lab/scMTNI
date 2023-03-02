
nseed=$1 # number of subsamples
maxReg=$2 # number of regulators 
cellfile=$3 # celltype_order.txt
indir=$4 # subsample/
nogid=$5

##################################################
# Merge separate networks from different gene sets
##################################################
echo "check if the results are ready"

samples=()
notdone=0
num_ogid=$((nogid-1))


for r in `seq 0 1 $((nseed-1))`
#for r in {0..1}
do
    #cd $indir/randseed${r}/
    
    a=$(ls ${indir}/randseed${r}/ogid*/Results.tar.gz | wc -l)
	if [[ $a != $num_ogid ]]
    then 
        echo randseed${r}
        notdone=1
    fi
    #cd -
done

echo "samples with error:" ${samples[@]}

#cd -

for r in `seq 0 1 $((nseed-1))`
#for r in {0..1}
do
    for og in `seq 0 1 $((nogid-1))`
	#for og in {0..1}
    do
        path0=$indir/randseed${r}/ogid${og}/
        cd $path0
        tar xzvf Results.tar.gz
		cd -
	done
done


for cell in `cat $cellfile`
do
   bash Scripts/Postprocessing/Makeconsensusnetwork1.sh $nseed $maxReg $cell $indir $nogid
done



