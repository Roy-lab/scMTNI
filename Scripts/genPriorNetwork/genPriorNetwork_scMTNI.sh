########################################################################################################################################
## Step 1: extract bam files
########################################################################################################################################
prefix=liger_sqrt_ncell50_k8_filterhumanbc_FBS;dataset=FBS
#clusterfile=/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/liger/sqrt_atacsum100_rnacell50/k10/ligerclusters.txt
#datadir=/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/data/$prefix/
#indir=/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/data/scATACseq/${dataset}/bams_bycluster/
#cellfile=${datadir}/celltype_order.txt
datadir=../../ExampleData/
cellfile=../../ExampleData/celltype_order.txt
clusterfile=../../ExampleData/LIGER/ligerclusters.txt
indir=../../ExampleData/bams_bycluster/

for cellcluster in `cat $cellfile`
do
    python PrepareBamForLigerclusters.py $cellcluster $clusterfile $indir
done


########################################################################################################################################
## Step 2: MACS2 peaks calling
## /mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/scripts/macs2Callpeaks.sh
########################################################################################################################################
outdir=${datadir}/macs2_peaks/
mkdir -p $outdir
cd $indir
for cluster in `cat $cellfile`
do
    macs2 callpeak -t ${indir}/${cluster}.bam -f BAM -g mm -n $outdir/${cluster}_peak -B 
done

########################################################################################################################################
## Step 3: Map peaks to genes 
## Script: /mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/scripts/mapPeaksTogenes.sh
########################################################################################################################################
motifs=${datadir}/motifs/all_motifs_sorted_clean.txt
outdir2=${datadir}/mapPeaksTogenes/
mkdir -p $outdir2
bedtools=./bedtools 
promoterfile=${datadir}/motifs/Mus_musculus.GRCm38.74.TSS.5000.bed
for cluster in `cat $cellfile`
do 
    $bedtools intersect -a ${outdir}/${cluster}_peak_peaks.narrowPeak -b $motifs -wb -sorted > ${outdir2}/motifs_in_${cluster} 
    $bedtools intersect -a ${outdir}/${cluster}_peak_peaks.narrowPeak -b $promoterfile -wb -sorted > ${outdir2}/TSS_in_${cluster}
done

########################################################################################################################################
#  Step 4: connect TF motfis to genes
# /mnt/dv/wid/projects5/Roy-singlecell/sridharanlab/scATACseq/from_Stefan/filtered_fragments/ESC_8/others/mapMot2Gene.py
########################################################################################################################################
outdir3=${datadir}/macs2_networks/
mkdir -p $outdir3
motif2tf=${datadir}/motifs/cisbp_motif2tf.txt
for cluster in `cat $cellfile`
do 
    python mapMot2Gene.py --mot2tf $motif2tf --mot2peak ${outdir2}/motifs_in_${cluster} --peak2gene ${outdir2}/TSS_in_${cluster} --outfile ${outdir3}/${cluster}_network.txt
done

####################################################################################
## Step 5: filter regulators and genes in prior network:
#################################################################################### 
outdir4=${datadir}/macs2_motifs/
mkdir -p $outdir4
for sample in `cat $cellfile`
do
    regfile=${datadir}/${sample}_allregulators.txt
    genefile=${datadir}/${sample}_allGenes.txt
    netfile0=${outdir3}/${sample}_network.txt
    netfile=${outdir3}/${sample}.txt
    awk -v c=$sample 'BEGIN { OFS ="\t" }  $1=$1"_"c, $2=$2"_"c' $netfile0 > ${netfile}

    outfile=${outdir4}/${sample}_network.txt
    python filterpriornetwork.py --regfile $regfile --genefile $genefile --netfile $netfile --outfile $outfile 
done

####################################################################################
## Step 6: get top 0.2*tfs*genes edges:
#################################################################################### 

Rscript --vanilla filtertop20Pedges.R $datadir $outdir4


####################################################################################
## Step 7: perform percentile ranking
#################################################################################### 
outdir4=${datadir}/macs2_motifs_top0.2/
outdir5=${datadir}/macs2_prior_percentile/
mkdir -p $outdir5
## "incr" implies that the higher the score the better. "decr" implies that the lower the score the better.
for sample in `cat $cellfile`
do
    networkfile=${outdir4}/${sample}_network.txt
    outfile=${outdir5}/${sample}_network.txt
    echo $networkfile $outfile
    rankEdges $networkfile $outfile incr
done
