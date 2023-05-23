# generate prior network using scATAC-seq data

####################################################################################
## Step 1: extract bam files
```datadir=../../ExampleData/
cellfile=../../ExampleData/celltype_order.txt
clusterfile=../../ExampleData/LIGER/ligerclusters.txt
indir=../../ExampleData/bams_bycluster/
for cellcluster in `cat $cellfile`
do
    python PrepareBamForLigerclusters.py $cellcluster $clusterfile $indir
done
```

####################################################################################
## Step 2: MACS2 peaks calling

```outdir=${datadir}/macs2_peaks/
mkdir -p $outdir
cd $indir
for cluster in `cat $cellfile`
do
    macs2 callpeak -t ${indir}/${cluster}.bam -f BAM -g mm -n $outdir/${cluster}_peak -B 
done
```

####################################################################################
## Step 3: Map peaks to genes 

```motifs=${datadir}/motifs/all_motifs_sorted_clean.txt
outdir2=${datadir}/mapPeaksTogenes/
mkdir -p $outdir2
bedtools=./bedtools
promoterfile=${datadir}/motifs/Mus_musculus.GRCm38.74.TSS.5000.bed
for cluster in `cat $cellfile`
do 
    $bedtools intersect -a ${outdir}/${cluster}_peak_peaks.narrowPeak -b $motifs -wb -sorted > ${outdir2}/motifs_in_${cluster} 
    $bedtools intersect -a ${outdir}/${cluster}_peak_peaks.narrowPeak -b $promoterfile -wb -sorted > ${outdir2}/TSS_in_${cluster}
done
```

####################################################################################
##  Step 4: connect TF motfis to genes
## Note: mapMot2Gene.py only works with narrowPeak format
## See example peak file format: /mnt/dv/wid/projects5/Roy-singlecell/shilu_work/scMTNI/code/scMTNI/ExampleData/macs2_peaks/cluster8_peak_peaks.narrowPeak
```
chr1	3425131	3425514	cluster8_peak_peak_1	110	.	8.01443	13.75227	11.07777	184
chr1	3444972	3445147	cluster8_peak_peak_2	45	.	5.08926	6.88859	4.58632	93
```

```outdir3=${datadir}/macs2_networks/
mkdir -p $outdir3
motif2tf=${datadir}/motifs/cisbp_motif2tf.txt
for cluster in `cat $cellfile`
do 
    python mapMot2Gene.py --mot2tf $motif2tf --mot2peak ${outdir2}/motifs_in_${cluster} --peak2gene ${outdir2}/TSS_in_${cluster} --outfile ${outdir3}/${cluster}_network.txt
done
```

####################################################################################
## Step 5: filter regulators and genes in prior network:

```outdir4=${datadir}/macs2_motifs/
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
```

####################################################################################
## Step 6: get top 0.2 * tfs * genes edges:
```
Rscript --vanilla filtertop20Pedges.R $datadir $outdir4
```

####################################################################################
## Step 7: perform percentile ranking
```outdir4=${datadir}/macs2_motifs_top0.2/
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
```
