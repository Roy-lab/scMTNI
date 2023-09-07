# generate prior network using scATAC-seq data

#### Note: All the files in ExampleData/ are subsamples of the original files, as a demo for the file format. The raw data is too large to upload. Please contact the Roy Lab for raw data if needed.

####################################################################################
## Step 1: extract bam files
```
datadir=../../ExampleData/
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

```
outdir=${datadir}/macs2_peaks/
mkdir -p $outdir
cd $indir
for cluster in `cat $cellfile`
do
    macs2 callpeak -t ${indir}/${cluster}.bam -f BAM -g mm -n $outdir/${cluster}_peak -B 
done
```

####################################################################################
## Step 3: Map peaks to genes 

#### See promoterfile example: ExampleData/motifs/Homo_sapiens.GRCh37.74.TSS.5000.bed
#### See example motifs file: ExampleData/motifs/all_motifs_sorted_clean.txt (This file is a subsample of the full file, as a demo for the file format. The raw motif file for mouse or human is too large to upload. The original motif files are available to download at Zenodo https://zenodo.org/record/8323399#.ZPlrKuzMLOE.
#### To generate the motif files, we obtained sequence-specific motifs from the Cis-BP database (http://cisbp.ccbr.utoronto.ca/) and used the script pwmmatch.exact.r available from the PIQ toolkit to identify significant motif instances genome-wide using the human genome assembly hg19 and mouse genome assembly mm10 respectively.


```
motifs=${datadir}/motifs/all_motifs_sorted_clean.txt
outdir2=${datadir}/mapPeaksTogenes/
mkdir -p $outdir2
bedtools=./bedtools
### promoter file for mouse:
promoterfile=${datadir}/motifs/Mus_musculus.GRCm38.74.TSS.5000.bed
### promoter file for human:
promoterfile=${datadir}/motifs/Homo_sapiens.GRCh37.74.TSS.5000.bed

for cluster in `cat $cellfile`
do 
    $bedtools intersect -a ${outdir}/${cluster}_peak_peaks.narrowPeak -b $motifs -wb -sorted > ${outdir2}/motifs_in_${cluster} 
    $bedtools intersect -a ${outdir}/${cluster}_peak_peaks.narrowPeak -b $promoterfile -wb -sorted > ${outdir2}/TSS_in_${cluster}
done
```

#### See example peak file format: /mnt/dv/wid/projects5/Roy-singlecell/shilu_work/scMTNI/code/scMTNI/ExampleData/macs2_peaks/cluster8_peak_peaks.narrowPeak
```
chr1    3425131 3425514 cluster8_peak_peak_1    110     .       8.01443 13.75227        11.07777        184
chr1    3444972 3445147 cluster8_peak_peak_2    45      .       5.08926 6.88859 4.58632 93
```

#### Example file of ${outdir2}/motifs_in_${cluster}
#### See ExampleData/mapPeaksTogenes/motifs_in_cluster8
```
chr1	10008	10009	cluster8_peak_peak_1	59	.	6.11179	8.98592	5.99756	74	chr1	9993	10009	M5970_1.02	5.638908	-
chr1	10008	10009	cluster8_peak_peak_1	59	.	6.11179	8.98592	5.99756	74	chr1	10001	10009	M1358_1.02	10.033302	+
chr1	10008	10015	cluster8_peak_peak_1	59	.	6.11179	8.98592	5.99756	74	chr1	10007	10015	M1358_1.02	10.033302	+
```

#### Example file of ${outdir2}/TSS_in_${cluster}
#### See ExampleData/mapPeaksTogenes/TSS_in_cluster8
```
chr1	10008	10162	cluster8_peak_peak_1	59	.	6.11179	8.98592	5.99756	74	chr1	6868	16868	TSS_ENST00000456328_DDX11L1	-1	+	DDX11L1
chr1	10008	10162	cluster8_peak_peak_1	59	.	6.11179	8.98592	5.99756	74	chr1	6871	16871	TSS_ENST00000515242_DDX11L1	-1	+	DDX11L1
chr1	10008	10162	cluster8_peak_peak_1	59	.	6.11179	8.98592	5.99756	74	chr1	6873	16873	TSS_ENST00000518655_DDX11L1	-1	+	DDX11L1
```

####################################################################################
##  Step 4: connect TF motfis to genes
#### Note: mapMot2Gene.py currently only works with narrowPeak format, column4 is the peak unique identifier, and the peak file contains 10 columns. It also requires the last column of the promoterfile to be gene name. It requires the second last column of the motifs file to be motif score, and the third last column of the motifs file to be motif name. 
#### See example peak file format: ExampleData/macs2_peaks/cluster8_peak_peaks.narrowPeak
#### See promoterfile example: ExampleData/motifs/Homo_sapiens.GRCh37.74.TSS.5000.bed
#### See example motifs file: ExampleData/motifs/all_motifs_sorted_clean.txt

```
outdir3=${datadir}/macs2_networks/
mkdir -p $outdir3
motif2tf=${datadir}/motifs/cisbp_motif2tf.txt
for cluster in `cat $cellfile`
do 
    python mapMot2Gene.py --mot2tf $motif2tf --mot2peak ${outdir2}/motifs_in_${cluster} --peak2gene ${outdir2}/TSS_in_${cluster} --outfile ${outdir3}/${cluster}_network.txt
done
```

####################################################################################
## Step 5: filter regulators and genes in prior network:

```
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
```

####################################################################################
## Step 6: get top 0.2 * tfs * genes edges:
```
Rscript --vanilla filtertop20Pedges.R $datadir $outdir4
```

####################################################################################
## Step 7: perform percentile ranking
```
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
```
