#!/bin/sh
Truefile=$1
Predfile=$2
outdir=$3
prefix=$4
regulators=$5
targets=$6
python Scripts/Evaluation/prepareAUClistfile_intersect.py --truefile $Truefile --predfile $Predfile --outdir $outdir --prefix $prefix --regulators $regulators --targets $targets

listfile=${outdir}/${prefix}_list.txt
outfile=${outdir}/${prefix}
AUC_jar=Scripts/Evaluation/auc.jar
echo java -jar ${AUC_jar} ${listfile} list
java -jar ${AUC_jar} ${listfile} list > $outfile
tail -n 2 $outfile > ${outdir}/${prefix}.out.txt

