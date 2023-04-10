DIR=Results_subsample/analysis/lda_TFcellbygene/k10

T=1
while [ $T -le 10 ]
do
	netfile=${DIR}/topic${T}_directed_giantcomponent.txt
	regfile=${DIR}/topic${T}_giantcomponent_min10_regdegree.txt
	regorderfile=${DIR}/topic${T}_min10_yorder.txt


	awk '{printf("%s\t%s\n",$1,$4)}' ${netfile}|sort |uniq -c |egrep -v regulator |awk '{printf("%s\t%s\t%d\n",$2,$3,$1)}'  > ${regfile}
	awk '{printf("%s\t%s\n",$1,$4)}' ${netfile} |sort |uniq -c |egrep -v regulator |awk '{printf("%s\t%s\t%d\n",$2,$3,$1)}' |cut -f1 |sort |uniq -c |sort -k1,1rg |awk '{printf("%s\n",$2)}' > ${regorderfile}
	echo "awk '{printf("%s\t%s\n",$1,$4)}' ${netfile} |sort |uniq -c |egrep -v regulator |awk '{printf("%s\t%s\t%d\n",$2,$3,$1)}' |cut -f1 |sort |uniq -c |sort -k1,1rg |awk '{printf("%s\n",$2)}' > ${regorderfile}"
	T=`expr ${T} + 1`

done

