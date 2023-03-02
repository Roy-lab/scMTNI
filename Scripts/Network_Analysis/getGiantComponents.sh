TOPICNETS=Results_subsample/analysis/lda_TFcellbygene/network_cf0.8/k10
OUTDIR=Results_subsample/analysis/lda_TFcellbygene/network_cf0.8/k10/

mkdir -p ${OUTDIR}
TOOL=Scripts/Network_Analysis/componentsanalysis_directed/getDirectedComponents

T=1

while [ $T -le 10 ]
do
	OUTSUFF=${OUTDIR}/topic${T}_directed
	if [ -f ${TOPICNETS}/networks_lda_k10_topic${T}.txt ]
	then
		echo "Now getting giant component for ${T}"
		${TOOL} ${TOPICNETS}/networks_lda_k10_topic${T}.txt ${OUTSUFF} 10
	fi
	T=`expr ${T} + 1`
done
