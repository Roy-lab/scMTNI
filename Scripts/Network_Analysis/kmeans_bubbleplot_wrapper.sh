

bubble_script=Scripts/Network_Analysis/bubbleplot.R

kstart=$1
kend=$2
cluster_dir=$3
out_dir=$4
cutoff=$5
prefix=$6

#for kmeans in {15..30}
for kmeans in $(seq ${kstart} ${kend})
	do
		#cluster_dir=Results_subsample/analysis/kmeansclustering_cf0.8/
		#out_dir=Results_subsample/analysis/kmeansclustering_cf0.8/

		clust_prefix=${prefix}_stability_kmeans_${cutoff}_k${kmeans}
		cluster_file=${cluster_dir}/${prefix}_stability_kmeans_${cutoff}_k${kmeans}.txt

		if [ -f "${cluster_file}" ]; then
				echo Rscript $bubble_script $kmeans $cluster_dir $out_dir $clust_prefix
				Rscript $bubble_script $kmeans $cluster_dir $out_dir $clust_prefix
				else echo ${cluster_file} does not exist
		fi

done


