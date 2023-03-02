source("Scripts/Network_Analysis/genNetall_multi.R");

giantcomponent="Results_subsample/analysis/lda_TFcellbygene/network_cf0.8/k10/"
k=10
filter="_cf0.8"
buenrostro_cellorder=c("cluster8", "cluster3","cluster2","cluster1","cluster6","cluster9","cluster10", "cluster7")

for (tp in 1:k)
{
	#netfile=paste0(giantcomponent,"topic",tp,"_directed_giantcomponent.txt")
	netfile=paste0(giantcomponent,"topic",tp,"_directed_giantcomponent.txt")
  	prefix=paste0("lda_topic",tp,filter)
  	#outname=paste0(giantcomponent,prefix,"_k",k,"_network.png")
  	#outname=paste0(giantcomponent,prefix,"_k",k,"_network_eigen.pdf")
  	outname=paste0(giantcomponent,prefix,"_k",k,"_network_8by11_min10.pdf")
	print(paste0("Now plotting ",outname))
	dtopic1=read.table(netfile,header=T,stringsAsFactors=F);
  	#genNetall_multi(prefix,dtopic1,outname,w=5000,h=3000,nr=2,ntop=20,cellorder)
  	genNetall_multi(prefix,dtopic1,outname,w=5000,h=3000,nr=3,ntop=10,buenrostro_cellorder)
	
}
