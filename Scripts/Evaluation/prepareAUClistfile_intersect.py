#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  3 22:59:02 2021

@author: shiluzhang
use intersection of tfs and targets in true network and predicted network, total possible = tfs*targets - self-loops;

"""


import os
import pandas as pd
import sys
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from itertools import product


def expand_grid(dictionary):
   return pd.DataFrame([row for row in product(*dictionary.values())], 
                       columns=dictionary.keys())


            
def main(args):
    tfile=args.truefile ##prefix='chromatin_marks_atac_expr_motif'
    pfile=args.predfile
    outdir=args.outdir
    prefix=args.prefix

    regulators=pd.read_table(args.regulators, sep='\n',header=None)
    regulators=regulators[0].values
    genes=pd.read_table(args.targets, sep='\n',header=None)
    genes=genes[0].values
    #tfile='/mnt/dv/wid/projects2/Roy-regnet-inference/CVN/data/chronis2017/SimulatedData/simulatedstructureReg30/simultopo_pgain0.1_pmain0.838_proot0.1/chromatin_marks_atac_expr_qmotif/mef_truenetwork_undirected.txt'
    Tnet = pd.read_table(tfile, sep='\t',header=None)
    if Tnet.iloc[0,0]=='TF':
        Tnet = pd.read_table(tfile, sep='\t',skiprows=1,header=None)
    Tnet['edge']=Tnet.iloc[:,0]+'->'+Tnet.iloc[:,1]
    Tnet['ytrue']=1
    
    #pfile='/mnt/dv/wid/projects2/Roy-regnet-inference/CVN/ChromNet/Results/GGMSimulatedData/simulatedstructureReg30/simultopo_pgain0.1_pmain0.838_proot0.1/chromatin_marks_atac_expr_qmotif/subsample//analysis_merged_top40//mef_consensus_edges.txt'
    Pnet=pd.read_table(pfile,header=None)
    Pnet['edge']=Pnet[0]+'->'+Pnet[1]
    Pnet.rename(columns={2:'score'},inplace=True)

    #tfs=list(set(Tnet.iloc[:,0]).intersection(set(Pnet[0])))
    #targets=list(set(Tnet.iloc[:,1]).intersection(set(Pnet[1])))
    tfs=list(set(Tnet[0]).intersection(regulators))
    targets=list(set(Tnet[1]).intersection(genes))

    ## filter True network with tfs and targets:
    id1=Tnet.index[Tnet[0].isin(tfs)]
    id2=Tnet.index[Tnet[1].isin(targets)]
    tid=id1.intersection(id2)
    if len(tid)==0:
        print("gold standard does not have tfs/genes present in the dataset")
        quit()
    Tnet_filter=Tnet.loc[tid,]
    idrm0=Tnet_filter.index[Tnet_filter[0]==Tnet_filter[1]]
    Tnet_filter.drop(idrm0,inplace=True)

    Tnet_filter.iloc[:,0:2].to_csv(outdir+'gold_filtrered.txt',index=None,header=None,sep='\t',mode='w',float_format='%g')
    tfs_filter=list(set(Tnet_filter[0]))
    if len(tfs_filter)>1:
        targets_filter=list(set(Tnet_filter[1]))
    else:
        targets_filter=genes
        print("Only 1 tf:"+tfs_filter[0])
    ## filter Predicted network with tfs and targets:
    idx1=Pnet.index[Pnet[0].isin(tfs_filter)]
    idx2=Pnet.index[Pnet[1].isin(targets_filter)]
    pid=idx1.intersection(idx2)
    if len(pid)==0:
        print("predicted network does not have tfs/genes present in the gold standard")
        quit()
    Pnet_filter=Pnet.loc[pid,]
    Pnet_filter.iloc[:,0:3].to_csv(outdir+'predicted_filtrered.txt',index=None,header=None,sep='\t',mode='w',float_format='%g')
     
    ## all possible for tfs*targets in filtered true network:
    dictionary={'tf':tfs_filter,'target':targets_filter}
    allpossible=expand_grid(dictionary)
    idrm=allpossible.index[allpossible['tf']==allpossible['target']]
    allpossible.drop(idrm,inplace=True)
    allpossible['edge']=allpossible['tf']+'->'+allpossible['target']

    # java AUC script:
    # list: prob(example == true)  true classification (1 == positive, 0 = negative)
    dp0=pd.merge(allpossible.loc[:,['edge']],Tnet_filter.loc[:,['edge','ytrue']],on='edge',how='left')
    dp0['ytrue']=dp0['ytrue'].fillna(0)
            
    dp=pd.merge(dp0, Pnet.loc[:,['edge','score']],on='edge',how='left')
    dp['score']=dp['score'].fillna(0)
    print("Number of True Edges:",int(sum(dp['ytrue']))," Number of total possible:",dp.shape[0])

    dp.sort_values(['score'],ascending=0,axis=0,inplace=True)
    dp.reset_index(drop=True,inplace=True)
    outfile=outdir+prefix+'_list.txt'
    print(outfile)
    dp.loc[:,['score','ytrue']].to_csv(outfile,index=None,header=None,sep='\t',mode='w',float_format='%g')
        
    precision, recall, thresholds = precision_recall_curve(dp['ytrue'], dp['score'])
    average_precision = average_precision_score(dp['ytrue'], dp['score'])
    print('Average precision-recall score: {0:0.6f}'.format(average_precision))    
    # plt.figure(figsize=(5,5))
    # plt.plot(recall, precision, marker='.')
    # plt.xlim(0, 1)  
    # plt.ylim(0, 1)  
    # plt.title("PR curve "+cell+" AUPR="+"{0:0.2f}".format(average_precision), fontsize = 10)
    # # axis labels
    # plt.xlabel('Recall')
    # plt.ylabel('Precision')
    # # show the legend
    # plt.legend()
    # # show the plot
    # plt.show()
    # plt.savefig('PRcurve_'+cell+'_960.pdf')

    
    
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(

        description=__doc__,

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--truefile',

                        help='predfile',

                        type=str,

                        default='')
    
    parser.add_argument('--predfile',

                        help='predfile',

                        type=str,

                        default='')
    parser.add_argument('--outdir',

                        help='outdir',

                        type=str,

                        default='')
    parser.add_argument('--prefix',

                        help='prefix',

                        type=str,

                        default='')
    parser.add_argument('--regulators',

                        help='regulators',

                        type=str,

                        default='')   
    parser.add_argument('--targets',

                        help='targets',

                        type=str,

                        default='')  
    args = parser.parse_args()
    main(args)


