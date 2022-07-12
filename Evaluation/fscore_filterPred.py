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


def main(args):
    tfile=args.gold ##prefix='chromatin_marks_atac_expr_motif'
    pfile=args.inferred
    outdir=args.outdir
    #outfile=args.outfile
    #tfile='/mnt/dv/wid/projects2/Roy-common/data/data_new/mouse/known_interactions/mesc/mESC_chipunion.txt'
    #pfile='/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/Results/scCVN_upto1transition/liger_sqrt_ncell50_k10_filterhumanbc/pg0.2_pm0.8_pr0.2_maxReg50_b2_bm4/subsample/analysis/cluster4/consensus_edges.txt'
    #regulators='/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/data/liger_sqrt_ncell50_k10_filterhumanbc/merged/merged_allregulators.txt'
    #targets='/mnt/dv/wid/projects5/Roy-singlecell/shilu_work/integrate_scrna_scatac/networkinference/data/liger_sqrt_ncell50_k10_filterhumanbc/merged/merged_allGenes.txt'
    regulators=pd.read_table(args.regulators, sep='\n',header=None)
    regulators=regulators[0].values
    genes=pd.read_table(args.targets, sep='\n',header=None)
    genes=genes[0].values
    Tnet = pd.read_table(tfile, sep='\t',header=None)
    if Tnet.iloc[0,0]=='TF':
        Tnet = pd.read_table(tfile, sep='\t',skiprows=1,header=None)
    Tnet['edge']=Tnet.iloc[:,0]+'->'+Tnet.iloc[:,1]
    Tnet['ytrue']=1

    ## Step1: filter gold standard network with tfs and targets in the dataset:
    tfs=list(set(Tnet[0]).intersection(regulators))
    targets=list(set(Tnet[1]).intersection(genes))
    id1=Tnet.index[Tnet[0].isin(tfs)]
    id2=Tnet.index[Tnet[1].isin(targets)]
    tid=id1.intersection(id2)
    if len(tid)==0:
        print("gold standard does not have tfs/genes present in the dataset")
        print('{0:.1f}\t{1:.1f}\t{2:.1f}'.format(0, 0, 0))
        quit()
    Tnet_filter=Tnet.loc[tid,]
    ## remove self-loops in the gold standard dataset:
    idrm0=Tnet_filter.index[Tnet_filter[0]==Tnet_filter[1]]
    Tnet_filter.drop(idrm0,inplace=True)
    Tnet_filter.iloc[:,0:2].to_csv(outdir+'gold_filtrered.txt',index=None,header=None,sep='\t',mode='w',float_format='%g')
    tfs_filter=list(set(Tnet_filter[0]))
    targets_filter=list(set(Tnet_filter[1]))
    # if len(tfs_filter)>1:
    #     targets_filter=list(set(Tnet_filter[1]))
    # else:
    #     targets_filter=genes
    #     print("Only 1 tf:"+tfs_filter[0])

    Pnet0=pd.read_table(pfile,header=None)
    Pnet0['edge']=Pnet0[0]+'->'+Pnet0[1]
    Pnet0.rename(columns={2:'score'},inplace=True)
    #if args.top_edges is not None:
    idx1=Pnet0.index[Pnet0[0].isin(tfs_filter)]
    idx2=Pnet0.index[Pnet0[1].isin(targets_filter)]
    pid=idx1.intersection(idx2)
    if len(pid)==0:
        print("predicted network does not have tfs/genes present in the gs dataset")
        print('{0:.1f}\t{1:.1f}\t{2:.1f}'.format(0, 0, 0))
        quit()
    Pnet0=Pnet0.loc[pid,:]
    Pnet0.sort_values(['score'],ascending=0,axis=0,inplace=True)
    print("Number of Edges in filtered predicted networks:",Pnet0.shape[0])
    Pnet0.iloc[:,0:3].to_csv(outdir+'prednet_filtrered.txt',index=None,header=None,sep='\t',mode='w',float_format='%g')

    ## Compure F-score for top edges compared to gold standard dataset:
    dout=[]
    for top_edges in [100,200,300,400,500,1000,2000]: #[5000,10000,50000,100000]:
        Pnet=Pnet0.iloc[:top_edges,:] #Pnet=Pnet0.iloc[:args.top_edges,:]
        ## filter Predicted network with tfs and targets in filtered gold standard dataset, but this will make different algorithms not comparable
        #idx1=Pnet.index[Pnet[0].isin(tfs_filter)]
        #idx2=Pnet.index[Pnet[1].isin(targets_filter)]
        #pid=idx1.intersection(idx2)
        #if len(pid)==0:
            #print("predicted network does not have tfs/genes present in the gold standard")
            #print('{0:.1f}\t{1:.1f}\t{2:.1f}'.format(0, 0, 0))
            #continue
        #Pnet_filter=Pnet.loc[pid,]
        #Pnet_filter.iloc[:,0:3].to_csv(outdir+'predicted_filtrered.txt',index=None,header=None,sep='\t',mode='w',float_format='%g')
        # fscore:
        dp=pd.merge(Tnet_filter.loc[:,['edge','ytrue']], Pnet.loc[:,['edge','score']],on='edge',how='inner')
        n_TP = dp.shape[0]
        precision = n_TP / Pnet.shape[0]
        recall = n_TP / Tnet_filter.shape[0]
        dpp=dp['edge'].str.split('->',expand=True)
        if n_TP>0:
            tfs_overlap=len(set(dpp[0]))
            targets_overlap=len(set(dpp[1]))
        else:
            tfs_overlap=0
            targets_overlap=0
        try:
            fscore = 2 * (precision * recall) / (precision + recall)
        except ZeroDivisionError:
            fscore = 0
        d = {'edge': top_edges,'TP':n_TP,'tfs_overlap':tfs_overlap,'targets_overlap':targets_overlap,'gs':Tnet_filter.shape[0],'gs_tfs':len(tfs_filter),'gs_targets':len(targets_filter),'predicted':Pnet.shape[0],'predicted_tfs':len(set(Pnet[0])),'predicted_targets':len(set(Pnet[1])),'precision': precision,'recall':recall,'fscore':fscore}
        df = pd.DataFrame(d, columns=['edge','TP','tfs_overlap','targets_overlap','gs','gs_tfs','gs_targets','predicted','predicted_tfs','predicted_targets','precision','recall','fscore'], index=[0])
        print("Number of Edges in gold standard:",Tnet_filter.shape[0]," Number of edges in predicted:",Pnet.shape[0],"intersection:",n_TP)
        #print("fscore:",fscore)
        print('{0:.4f}\t{1:.4f}\t{2:.4f}'.format(precision, recall, fscore))  
        dout.append(df)

    dout = pd.concat(dout)
    dout.to_csv(outdir+'fscore.txt',index=None,header=True,sep='\t',mode='w',float_format='%g')
    print(outdir+'fscore.txt')
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(

        description=__doc__,

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--gold',

                        help='gold',

                        type=str,

                        default='')

    parser.add_argument('--inferred',

                        help='inferred',

                        type=str,

                        default='')
    # parser.add_argument('--top_edges',
    #                     type=int, default=None,
    #                     help="Number of edges to filter to")
    # parser.add_argument('-k', '--top-edges',
    #                     type=int, default=None,
    #                     help="Number of edges to filter to")
    parser.add_argument('--outdir',

                        help='outdir',

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




