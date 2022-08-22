import os
import pandas as pd
import sys
import numpy as np
import argparse

def genOGid(frames,cells,indir,regdf,splitgene):
    ## 1) generate OGid file:
    gene=frames[0]
    for i in range(1,len(cells)):
        gene=gene.merge(frames[i],how='outer',on='Gene')
    gene.fillna('None',inplace=True)
    n=gene.shape[0]
    Gene_OGID='OG'+np.arange(1,n+1).astype(str).astype('object')+'_1'
    NAME=gene.iloc[:,1]
    for i in range(2,len(cells)+1):
        NAME=NAME+','+gene.iloc[:,i]
    ogid={'Gene_OGID':Gene_OGID, 'NAME':NAME}
    ogid = pd.DataFrame(ogid, columns=['Gene_OGID','NAME'])
    outname=indir+'testdata_ogids.txt'
    print('OGID file:',outname)
    ogid.to_csv(outname, index = None, header=True,  sep='\t', mode='w',float_format='%g')
    
    ## 2) TFs_OGs.txt (get OGID for regulators!!)
    regdata=regdf[0]
    for i in range(1,len(cells)):
        regdata=regdata.merge(regdf[i],how='outer',on='Gene')
    regulators=regdata['Gene']
    ## find OG id of regulators (OGID=index+1):
    tfs=gene[gene['Gene'].isin(regulators)==True].index+1
    ## double check if regulator IDs are correct:
    #print("Are regulator IDs correct?",all(gene.loc[tfs-1]['Gene']==regulators))
    tfs=pd.DataFrame({'TF':tfs})
    tffilename=indir+'TFs_OGs.txt'
    print('regulator OG id file:',tffilename)
    tfs.to_csv(tffilename, index = None, header=None,  sep='\t', mode='w',float_format='%g')
    
    ## 3) AllGenes.txt (OGID for all genes)
    filename=indir+'AllGenes.txt'
    print('gene OG id file:',filename)
    geneid=list(range(1,n+1))
    geneid=pd.DataFrame({'gene':geneid})
    geneid.to_csv(filename, index = None, header=None,  sep='\t', mode='w',float_format='%g')

    ## split AllGenes.txt into 50 genes per run
    if splitgene:
        j=0
        outdir=indir+'ogids/'
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        for i in range(1,n+1,splitgene):
            ogs=list(range(i,min(i+splitgene,n+1)))
            #print(ogs)
            ogs=pd.DataFrame({'gene':ogs})
            filename=outdir+'AllGenes'+str(j)+'.txt'
            ogs.to_csv(filename, index = None, header=None,  sep='\t', mode='w',float_format='%g')
            j+=1
        print("Split genes into ogs 0 ->",j)


def main(args):
    indir=args.indir+'/'
    outdir=args.outdir+'/'
    regfile=args.regfile
    splitgene=args.splitgene
    cellsuffix=args.addcellsuffix
    files=pd.read_table(args.filelist,header=None)
    filelist=files[1].values
    
    ## 0) celltype_order.txt
    files[0].to_csv(indir+'celltype_order.txt',sep='\n',mode='w',index=None,header=False,float_format='%g')
    cells=files[0].values
    frames=[]
    if len(regfile)>0:
        reg=pd.read_table(regfile,header=None)

    regdf=[]
    for i in range(len(cells)):
        cell=cells[i]
        infile=filelist[i]
        print(infile)
        data = pd.read_table(infile, sep='\t',index_col=0)
        allgenes = data.index.str.split('_'+cell).str.join('')
        if cellsuffix==1:
            d=pd.DataFrame({'Gene':allgenes,cell:allgenes+'_'+cell})
        else:
            d=pd.DataFrame({'Gene':allgenes,cell:allgenes})
        frames.append(d)
        
        ## 1) ${cell}_allGenes.txt
        filename=indir+cell+'_allGenes.txt'
        d[cell].to_csv(filename,sep='\n',mode='w',index=None,header=False,float_format='%g')
        print(cell,"Number of genes:",d.shape[0])  
        
        ## 2) ${cell}_allregulators.txt
        filename=indir+cell+'_allregulators.txt'
        if len(regfile)>0:
            dreg=reg.loc[reg[0].isin(allgenes)]
        else:
            dreg=allgenes.tolist()
            if 'Expression' in dreg:
            	dreg.remove('Expression')
            dreg=pd.DataFrame({0:dreg})
        
        ## double check if regulators in genes:
        print(cell,"Are filtered regulators all in data?",all(dreg[0].isin(allgenes)))
        print(cell,"Number of filtered regulators:",dreg.shape[0])    
        if cellsuffix==1:
            df=pd.DataFrame({'Gene':dreg[0].values,cell:dreg[0].values+'_'+cell})
        else:
            df=pd.DataFrame({'Gene':dreg[0].values,cell:dreg[0].values})
        regdf.append(df)
        df[cell].to_csv(filename,sep='\n',mode='w',index=None,header=False,float_format='%g')
        
        ## 3) add cell name to the Gene column of data and overwite the file
        if not all(data.index==allgenes+'_'+cell) and cellsuffix==1:
            data['Gene']=allgenes+'_'+cell
            cols=data.columns.tolist()
            cols=[cols[-1]]+cols[:-1]
            data1=data.reindex(columns=cols)
            outname=indir+cell+'.table'
            print("Overwrite data:"+outname)
            print(data1.iloc[:4,:4])
            data1.to_csv(outname,sep='\t',mode='w',index=None,header=True,float_format='%g')

    genOGid(frames,cells,indir,regdf,splitgene)
    ## generate config file: testdata_config.txt
    if args.motifs:
        config=pd.DataFrame({0:cells,1:indir+cells+'.table',2:outdir+cells,3:indir+cells+'_allregulators.txt',4:indir+cells+'_allGenes.txt',5:indir+cells+'_network.txt'})
        #config.to_csv(indir+'testdata_config_motifs.txt',sep='\t',mode='w',index=None,header=False,float_format='%g')
        config.to_csv(indir+'testdata_config.txt',sep='\t',mode='w',index=None,header=False,float_format='%g')
    else:
        config=pd.DataFrame({0:cells,1:indir+cells+'.table',2:outdir+cells,3:indir+cells+'_allregulators.txt',4:indir+cells+'_allGenes.txt',5:indir+'motif_empty.txt'})
        config.to_csv(indir+'testdata_config_noprior.txt',sep='\t',mode='w',index=None,header=False,float_format='%g')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(

        description=__doc__,

        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--filelist',

                        help='file list',

                        type=str,

                        default='')
    parser.add_argument('--regfile',

                        help='regulator file',

                        type=str,

                        default='')
    parser.add_argument('--indir',

                        help='directory of input data',

                        type=str,

                        default='')
    parser.add_argument('--outdir',

                        help='output directory',

                        type=str,

                        default='Results/')
    parser.add_argument('--motifs',

                        help='add motifs as prior or not',

                        type=bool,

                        default=0)
    parser.add_argument('--splitgene',

                        help='split all genes into X genes per run(0: no split, >0: X)',

                        type=int,

                        default=0)
    parser.add_argument('--addcellsuffix',

                        help='add cell name to the end of gene for each cell: gene_cellname',

                        type=int,

                        default=1)
    args = parser.parse_args()
    main(args)




