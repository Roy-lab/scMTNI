import sys
import pysam
import os
def readBarcodes(inname):
    cluster = {}
    f = open(inname,'r')
    for l in f:
        parts1 = l.strip().split('\t')
        parts2 = l.strip().split(',')
        if len(parts1) == 2:
            cluster[parts1[0]] = 'cluster'+parts1[1]
        else:
            cluster[parts2[0]] = 'cluster'+parts2[1]
    f.close()
    return cluster

def readData(c,cluster,outdir):  #readData(inname,outname,cluster,c):
    data   = {}
    revcluster = {}
    for b in cluster:
        cl = cluster[b]
        bset = revcluster.get(cl,set([]))
        bset.add(b)
        revcluster[cl] = bset

    #for c in revcluster:
    bset = revcluster[c]
    outname=outdir+'/'+c+'.bam' #c='cluster1'
    inname='../../ExampleData/bam/scATACseq.bam'
    print(outname)
    ncell=0
    f = pysam.AlignmentFile(inname,'rb')
    of = pysam.AlignmentFile(outname,'wb',template=f)
    cells=set()
    barcodes=set()
    for l in f:
        try:
            tags = l.get_tag('RG').split('-')
            cell=''
            for i in range(len(tags)):
                if tags[i] in ['GMP',"CMP","mono","HSC"]:
                    cell=tags[i]
                    break
            barcode=l.query_name.split(':')[0]
            cells.add(cell)
            #barcodes.add(barcode)
            key='atac_'+cell+'.'+barcode
            #print(key)
            if key in bset:
                of.write(l)
                ncell+=1
        except:
                pass

    f.close()
    of.close()
    print(c+":"+str(ncell))
    print(cells)

if __name__ == '__main__':
    clusterfile=sys.argv[2]
    outdir=sys.argv[3]
    os.makedirs(outdir,exist_ok=True)
    cluster = readBarcodes(clusterfile)
    readData(sys.argv[1],cluster,outdir) ##[2],sys.argv[3],cluster)

