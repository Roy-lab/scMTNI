import sys
import argparse

parser = argparse.ArgumentParser('map motifs to genes through peaks')
parser.add_argument("--mot2tf",    type=str, help="Motif to TFs",        required=True)
parser.add_argument("--mot2peak",  type=str, help="Motif to Peaks",      required=True)
parser.add_argument("--peak2gene", type=str, help="Peak to Genes",       required=True)
parser.add_argument("--outfile",   type=str, help="Output (TF to Gene)", required=True)
args = parser.parse_args()

def readMot2TF(inname):
	mot2tf = {}
	f = open(inname,'r')
	for l in f:
		#M6519_1.02	Tgif2::Tgif1
		parts = l.strip().split('\t')
		mot2tf[parts[0]] = parts[1].split('::')
	f.close()
	return mot2tf

def readPeak2Gene(inname):
	peak2gene = {}
	f = open(inname,'r')
	for l in f:
		#chr1	3460841	3460939	ESC.8_peak_11	37	.	4.47796	5.91429	3.78518	77	chr1	3456586	3476586	TSS_ENSMUST00000161581_Gm1992	-1	+	Gm1992
		#0	chr1	
		#1	3460841	
		#2	3460939	
		#3	ESC.8_peak_11	
		#4	37	
		#5	.	
		#6	4.47796	
		#7	5.91429	
		#8	3.78518	
		#9	77	
		#10	chr1	
		#11	3456586	
		#12	3476586	
		#13	TSS_ENSMUST00000161581_Gm1992	
		#14	-1	
		#15	+	
		#16	Gm1992
		parts = l.strip().split('\t')
		peak = parts[3]
		gene = parts[-1]
		genes = peak2gene.get(peak,set([]))
		genes.add(gene)
		peak2gene[peak] = genes
	f.close()
	return peak2gene

def readMot2Peak(inname):
	mot2peak = []
	f = open(inname,'r')
	for l in f:
		#chr1	3007524	3007532	ESC.8_peak_1	20	.	3.48286	4.05357	2.09937	23	chr1	3007524	3007532	M6466_1.028.776775	+
		#0	chr1	
		#1	3007524	
		#2	3007532	
		#3	ESC.8_peak_1	
		#4	20	
		#5	.	
		#6	3.48286	
		#7	4.05357	
		#8	2.09937	
		#9	23	
		#10	chr1	
		#11	3007524	
		#12	3007532	
		#13	M6466_1.02
		#14	8.776775	
		#15	+
		parts = l.strip().split('\t')
		peak = parts[3]
		mot  = parts[-3] #parts[13]
		val  = float(parts[-2])  # parts[14]
		mot2peak.append((mot,peak,val))
	f.close()
	return mot2peak

def makeNet(mot2peak,mot2tf,peak2gene):
	net = {}
	for (mot,peak,v) in mot2peak:
		if mot not in mot2tf:
			continue
		if peak not in peak2gene:
			continue
		tfs = mot2tf[mot]
		genes = peak2gene[peak]
		for tf in tfs:
			for gene in genes:
				if (tf,gene) not in net:
					net[(tf,gene)] = v
				else:
					if v > net[(tf,gene)]:
						net[(tf,gene)] = v
	return net

def writeNet(outname,net):
	f = open(outname,'w')
	for (tf,gene) in net:
		f.write('%s\t%s\t%f\n' % (tf,gene,net[(tf,gene)]))
	f.close()

if __name__ == '__main__':
	mot2tf = readMot2TF(args.mot2tf)
	mot2peak = readMot2Peak(args.mot2peak)
	peak2gene = readPeak2Gene(args.peak2gene)
	net = makeNet(mot2peak,mot2tf,peak2gene)
	writeNet(args.outfile,net)
