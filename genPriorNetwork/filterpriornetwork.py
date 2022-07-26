#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 21:56:10 2021
Shilu
"""

import sys
import argparse

parser = argparse.ArgumentParser('filter regulators and genes in the prior network')
parser.add_argument("--regfile",    type=str, help="regulator list file",        required=True)
parser.add_argument("--genefile",  type=str, help="gene list file",      required=True)
parser.add_argument("--netfile", type=str, help="network file that needs filtering",       required=True)
parser.add_argument("--outfile",   type=str, help="Output (TF to Gene) filtered network", required=True)
args = parser.parse_args()

def readregulators(regfile):
    regulators = set()
    f = open(regfile,'r')
    for l in f:
        parts = l.strip()
        regulators.add(parts)
    f.close()
    return regulators

def readgenes(genefile):
    genes=set()
    f = open(genefile,'r')
    for l in f:
        parts = l.strip()
        genes.add(parts)
    f.close()
    return genes


def filterNet(netfile,regulators,genes):
    network = {}
    f = open(netfile,'r')
    for l in f:
        parts = l.strip().split('\t')
        tf=parts[0]
        gene=parts[1]
        weight=float(parts[2])
        if tf not in regulators:
            continue
        if gene not in genes:
            continue
        network[(tf,gene)]=weight
    f.close()
    return network


def writeNet(outname,network):
    f = open(outname,'w')
    for (tf,gene) in network:
        f.write('%s\t%s\t%f\n' % (tf,gene,network[(tf,gene)]))
    f.close()

if __name__ == '__main__':
    regulators = readregulators(args.regfile)
    genes = readgenes(args.genefile)
    network = filterNet(args.netfile,regulators,genes)
    writeNet(args.outfile,network)



