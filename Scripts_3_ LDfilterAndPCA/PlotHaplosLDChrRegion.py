#!/usr/bin/env python3
# _*_ coding: Utf-8 _*_
# coding: utf-8

import pandas as pd
import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib import collections  as mc
from matplotlib import colors
#from matplotlib import figure as figure
import re
import csv
from collections import defaultdict

# Import data
#Import chromosome lengths
chroms = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/StandardFiles/HAv3_1_Chromosomes.list", sep = "\t", header = None)
chroms.columns = ['Accession','Length']
chroms['chrNb'] = chroms.index + 1
chroms['LengthMb'] = chroms['Length'] / 1000000

#List of Chromosomes and chromosome names
chromsList = list(range(1,17))
chromsList.reverse()
ybins = len(chromsList)

chromlengthList = chroms.LengthMb.tolist()[0:16]
chromlengthList.reverse()

chrNames = list(range(1,ybins+1))
chrNames = ["Chr " + str(s) for s in chrNames]
chrNames.reverse()

#Import haploblocks data
haploBlocks = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/HaplotypeBlocksPlink/haplotypeBlocks5000.blocks.cols", sep = " ")
haploBlocks['BP'] = haploBlocks['KB'] * 1000

#Import SNP selected on LD

LD02Chr = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/PCA7millionSNPsWindowChr/bimFiles/LD02_Chromosomes_pruned.bim", sep = "\t", header=None)
LD02Chr.columns = ['Chr','SNP','Val','Pos','Ref','Alt']
LD02Chr['PosMb'] = LD02Chr['Pos'] / 1000000

LD05Chr = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/PCA7millionSNPsWindowChr/bimFiles/LD05_Chromosomes_pruned.bim", sep = "\t", header=None)
LD05Chr.columns = ['Chr','SNP','Val','Pos','Ref','Alt']
LD05Chr['PosMb'] = LD05Chr['Pos'] / 1000000

LD09Chr = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/PCA7millionSNPsWindowChr/bimFiles/LD09_Chromosomes_pruned.bim", sep = "\t", header=None)
LD09Chr.columns = ['Chr','SNP','Val','Pos','Ref','Alt']
LD09Chr['PosMb'] = LD09Chr['Pos'] / 1000000

LD02kb100 = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/PCA7millionSNPsWindow100kb/bimFiles/LD02_Chromosomes_pruned.bim", sep = "\t", header=None)
LD02kb100.columns = ['Chr','SNP','Val','Pos','Ref','Alt']
LD02kb100['PosMb'] = LD02kb100['Pos'] / 1000000

LD05kb100 = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/PCA7millionSNPsWindow100kb/bimFiles/LD05_Chromosomes_pruned.bim", sep = "\t", header=None)
LD05kb100.columns = ['Chr','SNP','Val','Pos','Ref','Alt']
LD05kb100['PosMb'] = LD05kb100['Pos'] / 1000000

LD09kb100 = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/PCA7millionSNPsWindow100kb/bimFiles/LD09_Chromosomes_pruned.bim", sep = "\t", header=None)
LD09kb100.columns = ['Chr','SNP','Val','Pos','Ref','Alt']
LD09kb100['PosMb'] = LD09kb100['Pos'] / 1000000

LD02kb10 = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/PCA7millionSNPsWindow10kb/bimFiles/LD02_Chromosomes_pruned.bim", sep = "\t", header=None)
LD02kb10.columns = ['Chr','SNP','Val','Pos','Ref','Alt']
LD02kb10['PosMb'] = LD02kb10['Pos'] / 1000000

LD05kb10 = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/PCA7millionSNPsWindow10kb/bimFiles/LD05_Chromosomes_pruned.bim", sep = "\t", header=None)
LD05kb10.columns = ['Chr','SNP','Val','Pos','Ref','Alt']
LD05kb10['PosMb'] = LD05kb10['Pos'] / 1000000

LD09kb10 = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/PCA7millionSNPsWindow10kb/bimFiles/LD09_Chromosomes_pruned.bim", sep = "\t", header=None)
LD09kb10.columns = ['Chr','SNP','Val','Pos','Ref','Alt']
LD09kb10['PosMb'] = LD09kb10['Pos'] / 1000000

"""
Edit the chromosome and region below
"""
# Select region
chrom = 11
Start = 3
End = 8

"""
End edit
"""

LD02ChrChr = LD02Chr[(LD02Chr.Chr == chrom) & (LD02Chr.PosMb > Start) & (LD02Chr.PosMb < End)]
LD05ChrChr = LD05Chr[(LD05Chr.Chr == chrom) & (LD05Chr.PosMb > Start) & (LD05Chr.PosMb < End)]
LD09ChrChr = LD09Chr[(LD09Chr.Chr == chrom) & (LD09Chr.PosMb > Start) & (LD09Chr.PosMb < End)]

LD02kb100Chr = LD02kb100[(LD02kb100.Chr == chrom) & (LD02kb100.PosMb > Start) & (LD02kb100.PosMb < End)]
LD05kb100Chr = LD05kb100[(LD05kb100.Chr == chrom) & (LD05kb100.PosMb > Start) & (LD05kb100.PosMb < End)]
LD09kb100Chr = LD09kb100[(LD09kb100.Chr == chrom) & (LD09kb100.PosMb > Start) & (LD09kb100.PosMb < End)]

LD02kb10Chr = LD02kb10[(LD02kb10.Chr == chrom) & (LD02kb10.PosMb > Start) & (LD02kb10.PosMb < End)]
LD05kb10Chr = LD05kb10[(LD05kb10.Chr == chrom) & (LD05kb10.PosMb > Start) & (LD05kb10.PosMb < End)]
LD09kb10Chr = LD09kb10[(LD09kb10.Chr == chrom) & (LD09kb10.PosMb > Start) & (LD09kb10.PosMb < End)]

#Plotting

binWidth = 0.1
nbBins = math.ceil(chromlengthList[chrom] / binWidth)

fig, ax = plt.subplots(figsize=(16,8))

#ax.set_facecolor("whitesmoke")

#plot haplotype blocks
haploBlocksChr = haploBlocks[(haploBlocks.CHR == chrom) & (haploBlocks.KB >= 100)]
for index, row  in haploBlocksChr.iterrows():
    start = row["BP1"] / 1000000
    end = row["BP2"] / 1000000
    ax.axvspan(start, end, ymin=0.0, ymax=1.0, alpha=0.99, color='lightgrey')

kwargs = dict(alpha = 1, bins=nbBins, density=False, stacked=True, histtype='step', linewidth=2)

ax.hist(LD09ChrChr['PosMb'], **kwargs, color = 'darkblue')
ax.hist(LD05ChrChr['PosMb'], **kwargs, color = 'blue')
ax.hist(LD02ChrChr['PosMb'], **kwargs, color = 'dodgerblue')

ax.hist(LD09kb100Chr['PosMb'], **kwargs, color = 'darkred')
ax.hist(LD05kb100Chr['PosMb'], **kwargs, color = 'red')
ax.hist(LD02kb100Chr['PosMb'], **kwargs, color = 'darkorange')

ax.hist(LD09kb10Chr['PosMb'], **kwargs, color = 'darkgreen')
ax.hist(LD05kb10Chr['PosMb'], **kwargs, color = 'limegreen')
ax.hist(LD02kb10Chr['PosMb'], **kwargs, color = 'lime')

ax.set_xlim(Start, End)

classes = ['LD 0.9: window = chromosome','LD 0.5: window = chromosome','LD 0.2: window = chromosome', \
            'LD 0.9: window = 100 kb','LD 0.5: window = 100 kb','LD 0.2: window = 100 kb', \
            'LD 0.9: window = 10 kb','LD 0.5: window = 10 kb','LD 0.2: window = 10 kb']
class_colours= ['darkblue','blue','dodgerblue','darkred','red','darkorange','darkgreen','limegreen','lime']

patches = [ plt.plot([],[], marker="o", ms=10, ls="", mec=None, color=class_colours[i],
            label="{:s}".format(classes[i]) )[0]  for i in range(len(classes)) ]
ax.legend(handles=patches,
           loc='upper center', ncol=1, facecolor="white", numpoints=1, fontsize = 10)

plt.xlabel("Chromosome " + str(chrom) + " (Mb)",fontsize=20)
plt.ylabel("Number of SNPs", fontsize = 20)
ax.set_title("Effect of window size on LD pruning in haplotype blocks",fontsize=30)

plt.savefig('/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/HaplotypeBlocksPlink/PlotChrs11_3_8_Mb.pdf', format="pdf")#, dpi=600)
