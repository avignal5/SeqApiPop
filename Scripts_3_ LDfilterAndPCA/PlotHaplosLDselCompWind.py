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

"""
Import chromosome lengths. In the present case, the file also used for the mapping and other scripts, with:
column 1 = accession number for the chromosomes, as found in the vcf files
column 2 = chromosome length in bp
"""
#Import chromosome lengths
chroms = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/StandardFiles/HAv3_1_Chromosomes.list", sep = "\t", header = None)
chroms.columns = ['Accession','Length']
chroms['chrNb'] = chroms.index + 1
chroms['LengthMb'] = chroms['Length'] / 1000000

#List of Chromosomes, as numbers such as in plink
chromsList = list(range(1,17))
chromsList.reverse()
ybins = len(chromsList)

#List of the corresponding chromosome lengths
chromlengthList = chroms.LengthMb.tolist()[0:16]
chromlengthList.reverse()

#List of chromosomes as Chr&, Chr2 ... Chr16, for the figure legend
chrNames = list(range(1,ybins+1))
chrNames = ["Chr " + str(s) for s in chrNames]
chrNames.reverse()

"""
Import the haploblock data, from a file corresponding to the 5 first columns of the output file "haplotypeBlocks5000.blocks.det"
from the plink "--block" function
To highlight the haploblocks of a size larger than the setting in the "plot haplotype blocks" section below
"""
#Import haploblocks data
haploBlocks = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/HaplotypeBlocksPlink/haplotypeBlocks5000.blocks.cols", sep = " ")
haploBlocks['BP'] = haploBlocks['KB'] * 1000

"""
Import SNP selected on LD: the *.bim outputs corresponding to the various parameters used in plink "--indep-pairwise" function
"""
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
Figure
"""

#Figure title
title = "Influence of window size on LD selection"
#Text
comments = ""

#Figure settings
space = (1 / (ybins + 1)) * 1 / 6
figHeight = (ybins + 1) * 3
figWidth = (max(chromlengthList)) * 1.4 # + space * 4

fig = plt.figure(figsize=(figWidth, figHeight))
#plt.axes([0, 10, 0, 10])
fig.suptitle(title, fontsize=96)
fig.text(.5,.2,comments, wrap=True, fontsize=24)
fig.text(.35,.01, "Position (Mb)", fontsize=40)
fig.text(.01,.5, "SNP density", fontsize=40, rotation='vertical',)

bottom = space

#Set bins step_size
binWidth = 0.05

#Plot chromosomes
count=0
images=[]
for chr_i in chromsList:
	#Define size and draw box
	#left = space * 4
	left = space * 7
	width = (chromlengthList[count]) * 1.4 /figWidth * 0.90
	height = (1/ybins) * 4.5/6
	ax = plt.axes([left, bottom, width, height])
	bottom = bottom + height + space

	#data
	LD02ChrChr = LD02Chr[(LD02Chr.Chr == chr_i)]
	LD05ChrChr = LD05Chr[(LD05Chr.Chr == chr_i)]
	LD09ChrChr = LD09Chr[(LD09Chr.Chr == chr_i)]

	LD02kb100Chr = LD02kb100[(LD02kb100.Chr == chr_i)]
	LD05kb100Chr = LD05kb100[(LD05kb100.Chr == chr_i)]
	LD09kb100Chr = LD09kb100[(LD09kb100.Chr == chr_i)]

	LD02kb10Chr = LD02kb10[(LD02kb10.Chr == chr_i)]
	LD05kb10Chr = LD05kb10[(LD05kb10.Chr == chr_i)]
	LD09kb10Chr = LD09kb10[(LD09kb10.Chr == chr_i)]

	#Histogram bins
	nbBins = math.ceil(chromlengthList[count] / binWidth)

	#plot SNP densities
	kwargs = dict(alpha = 1, bins=nbBins, density=False, stacked=True, histtype='step')

	images.append(ax.hist(LD09ChrChr['PosMb'], **kwargs, color = 'darkblue'))
	images.append(ax.hist(LD05ChrChr['PosMb'], **kwargs, color = 'blue'))
	images.append(ax.hist(LD02ChrChr['PosMb'], **kwargs, color = 'dodgerblue'))

	images.append(ax.hist(LD09kb100Chr['PosMb'], **kwargs, color = 'darkred'))
	images.append(ax.hist(LD05kb100Chr['PosMb'], **kwargs, color = 'red'))
	images.append(ax.hist(LD02kb100Chr['PosMb'], **kwargs, color = 'darkorange'))

	images.append(ax.hist(LD09kb10Chr['PosMb'], **kwargs, color = 'darkgreen'))
	images.append(ax.hist(LD05kb10Chr['PosMb'], **kwargs, color = 'limegreen'))
	images.append(ax.hist(LD02kb10Chr['PosMb'], **kwargs, color = 'lime'))

	ax.set_xlim(-0.1, chromlengthList[count] +0.1)

	#Labels, ticks, grids
	ax.tick_params(axis='y', labelsize=15)
	ax.set_ylabel(chrNames[count], fontsize=30)
	major_ticks = np.arange(0,chromlengthList[count],1,dtype='int16')
	ax.set_xticks(major_ticks)
	ax.set_xticklabels(major_ticks, fontsize=24)
	#ax.set_yticklabels(, fontsize=24)
	minor_ticks = np.arange(0,chromlengthList[count],0.250000)
	ax.set_xticks(minor_ticks, minor=True)
	ax.grid(which='minor', alpha=0.2)
	ax.grid(which='major', alpha=0.5)

	#plot haplotype blocks
	haploBlocksChr = haploBlocks[(haploBlocks.CHR == chr_i) & (haploBlocks.KB >= 100)]
	for index, row  in haploBlocksChr.iterrows():
		start = row["BP1"] / 1000000
		end = row["BP2"] / 1000000
		ax.axvspan(start, end, ymin=0.0, ymax=1.0, alpha=0.3, color='grey')

	count = count + 1

#plt.figlegend([a, b, c, d, e], ['AllSNPs','LD 0.9', 'LD 0.5', 'LD 0.3', 'LD 0.1'], loc=(0.85, 0.65))

classes = ['LD 0.9: window = chromosome','LD 0.5: window = chromosome','LD 0.2: window = chromosome', \
			'LD 0.9: window = 100 kb','LD 0.5: window = 100 kb','LD 0.2: window = 100 kb', \
			'LD 0.9: window = 10 kb','LD 0.5: window = 10 kb','LD 0.2: window = 10 kb']
class_colours= ['darkblue','blue','dodgerblue','darkred','red','darkorange','darkgreen','limegreen','lime']

patches = [ plt.plot([],[], marker="o", ms=15, ls="", mec=None, color=class_colours[i],
            label="{:s}".format(classes[i]) )[0]  for i in range(len(classes)) ]
plt.figlegend(handles=patches,
           loc='lower right', ncol=1, facecolor="white", numpoints=1, fontsize = 36)


plt.savefig('/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/HaplotypeBlocksPlink/PlotAllChrsCompareWindBin05.pdf', format="pdf")#, dpi=600)
#plt.savefig('/Users/avignal/Documents/Stats/2016_PacificBee/TandemRepeatFinder/chrs.pdf', format="pdf")

#plt.show()
