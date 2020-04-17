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

selStats2 = pd.read_csv("/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/SNPselectionStats2.txt",header=None, sep="\t")
selStats2.columns=['Window','LDValue','NbSNPs']
selStats2['Nb_k_SNPs'] = selStats2['NbSNPs'] / 1000
selStats2 = selStats2.set_index(['LDValue','Window'])
selStats2 = selStats2.unstack('Window')
windows = selStats2.columns.levels[1].tolist()
LDvalues = selStats2.index.values.tolist()
barWidth = 0.15

y1 = np.arange(len(LDvalues)).tolist()
y2 = [x + barWidth for x in y1]
y3 = [x + barWidth for x in y2]
y4 = [x + barWidth for x in y3]
y5 = [x + barWidth for x in y4]
y6 = [x + barWidth for x in y5]

fig, ax = plt.subplots(figsize=(16,8))
plt.bar(y1,selStats2['Nb_k_SNPs']['10kb'].tolist(), width=barWidth, color = 'red', edgecolor='white', label='10 kb')
plt.bar(y2,selStats2['Nb_k_SNPs']['50kb'].tolist(), width=barWidth, color = 'blue', edgecolor='white', label='50 kb')
plt.bar(y3,selStats2['Nb_k_SNPs']['100kb'].tolist(), width=barWidth, color = 'green', edgecolor='white', label='100 kb')
plt.bar(y4,selStats2['Nb_k_SNPs']['500kb'].tolist(), width=barWidth, color = 'darkred', edgecolor='white', label='500 kb')
plt.bar(y5,selStats2['Nb_k_SNPs']['1Mb'].tolist(), width=barWidth, color = 'dodgerblue', edgecolor='white', label='1 Mb')
plt.bar(y6,selStats2['Nb_k_SNPs']['MillionSNPs'].tolist(), width=barWidth, color = 'limegreen', edgecolor='white', label='Million SNPs')

plt.xticks([r + 2.5 * barWidth for r in range(len(LDvalues))], LDvalues)
plt.legend(fontsize = 10, title = "Window", loc = 'upper left')

plt.xlabel("Linkage Desequilibrium threshold value", fontsize = 20)
plt.ylabel("Number of SNPs (X 1000)", fontsize = 20)
ax.set_title("SNPs selected by LD values and window sizes",fontsize=30)

fig.savefig("/Users/avignal/Documents/GitHub/SeqApiPop/Figures_3_LD_PCA/NbSNPsFunctionLD_Window.pdf", bbox_inches='tight')
