#!/usr/bin/env python3
# _*_ coding: Utf-8 _*_
# coding: utf-8


import numpy as np
import matplotlib.patches as mpatches
from matplotlib import colors
from matplotlib import cm
import re
import csv
from collections import defaultdict
from IPython.display import display
from IPython.display import set_matplotlib_formats
set_matplotlib_formats('png', 'pdf')

"""
Edit path; the 6 eigenvec files to plot; the output file name
"""
path = "/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/PCAsAll/"

# Read data description
dataSamples = pd.read_excel('/Users/avignal/Documents/Stats/2019_SeqApiPop_HAv3_1/ListsSamples/SequenceGroups870.xlsx',sheet_name='869Samples')
dataSamples = dataSamples.loc[:,['name','RepresentsHive','CategPCA']]
dataSamples = dataSamples.set_index('name')

#The 6 PCAs to plot
PCA1 = "PCA_MillionSNPs_LD005"
PCA2 = "PCA_1Mb_LD005"
PCA3 = "PCA_500kb_LD005"
PCA4 = "PCA_100kb_LD005"
PCA5 = "PCA_50kb_LD005"
PCA6 = "PCA_10kb_LD005"

outName = "LD005_AllWindows"

"""
End Edits
"""

"""
EIGENVECTORS DATABASE
"""
def importEigenvec(path,fileName):
    eigenvec = pd.read_csv(path + fileName + ".eigenvec", sep = " ", header=None)
    eigenvec.columns = ['name','name2','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','PC16','PC17','PC18','PC19','PC20']
    eigenvec = eigenvec.drop(columns=['name2'])
    eigenvec['Color'] = '#88888800'
    eigenvec['Selection'] = fileName
    eigenvec = eigenvec.set_index('name')
    data = eigenvec.join(dataSamples)
    return(data)

data1 = importEigenvec(path=path,fileName=PCA1)
data2 = importEigenvec(path=path,fileName=PCA2)
data3 = importEigenvec(path=path,fileName=PCA3)
data4 = importEigenvec(path=path,fileName=PCA4)
data5 = importEigenvec(path=path,fileName=PCA5)
data6 = importEigenvec(path=path,fileName=PCA6)

data = pd.concat([data1,data2,data3,data4,data5,data6])

"""
ASSIGN COLORS
"""
classes = []
class_colours = []

#melliferas
data.loc[data.CategPCA.str.contains('Ouessant',regex=True),'Color'] = col = '#000000ff'
classes.append("Mellifera Ouessant")
class_colours.append(col)
data.loc[data.CategPCA.str.contains('UK',regex=True),'Color'] = col  = '#1313ecff'
classes.append("Mellifera Colonsay")
class_colours.append(col)
data.loc[data.CategPCA.str.contains('Spain',regex=True),'Color'] = col  = '#4d79ffff'
classes.append("Mellifera Spain")
class_colours.append(col)
data.loc[data.CategPCA.str.contains('Savoy',regex=True),'Color'] = col  = '#1affffff'
classes.append("Mellifera Savoy")
class_colours.append(col)
data.loc[data.CategPCA.str.contains('Porquerolles',regex=True),'Color'] = col  = '#ac00e6ff'
classes.append("Mellifera Porquerolles")
class_colours.append(col)
data.loc[data.CategPCA.str.contains('Solliès',regex=True),'Color'] = col  = '#df80ffff'
classes.append("Mellifera Solliès")
class_colours.append(col)

#ligusticas and carnicas
data.loc[data.CategPCA.str.contains('Italy',regex=True),'Color'] = col  = '#ffff00ff'
classes.append("Ligustica Italy")
class_colours.append(col)
data.loc[data.CategPCA.str.contains('Slovenia',regex=True),'Color'] = col  = '#ffcc00ff'
classes.append("Carnica Slovenia")
class_colours.append(col)
data.loc[data.CategPCA.str.contains('CarGermany',regex=True),'Color'] = col  = '#ff9933ff'
classes.append("Carnica Germany")
class_colours.append(col)
data.loc[data.CategPCA.str.contains('CarSwitzerland',regex=True),'Color'] = col  = '#b35900'
classes.append("Carnica Switzerland")
class_colours.append(col)
data.loc[data.CategPCA.str.contains('CarFrance',regex=True),'Color'] = col  = '#663300ff'
classes.append("Carnica France")
class_colours.append(col)
data.loc[data.CategPCA.str.contains('Poland',regex=True),'Color'] = col  = '#8f8f3d'
classes.append("Carnica Poland")
class_colours.append(col)

#caucasicas
data.loc[data.CategPCA.str.contains('CauFrance',regex=True),'Color'] = col  = '#00ff55ff'
classes.append("Caucasica")
class_colours.append(col)

"""
EIGENVALUES DATABASE
"""
def importEigenval(path,fileName):
    eigenval = pd.read_csv(path + fileName + ".eigenval", sep = " ", header=None)
    eigenval.columns = ['Var']
    sumVar = eigenval.Var.sum()
    eigenval['PercentVar'] = eigenval['Var'] / sumVar * 100
    eigenval['Selection'] = fileName
    return(eigenval)

dataVal1 = importEigenval(path=path,fileName=PCA1)
dataVal2 = importEigenval(path=path,fileName=PCA2)
dataVal3 = importEigenval(path=path,fileName=PCA3)
dataVal4 = importEigenval(path=path,fileName=PCA4)
dataVal5 = importEigenval(path=path,fileName=PCA5)
dataVal6 = importEigenval(path=path,fileName=PCA6)

dataVal = pd.concat([dataVal1,dataVal2,dataVal3,dataVal4,dataVal5,dataVal6])

fig, axes = plt.subplots(2, 3, figsize=[36, 24])
y_pos = list(range(1,21))

def plotPCAs(axe,fileName):
    axe.scatter(x="PC1",y="PC2",data=data[(data.Selection == fileName)], color="Color")
    insetVals = inset_axes(axe,
                    width="20%",
                    height='20%',
                    loc='lower center')

    plt.bar(y_pos,dataVal[(dataVal.Selection == fileName)]['PercentVar'])
    plt.xticks([])
    plt.yticks([])

    #Window size and LD threshold
    window = fileName.split('_')[1]
    LDthresh = fileName.split('_')[2]

    #Obtain the number of SNPs
    with open(path + fileName + ".log") as csvFile:
        SNPnb = csv.reader(csvFile, delimiter = " ")
        for row in SNPnb:
            if row:
                if len(row) > 8 and re.search('filters',row[6]):
                    if window != "All":
                        axe.set_title("Window: " + window + ", threshold: " + LDthresh +"\n"+ f"{int(row[0]):,}" + " SNPs retained", fontsize=30)
                    else:
                        axe.set_title("All SNPs" + "\n"+ f"{int(row[0]):,}" + " SNPs retained", fontsize=30)

plotPCAs(axes[0,0],PCA1)
plotPCAs(axes[0,1],PCA2)
plotPCAs(axes[0,2],PCA3)
plotPCAs(axes[1,0],PCA4)
plotPCAs(axes[1,1],PCA5)
plotPCAs(axes[1,2],PCA6)

patches = [ plt.plot([],[], marker="o", ms=25, ls="", mec=None, color=class_colours[i],
    label="{:s}".format(classes[i]) )[0]  for i in range(len(classes)) ]
plt.figlegend(handles=patches,loc='lower center', ncol=5, facecolor="white", numpoints=1, fontsize = 25)

plt.savefig(path + outName + ".pdf", format = 'pdf')
