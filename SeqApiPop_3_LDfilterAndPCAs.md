# SeqApiPop analyses: filtering on LD and PCA

<!-- TOC depthFrom:2 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [1. Defining haplotype blocks](#1-defining-haplotype-blocks)
	- [1.1 Obtaining 5 first columns for further Analysis](#11-obtaining-5-first-columns-for-further-analysis)
- [2. Preleminary filters on missing genotype data for SNPs and for samples:](#2-preleminary-filters-on-missing-genotype-data-for-snps-and-for-samples)
- [3. Filters on linkage desiquilibrium (LD), using different LD values and windows in plink and PCA](#3-filters-on-linkage-desiquilibrium-ld-using-different-ld-values-and-windows-in-plink-and-pca)
	- [3.1. PCA on complete dataset](#31-pca-on-complete-dataset)
	- [3.2. PCAs after LD filters : Use a window of the size of the largest chromosome ~ 1,000,000 SNPs](#32-pcas-after-ld-filters-use-a-window-of-the-size-of-the-largest-chromosome-1000000-snps)
		- [3.2.1. LD = 0.9](#321-ld-09)
		- [3.2.2. LD = 0.5](#322-ld-05)
		- [3.2.3. LD = 0.3](#323-ld-03)
		- [3.2.4. LD = 0.2](#324-ld-02)
		- [3.2.5. LD = 0.15](#325-ld-015)
		- [3.2.6. LD = 0.1](#326-ld-01)
		- [3.2.7. LD = 0.05](#327-ld-005)
	- [3.3. The same LD values were also run with smaller size windows](#33-the-same-ld-values-were-also-run-with-smaller-size-windows)
		- [3.3.1. For instance 1 Mb window and LD = 0.05:](#331-for-instance-1-mb-window-and-ld-005)
	- [3.4. Figure showing the number of SNPs selected for different combinations of LD thresholds and window sizes.](#34-figure-showing-the-number-of-snps-selected-for-different-combinations-of-ld-thresholds-and-window-sizes)
		- [3.4.1 Make a table fot the plot](#341-make-a-table-fot-the-plot)
		- [Plot data](#plot-data)
	- [3.5. Figures showing the influence of LD thresholds and window sizes on the selection of SNPs in haplotype blocks](#35-figures-showing-the-influence-of-ld-thresholds-and-window-sizes-on-the-selection-of-snps-in-haplotype-blocks)

<!-- /TOC -->

## 1. Defining haplotype blocks

```bash
#! /bin/bash

#haplotypeBlocks5000.sh

module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module

NAME=haplotypeBlocks5000

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
  --out ${NAME} \
  --blocks no-pheno-req no-small-max-span \
  --blocks-max-kb 5000
```

### 1.1 Obtaining 5 first columns for further Analysis
```bash
 awk '{print $1,$2,$3,$4,$5}' haplotypeBlocks5000.blocks.det > haplotypeBlocks5000.blocks.cols
```

## 2. Preleminary filters on missing genotype data for SNPs and for samples:

* --maf filters out all variants with minor allele frequency below the provided threshold (default 0.01)
* --geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed
* --mind does the same for samples.
* --indep-pairwise 500000 50000 0.9 : LD filters out variants with LD > 0.9, on a window of 500000 SNPs, a sliding window of 50000.


File with the 7023689  variants:
/work/project/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/plinkAnalyses/MetaGenotypesCalled870_raw_snps_allfilter_plink


More stringent on missing data in individuals:

```bash
#! /bin/bash

#convertToBed.bash

module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module

VCFin=/work/project/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/plinkAnalyses/MetaGenotypesCalled870_raw_snps_allfilter_plink.vcf
VCFout=/work/project/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/plinkAnalyses/MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno
plink --vcf ${VCFin} \
  --keep-allele-order \
  --a2-allele ${VCFin} 4 3 '#' \
  --allow-no-sex \
  --allow-extra-chr \
  --chr-set 16 \
  --set-missing-var-ids @:#[HAV3.1]\$1\$2 \
  --chr 1-16 \
  --mind 0.1 \
  --geno 0.1 \
  --out ${VCFout} \
  --make-bed \
  --missing
```

* 11075 variants removed due to missing genotype data (--geno)
* 15 samples removed due to missing genotype data (--mind).

**Samples removed: frequency of missing genotypes from the \*.imiss plink file:**

|ID | N_MISS | N_GENO | F_MISS|
|:---|---:|---:|---:|
|AOC4 | 1270404 | 7023689 | 0.1809|
|BR12 | 1269182 | 7023689 | 0.1807|
|BR1A | 1253619 | 7023689 | 0.1785|
|ESP9 | 6208279 | 7023689 | 0.8839|
|JFM21 | 725846 | 7023689 | 0.1033|
|JFM24 | 817509 | 7023689 | 0.1164|
|JFM3 | 875208 | 7023689 | 0.1246|
|JFM5 | 830181 | 7023689 | 0.1182|
|KF21 | 722607 | 7023689 | 0.1029|
|OUE8 | 831427 | 7023689 | 0.1184|
|PM1 | 969888 | 7023689 | 0.1381|
|SavB1 | 823422 | 7023689 | 0.1172|
|SavB3 | 706024 | 7023689 | 0.1005|
|XC3 | 821334 | 7023689 | 0.1169|
|XC4 | 747325 | 7023689 | 0.1064|


## 3. Filters on linkage desiquilibrium (LD), using different LD values and windows in plink and PCA


### 3.1. PCA on complete dataset

```bash
#! /bin/bash

#7millionSNPs.sh

module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module

NAME=7millionSNPs

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
  --out PCA_${NAME} \
  --pca
```

### 3.2. PCAs after LD filters : Use a window of the size of the largest chromosome ~ 1,000,000 SNPs

#### 3.2.1. LD = 0.9

```bash
#! /bin/bash

#LD09_Chromosomes.sh

module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module

NAME=LD09_Chromosomes

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
  --out ${NAME} \
  --indep-pairwise 1000000 100000 0.9

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
    --out ${NAME}_pruned \
	--extract ${NAME}.prune.in \
	--make-bed

plink --bfile ${NAME}_pruned \
  --out PCA_${NAME} \
  --pca
```

#### 3.2.2. LD = 0.5

```bash
#! /bin/bash

#LD05_Chromosomes.sh

module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module

NAME=LD05_Chromosomes

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
  --out ${NAME} \
  --indep-pairwise 1000000 100000 0.5

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
    --out ${NAME}_pruned \
	--extract ${NAME}.prune.in \
	--make-bed

plink --bfile ${NAME}_pruned \
  --out PCA_${NAME} \
  --pca
```

#### 3.2.3. LD = 0.3

```bash
#! /bin/bash

#LD03_Chromosomes.sh

module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module

NAME=LD03_Chromosomes

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
  --out ${NAME} \
  --indep-pairwise 1000000 100000 0.3

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
    --out ${NAME}_pruned \
	--extract ${NAME}.prune.in \
	--make-bed

plink --bfile ${NAME}_pruned \
  --out PCA_${NAME} \
  --pca
```

#### 3.2.4. LD = 0.2

```bash
#! /bin/bash

#LD02_Chromosomes.sh

module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module

NAME=LD02_Chromosomes

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
  --out ${NAME} \
  --indep-pairwise 1000000 100000 0.2

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
    --out ${NAME}_pruned \
	--extract ${NAME}.prune.in \
	--make-bed

plink --bfile ${NAME}_pruned \
  --out PCA_${NAME} \
  --pca
```

#### 3.2.5. LD = 0.15

```bash
#! /bin/bash

#LD015_Chromosomes.sh

module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module

NAME=LD015_Chromosomes

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
  --out ${NAME} \
  --indep-pairwise 1000000 100000 0.15

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
    --out ${NAME}_pruned \
	--extract ${NAME}.prune.in \
	--make-bed

plink --bfile ${NAME}_pruned \
  --out PCA_${NAME} \
  --pca
```

#### 3.2.6. LD = 0.1

```bash
#! /bin/bash

#LD01_Chromosomes.sh

module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module

NAME=LD01_Chromosomes

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
  --out ${NAME} \
  --indep-pairwise 1000000 100000 0.1

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
    --out ${NAME}_pruned \
	--extract ${NAME}.prune.in \
	--make-bed

plink --bfile ${NAME}_pruned \
  --out PCA_${NAME} \
  --pca
```

#### 3.2.7. LD = 0.05

```bash
#! /bin/bash

#LD005_Chromosomes.sh

module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module

NAME=LD005_Chromosomes

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
  --out ${NAME} \
  --indep-pairwise 1000000 100000 0.05

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
    --out ${NAME}_pruned \
	--extract ${NAME}.prune.in \
	--make-bed

plink --bfile ${NAME}_pruned \
  --out PCA_${NAME} \
  --pca
```

### 3.3. The same LD values were also run with smaller size windows
* 1 Mb
* 500 kb
* 100 kb
* 50 kb
* 10 kb

#### 3.3.1. For instance 1 Mb window and LD = 0.05:

```bash
#! /bin/bash

#LD005_Chromosomes.sh

module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module

NAME=LD005_Chromosomes

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
  --out ${NAME} \
  --indep-pairwise 1000 kb 100 0.05

plink --bfile ../MetaGenotypesCalled870_raw_snps_allfilter_plink_missIndGeno \
    --out ${NAME}_pruned \
        --extract ${NAME}.prune.in \
        --make-bed

plink --bfile ${NAME}_pruned \
  --out PCA_${NAME} \
  --pca
```

### 3.4. Figure showing the number of SNPs selected for different combinations of LD thresholds and window sizes.
#### 3.4.1 Make a table fot the plot

* Table : LD threshold, window size, SNPs retained
* Done by extracting the data from the logfiles of the plink --pca runs

```bash
#!/bin/bash
grep 'variants loaded' PCAsWind*/PCA*log | sed s'/variants loaded from .bim file.//' | sed 's/:/\t/' | sed 's/_/\t/g' | sed 's/\//\t/' | sed 's/PCAsWindow//' | grep -v 7millionSNPs
 | awk 'BEGIN{OFS="\t"}{print $1,$3,$5}' > SNPselectionStats2.txt
```

#### Plot data
* script PlotHaplosLDChrRegion.py

![SNPs selected by LD values and window sizes.](/Figures_3_LD_PCA/NbSNPsFunctionLD_Window.png)
![downloadable pdf version](/Figures_3_LD_PCA/NbSNPsFunctionLD_Window.pdf)

### 3.5. Figures showing the influence of LD thresholds and window sizes on the selection of SNPs in haplotype blocks
* Areas shaded in grey are the haplotype blocks > 100 as detected by plink.
* Whole genome: ![SNP numbers along the genome after LD pruning.](Figures_3_LD_PCA/PlotAllChrsCompareWindBin05.pdf)
  - Figure generated with PlotHaplosLDselCompWind.py
* And a selection of 4 haplotype haploBlocks
  - Figures generated with PlotHaplosLDChrRegion.py

![Chromosome 2](/Figures_3_LDfilterAndPCA/PlotChrs2_0_3_Mb.pdf)
![Chromosome 4](/Figures_3_LDfilterAndPCA/PlotChrs4_0_3_Mb.pdf)
![Chromosome 7](/Figures_3_LDfilterAndPCA/PlotChrs7_3_8_Mb.pdf)
![Chromosome 11](/Figures_3_LDfilterAndPCA/PlotChrs11_3_8_Mb.pdf)

* SNP counts in large haplotypeblocks are lower.
* Larger window sizes for LD pruning eliminate more SNPs in large haplotype blocks
