# SeqApiPop analyses: PCA

## Filters:

* --maf filters out all variants with minor allele frequency below the provided threshold (default 0.01)
* --geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed
* --mind does the same for samples.
* --indep-pairwise 500000 50000 0.9 : LD filters out variants with LD > 0.9, on a window of 500000 SNPs, a sliding window of 50000.


File with the 7023689  variants:
/work/project/cytogen/Alain/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/plinkAnalyses/MetaGenotypesCalled870_raw_snps_allfilter_plink


More stringent on missing data in individuals:


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

* 11075 variants removed due to missing genotype data (--geno)
* 15 samples removed due to missing genotype data (--mind).

**From the *.imiss plink file:**

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




### LD filters : Use a window of the size of the largest chromosome ~ 1,000,000 SNPs

#### LD = 0.9

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

#### LD = 0.5

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

#### LD = 0.3

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

#### LD = 0.2

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

#### LD = 0.15

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

#### LD = 0.1

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

#### LD = 0.05

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
