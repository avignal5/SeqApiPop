# SeqApiPop analyses: Treemix

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [SeqApiPop analyses: Treemix](#seqapipop-analyses-treemix)
	- [1. Files](#1-files)
		- [Lists for selection:](#lists-for-selection)
		- [All samples with > 80 pure backgrounds, plus Corsica:](#all-samples-with-80-pure-backgrounds-plus-corsica)
		- [All samples with > 90 pure backgrounds, plus Corsica:](#all-samples-with-90-pure-backgrounds-plus-corsica)
		- [Select samples](#select-samples)
- [! /bin/bash](#-binbash)
- [! /bin/bash](#-binbash)
		- [Calculate frequencies and format for TreeMix](#calculate-frequencies-and-format-for-treemix)
- [! /bin/bash](#-binbash)
		- [Run Treemix](#run-treemix)
- [!/bin/bash](#binbash)
		- [Estimate the number of migrations with OptM](#estimate-the-number-of-migrations-with-optm)

<!-- /TOC -->

## 1. Files

### Lists for selection:
```bash
Selection from the 629 samples, after Admixture analysis
```

### All samples with > 80 pure backgrounds, plus Corsica:

```bash
$ head SeqApiPop_390_9pops.list
AOC10 AOC10
AOC11 AOC11
AOC12 AOC12
AOC14 AOC14

$ head SeqApiPop_390_9pops.fam
Corsica AOC10 0 0 0 -9
Corsica AOC11 0 0 0 -9
Corsica AOC12 0 0 0 -9
Corsica AOC14 0 0 0 -9
```

### All samples with > 90 pure backgrounds, plus Corsica:

```bash
$ head SeqApiPop_324_9pops.list
AOC10 AOC10
AOC11 AOC11
AOC12 AOC12
AOC14 AOC14

$ head SeqApiPop_324_9pops.fam
AOC10 AOC10 0 0 0 -9
AOC11 AOC11 0 0 0 -9
AOC12 AOC12 0 0 0 -9
AOC14 AOC14 0 0 0 -9

```

### Select samples

selectRefPopInds90.bash

* All SNPs
```bash
#! /bin/bash
module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module
plink --bfile ../SeqApiPop_629 --out SeqApiPop_324_9pops_AllSNPs --keep SeqApiPop_324_9pops.list --make-bed
```

* SNPs selected on MAF and LD
```bash
#! /bin/bash
module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module
plink --bfile ../_maf001_LD03_prune --out SeqApiPop_324_9pops_AllSNPs --keep SeqApiPop_324_9pops.list --make-bed
```

### Calculate frequencies and format for TreeMix

calculateFrequencies90.bash

```bash
#! /bin/bash

module load -f /work/project/cytogen/Alain/seqapipopOnHAV3_AV/program_module

plink --bfile SeqApiPop_324_9pops --freq --missing --family --out SeqApiPop_324_9pops
gzip SeqApiPop_324_9pops.frq.strat
```

```bash
sbatch --mem=8G --wrap="..plink2treemix.py SeqApiPop_324_9pops.frq.strat.gz SeqApiPop_324_9pops.frq.gz"
```

### Run Treemix

For instance, in directory bootstrap90, for the 324 samples with > 90 % pure backgrounds

launchTreemix.bash

```bash
#!/bin/bash

module load bioinfo/treemix-1.13

for i in $(seq 0 9)
do
    for j in $(seq 0 99)
    do
        sbatch --wrap="treemix -i ../SeqApiPop_324_9pops.frq.gz -bootstrap -seed ${RANDOM} -k 500 -m ${i} -o outstemM${i}_rep${j}"
    done
done
```

### Estimate the number of migrations with OptM

```R
Bootstraps90.optm = optM("/Users/avignal/GenotoulBigWork/seqapipopOnHAV3_1/seqApiPopVcfFilteredSonia/plinkAnalyses/WindowSNPs/TreeMix/bootstraps90")

plot_optM(Bootstraps80.optm, method = "Evanno")
```
