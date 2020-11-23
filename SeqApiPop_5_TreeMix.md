# SeqApiPop analyses: Treemix

<!-- TOC depthFrom:1 depthTo:6 withLinks:1 updateOnSave:1 orderedList:0 -->

- [SeqApiPop analyses: Treemix](#seqapipop-analyses-treemix)
	- [1. Files](#1-files)
		- [Lists for selection:](#lists-for-selection)
		- [All samples with > 80 pure backgrounds, plus Corsica:](#all-samples-with-80-pure-backgrounds-plus-corsica)
		- [All samples with > 90 pure backgrounds, plus Corsica:](#all-samples-with-90-pure-backgrounds-plus-corsica)
		- [Select samples with plinkAnalyses](#select-samples-with-plinkanalyses)
- [! /bin/bash](#-binbash)
- [! /bin/bash](#-binbash)

<!-- /TOC -->

## 1. Files

### Lists for selection:

Selection from the 629 samples, after Admixture analysis

### All samples with > 80 pure backgrounds, plus Corsica:

``bash
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

``bash
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

### Select samples with plinkAnalyses

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
