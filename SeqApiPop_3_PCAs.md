# SeqApiPop analyses: PCA

## Filters:

* --maf filters out all variants with minor allele frequency below the provided threshold (default 0.01)
* --geno filters out all variants with missing call rates exceeding the provided value (default 0.1) to be removed
* --mind does the same for samples.
* --indep-pairwise 500000 50000 0.9 : LD filters out variants with LD > 0.9, on a window of 500000 SNPs, a sliding window of 50000.



Use a window of the size of the largest chromosome ~ 500,000 SNPs

NC_001566.1     287
NC_037638.1     913023
NC_037639.1     534733
NC_037640.1     442882
NC_037641.1     440141
NC_037642.1     462122
NC_037643.1     577596
NC_037644.1     463575
NC_037645.1     397891
NC_037646.1     378566
NC_037647.1     355296
NC_037648.1     441395
NC_037649.1     388488
NC_037650.1     378466
NC_037651.1     330298
NC_037652.1     285750
NC_037653.1     233467
Sum     7023976
