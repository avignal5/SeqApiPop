# Update on the filtering scripts
The scripts here are those that were used in the paper:
Wragg, D., Eynard, S. E., Basso, B., Canale-Tabet, K., Labarthe, E., Bouchez, O., Bienefeld, K., Bieńkowska, M., Costa, C., Gregorc, A., Kryger, P., Parejo, M., Pinto, M. A., Bidanel, J.-P., Servin, B., Le Conte, Y., & Vignal, A. (2021). Complex population structure and haplotype patterns in Western Europe honey bee from sequencing a large panel of haploid drones. Molecular Ecology Resources. 2022 Jun 27;1755-0998.13665. https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13665

For an updated version, see https://github.com/seynard/vcf_cleanup


## R packages required
The scripts will work on a slurm cluster. In theory, they will install the required R packages if missing, but this will not always work, according to the version of R used (see module load in the run_vcfcleanup_3.sh script). The best is to install the required packages ('data.table','VennDiagram','reshape2','RColorBrewer','grDevices','ggplot2','viridis') before running the script.

## To edit in run_vcfcleanup_3.sh
* run='diagnostic'
Will only make the plots
* run='filter_all'
Will produce a vcf file once all filters are passed
* run='filter_sequential'
Will produce a vcf file at each filtering stage (not usually useful)

oui oui les ~ entre gorupes et _ intra groupe
## Grouping the filters
* in:
- kept_above_threshold="MQ_QUAL_QD ~GQ~GQ"
- kept_below_threshold="FS_SOR_allele ~miss_het ~GQfiltered"
* Filter joined by a _ are done together. Groups of filters are separated by ~
