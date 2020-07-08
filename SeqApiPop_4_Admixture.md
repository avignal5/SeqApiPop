# SeqApiPop analyses: admixture


```bash
#RefPopsAllSNP

PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
IN="/work/project/cytogen/Alain/metaAdmixture/AdmixtureZzDefinitive/MetaGenotypes_snps_OK_NoPL_NoIndel_620_qc2"
DIR="/work/project/cytogen/Alain/metaAdmixture/AdmixtureZzDefinitive/RefPopsAllSNP"
SAMPLES="RefPops.list"
OUT="RefPopsInds"


#awk '{print $1"\t"$1}' ${DIR}/${SAMPLES} > ${DIR}/list.temp #list for plink --keep (couldn't get it to work)

#awk '{print $1}' ${DIR}/${SAMPLES} > ${DIR}/list.temp

#bcftools view -O vcf -o ${DIR}/RefPopsInds.vcf -S ${DIR}/list.temp ${IN}.vcf

##${PLINK}/plink --file ${IN} \
#  --keep ${DIR}/list.temp \
#  --out ${DIR}/${OUT} \
#  --make-bed

##${PLINK}/plink --vcf ${DIR}/RefPopsInds.vcf \
#       --make-bed



for K in 1 2 3 4 5 6 7 8 9 ;
do
/usr/local/bioinfo/src/ADMIXTURE/admixture_linux-1.23/admixture --cv ${DIR}/${OUT}.bed ${K} | tee log${K}.out
done

#rm list.temp
```

## Other script
```bash
#RefPopsPlusSav

PLINK="/usr/local/bioinfo/src/plink/plink-1.90-x86_64"
IN="/work/project/cytogen/Alain/metaAdmixture/AdmixtureZzDefinitive/MetaGenotypes_snps_OK_NoPL_NoIndel_619_qc2_1KLD01"
DIR="/work/project/cytogen/Alain/metaAdmixture/AdmixtureZzDefinitive/RefPopsPlusSav"
SAMPLES="RefPopsSav.list"
OUT="RefPopsSavInds"

awk '{print $1"\t"$1}' ${DIR}/${SAMPLES} > ${DIR}/list.temp

${PLINK}/plink --file ${IN} \
  --keep ${DIR}/list.temp \
  --out ${DIR}/${OUT} \
  --make-bed

for K in 1 2 3 4 5 6 7 8 9 ;
do
/usr/local/bioinfo/src/ADMIXTURE/admixture_linux-1.23/admixture --cv ${DIR}/${OUT}.bed ${K} | tee log${K}.out
done

rm list.temp
```


## New script for slurm

### Analyses on LD03 window = chromosome

```bash
#!/bin/bash

#admixtureAnalysis.sh

module load -f  /home/gencel/vignal/save/000_ProgramModules/program_module

IN=SeqApiPop_628_MillionSNPs_LD03.bed


for K in 1 2 3 4 5 6 7 8 9 ;
do
sbatch --cpus-per-task=1 --mem-per-cpu=4G \
    -J ${K}admixt -o ${IN}.o -e ${IN}.e \
    --wrap="admixture --cv ../plinkFilesAll/${IN} ${K} | tee ${IN}.log${K}"
done
```
