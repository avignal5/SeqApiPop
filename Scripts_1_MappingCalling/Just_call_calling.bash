#!/bin/bash

# Usage:	Just_call_calling.bash

N=1
PLOIDY=2

#Sometimes, the mapping scripts stop before the merging of the per-chromosome bootstraped bam files
#This script will call the script that does the merging, sorting and calling => *.g.vcf files
#The lines in the file samples_list are paths to the directories containing the mapping, etc directories

#------------------------------------------------------------------
#EDIT FOR ACCESS TO SEQUENCE FILES
# First download links from NG6 to genotoul and update the path here
SAMPLE_FILE=/work/genphyse/cytogen/Alain/seqapipopOnHAV3_1/MapSardinians/samples_list
#OUT=/work/genphyse/cytogen/Alain/seqapipopOnHAV3_1/MapSardinians
REF=/home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/GCF_003254395.2_Amel_HAv3.1_genomic.fna
CHROMOSOME=/home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/HAv3_1_Chromosomes.list
UNKNOWN=/home/gencel/vignal/save/Genomes/Abeille/HAv3_1_indexes/HAv3_1_Unknown.list
SCRIPT=/genphyse/cytogen/seqapipop/ScriptsHAV3Mapping/mapping_calling
#------------------------------------------------------------------
#mkdir -p ${OUT}/logs/mapping	# /work/ktabetaoul/Project_seqapipop/RESULTATS/logs/mapping
#cd ${OUT}

#EDIT FOR mapping.sh
# This is the loop to send each sample to the cluster to be mapped
while read line
do
	ID=`basename ${line}`	# give the name of the sample-population to map
#	IN1=`dirname ${line}`	# give the path where find the sample fastq
	OUT=${line}/calling
	IN=${line}/bootstraping
	mkdir -p ${OUT}/logs

	sbatch --cpus-per-task=1 --mem-per-cpu=6G \
		-J ${ID}_mapping -o ${OUT}/logs/${ID}_mapping.o -e ${OUT}/logs/${ID}_mapping.e \
		${SCRIPT}/callingAV_2019_Dec.sh -s ${ID} -i ${IN} -o ${OUT} -p ${PLOIDY} -n ${N} -R ${REF}
		
done < ${SAMPLE_FILE}


# end of file
