#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mail-user tprieto@nygenome.org
#SBATCH --mail-type FAIL
#SBATCH --cpus-per-task 12
#SBATCH -t 100:00:00
#SBATCH --mem 50G

source ReadConfig.sh $1

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" ${ORIDIR}/${SAMPLELIST})

module purge
module load java/1.9
module load miniconda2/4.4.10
source activate /gpfs/commons/groups/landau_lab/tprieto/conda/mgatk

# It is not working when there are two forward slash on the path name!!!
WORKDIR="${WORKDIR}${SAMPLE}/outs"

# I used the cellranger barcodes only for the ASAP sample (I think at the end I used it for everything??(22Jul))
# The tenth column within the singlecell.csv cell range output file contains 
# a binary indicator of whether barcode is associated with a cell
awk 'BEGIN { FS = "," }{print $1"\t"$10}' ${WORKDIR}/singlecell.csv  | \
	awk '{if ($2==1){print $1}}' > ${WORKDIR}/barcodes.${SAMPLE}.txt


mgatk tenx -i ${WORKDIR}/possorted_bam.bam \
	--mito-genome 'hg38' \
	-o ${WORKDIR}/mgatk \
	--ncores 12 \
  	-bt CB \
	-b ${WORKDIR}/barcodes.${SAMPLE}.txt \
	--alignment-quality 20 \
	--emit-base-qualities \
        --base-qual 20

source deactivate
