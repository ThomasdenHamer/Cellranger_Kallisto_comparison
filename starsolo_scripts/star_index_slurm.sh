#!/bin/bash

#SBATCH --job-name=STARindex
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=48000
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8


#srun mkdir /tmp/$USER/$SLURM_JOB_ID
#srun chmod 700 /tmp/$USER
#cd /tmp/$USER/$SLURM_JOB_ID

echo "Generating STAR index"

cmd="""srun /exports/sascstudent/thomas/STAR/STAR-2.7.10b/source/STAR --runThreadN 8 --runMode genomeGenerate \
	--genomeDir /exports/sascstudent/thomas/ref_indices/GRCh38_STAR/ \
	--genomeFastaFiles /exports/sascstudent/thomas/ref_indices/GRCh38_cellranger_2020A/refdata-gex-GRCh38-2020-A/fasta/genome.fa \
	--sjdbGTFfile /exports/sascstudent/thomas/ref_indices/GRCh38_cellranger_2020A/refdata-gex-GRCh38-2020-A/genes/genes.gtf
"""
echo "Running $cmd"
eval $cmd
echo "Index complete. Ending script"

