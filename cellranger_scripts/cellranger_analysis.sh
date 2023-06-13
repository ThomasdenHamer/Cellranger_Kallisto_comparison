#!/bin/bash

#SBATCH --job-name=cellranger_analysis
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96000
#SBATCH --time=48:00:00


#FASTQS=/exports/sascstudent/thomas/inputdata/corrected_pbmc3k
FASTFOL=/exports/sascstudent/thomas/hiPSC_susana/
TRANSCRIPTOME=/exports/sascstudent/thomas/ref_indices/GRCh38_cellranger_2020A/refdata-gex-GRCh38-2020-A
OUTPUTFOL="/exports/sascstudent/thomas/output/cellranger/$1"

#Make local folders
LOCALTMP=/tmp/$USER/$SLURM_JOB_ID
srun mkdir -p $LOCALTMP
srun chmod 700 /tmp/$USER
srun mkdir -p $LOCALTMP/infasta
srun mkdir -p $LOCALTMP/inref
srun mkdir -p $LOCALTMP/output

#Copy input files
FASTAS=$(find $FASTFOL -iname "$2*.fastq.gz" | sort | tr "\n" " " | sed 's/ *$//g')
echo "Copying of the FASTQ input files"
echo "$FASTAS"
srun cp -fr --target-directory "$LOCALTMP/infasta/" $FASTAS
echo "Copying the input Transcriptome"
srun cp -fr $TRANSCRIPTOME/* $LOCALTMP/inref

#Cellranger run
echo "Start of the Cellranger run"
cd $LOCALTMP/output/
cmd="""srun /exports/sascstudent/thomas/cellranger-7.0.1/cellranger count --id=$2 \
		      --transcriptome=$LOCALTMP/inref \
		      --fastqs=$LOCALTMP/infasta \
		      --localcores=$SLURM_CPUS_PER_TASK \
		      --localmem=88 \
		      --include-introns=false"""
echo "Running $cmd"
eval $cmd
echo "Done Running"

#copy back results
echo "Copying back the output to $OUTPUTFOL"
srun rsync -r -q $LOCALTMP/output/ $OUTPUTFOL
echo "Copy back complete"

#Remove tmp folder
echo "Remove tmp folder"
rm -Rf $LOCALTMP/
echo "$SBATCH_JOB_NAME has completed"

