#!/bin/bash

#SBATCH --job-name=kbpython
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#SBATCH --mem=120000
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8


module purge
module add tools/miniconda/python3.8/4.9.2

#User variables
CONDAPATH=/exports/sascstudent/thomas/conda_envs/kbpython

#REFERENCEFOL=/exports/sascstudent/thomas/ref_genomes/GRCh37.75_homo_sapiens_ensemble
#REFERENCEFASTA=Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
#REFERENCEGTF=Homo_sapiens.GRCh37.75.gtf.gz

#REFERENCEFOL=/exports/sascstudent/thomas/ref_genomes/mus_musculus_GRCm38_98
#REFERENCEFASTA=Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
#REFERENCEGTF=Mus_musculus.GRCm38.98.gtf.gz

REFERENCEFOL=/exports/sascstudent/thomas/ref_indices/GRCh38_cellranger_2020A/refdata-gex-GRCh38-2020-A
REFERENCEFASTA=genome.fa.gz
#REFERENCEGTF=genes.gtf.gz
REFERENCEGTF=paralogous_intron.gtf

OUTPUTFOL="/exports/sascstudent/thomas/ref_indices/$1"

REFTYPE="intron"

#Check user var

if [ -z "$1" ]
then
	echo "No Ouput folder given"
	exit 1
fi

#Conda activate
conda activate $CONDAPATH


#Setup the tmp folders on the cluster
LOCALTMP=/tmp/$USER/$SLURM_JOB_ID
srun mkdir -p $LOCALTMP
srun chmod 700 /tmp/$USER
srun mkdir -p $LOCALTMP/outputref
srun mkdir -p $LOCALTMP/inref


#Copy input files
echo "Copy input files"
echo "Copying Reference fasta file from $REFERENCEFOL"
srun cp $REFERENCEFOL/fasta/$REFERENCEFASTA $LOCALTMP/inref/
echo "Copying Reference GTF file from $REFERENCEFOL"
#srun cp $REFERENCEFOL/$REFERENCEGTF $LOCALTMP/inref/
srun cp /home/tjjdenhamer/$REFERENCEGTF $LOCALTMP/inref/
echo "Done copying input files"

echo


#Indexing of the reference genome
if [[ $REFTYPE == "default" ]]
then
	echo "Generating a default Kallisto reference index"
	ref_cmd="srun kb ref -i $LOCALTMP/outputref/transcriptome.idx -g $LOCALTMP/outputref/tr2g.txt -f1 $LOCALTMP/outputref/cdna.fa --tmp $LOCALTMP/tempref $LOCALTMP/inref/$REFERENCEFASTA $LOCALTMP/inref/$REFERENCEGTF"
fi

if [[ $REFTYPE == "intron" ]]
then
	echo "Generating an intron splice aware Kallisto reference index"
	ref_cmd="srun kb ref -i $LOCALTMP/outputref/transcriptome.idx -g $LOCALTMP/outputref/tr2g.txt -f1 $LOCALTMP/outputref/cdna.fa --tmp $LOCALTMP/tempref -f2 $LOCALTMP/outputref/intron.fa -c1 $LOCALTMP/outputref/cdna_tr2c.txt -c2 $LOCALTMP/outputref/intron_tr2c.txt --workflow lamanno $LOCALTMP/inref/$REFERENCEFASTA $LOCALTMP/inref/$REFERENCEGTF"
fi
echo "Running kb reference command: $ref_cmd"
eval $ref_cmd
echo "Kallisto reference generation done"
echo ""

#Copy back the results of count command
echo "Copying back the output data to $OUTPUTFOL"
srun rsync -r -q $LOCALTMP/outputref/ $OUTPUTFOL/
echo "copy back complete"


#Remove TMP folder
echo "Removing local temp folder"
rm -Rf $LOCALTMP/
echo "$SBATCH_JOB_NAME has finished"
