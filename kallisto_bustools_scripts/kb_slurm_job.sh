#!/bin/bash

#SBATCH --job-name=kbpython
#SBATCH --ntasks=1
#SBATCH --time=8:00:00
#SBATCH --mem=120000
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16


module purge
module add tools/miniconda/python3.8/4.9.2

#User variables
CONDAPATH=/exports/sascstudent/thomas/conda_envs/kbpython

INPUTFASTA=/exports/sascstudent/thomas/inputdata/corrected_pbmc3k
#INPUTFASTA=/exports/sascstudent/thomas/inputdata/mouse_brain1k/nuclei_900_fastqs

REFERENCEFOL=/exports/sascstudent/thomas/ref_genomes/GRCh37.75_homo_sapiens_ensemble
REFERENCEFASTA=Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz
REFERENCEGTF=Homo_sapiens.GRCh37.75.gtf.gz

#REFERENCEFOL=/exports/sascstudent/thomas/ref_genomes/mus_musculus_GRCm38_98
#REFERENCEFASTA=Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
#REFERENCEGTF=Mus_musculus.GRCm38.98.gtf.gz

OUTPUTFOL="/exports/sascstudent/thomas/output/bustools/$1"

TECHNOLOGY="10xv1"
WHITELIST=/exports/sascstudent/thomas/cellranger_v1_whitelist.txt

REFTYPE="intron"
COUNTWORKFLOW="nucleus"

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
srun mkdir -p $LOCALTMP/outputcount
srun mkdir -p $LOCALTMP/infasta
srun mkdir -p $LOCALTMP/inref


#Copy input files
echo "Copy input files"
echo "Copying the FASTA/Q files from $INPUTFASTA"
srun rsync -r -q $INPUTFASTA/* $LOCALTMP/infasta/
echo "Copying Reference fasta file from $REFERENCEFOL"
srun cp $REFERENCEFOL/$REFERENCEFASTA $LOCALTMP/inref/
echo "Copying Reference GTF file from $REFERENCEFOL"
srun cp $REFERENCEFOL/$REFERENCEGTF $LOCALTMP/inref/
if [[ $TECHNOLOGY == "10xv1" ]]
then
	echo "Copying 10XV1 Whitelist"
	cp $WHITELIST $LOCALTMP
	TECHNOLOGY=1,0,14:2,0,10:0,0,0
fi
echo "Done copying input files"


#Indexing of the reference genome
#nucl_ref_cmd="srun kb ref -i $LOCALTMP/outputref/transcriptome.idx -g $LOCALTMP/outputref/transcripts_to_genes.txt -f1 $LOCALTMP/outputref/cdna.fa -f2 $LOCALTMP/outputref/intron.fa -c1 $LOCALTMP/outputref/cdna_t2c.txt -c2 $LOCALTMP/outputref/intron_t2c --workflow lamanno --tmp $LOCALTMP/tempref  $LOCALTMP/inref/$REFERENCEFASTA $LOCALTMP/inref/$REFERENCEGTF"

if [[ $REFTYPE == "default" ]]
then
	echo "Generating a default Kallisto reference index"
	ref_cmd="srun kb ref -i $LOCALTMP/outputref/transcriptome.idx -g $LOCALTMP/outputref/tr2g.txt -f1 $LOCALTMP/outputref/cdna.fa --tmp $LOCALTMP/tempref $LOCALTMP/inref/$REFERENCEFASTA $LOCALTMP/inref/$REFERENCEGTF"
fi

if [[ $REFTYPE == "intron" ]]
then
	echo "Generating a intron splice aware Kallisto reference index"
	ref_cmd="srun kb ref -i $LOCALTMP/outputref/transcriptome.idx -g $LOCALTMP/outputref/tr2g.txt -f1 $LOCALTMP/outputref/cdna.fa --tmp $LOCALTMP/tempref -f2 $LOCALTMP/outputref/intron.fa -c1 $LOCALTMP/outputref/cdna_tr2c.txt -c2 $LOCALTMP/outputref/intron_tr2c.txt --workflow lamanno $LOCALTMP/inref/$REFERENCEFASTA $LOCALTMP/inref/$REFERENCEGTF"
fi
echo "Running kb reference command: $ref_cmd"
eval $ref_cmd
echo "Kallisto reference generation done"
echo ""

#FASTA files in order for 10xv1
#FASTAS=$(ls $LOCALTMP/infasta/* | sort -t '-' -k 3)
#FASTAS=$(ls $LOCALTMP/infasta/* | egrep "*R*.fastq.gz")
FASTAS=$(find $LOCALTMP/infasta/* -path '*R*.fastq.gz' | sort)
#echo $FASTAS
#Quantification of the dataset using the reference index

if [[ $TECHNOLOGY == *":"* ]]
then
	echo "Custom technology detected. Using provided whitelist"

	if [[ $COUNTWORKFLOW == "default" ]]
	then
		echo "Running default kb count workflow"
		count_cmd="""srun kb count -o $LOCALTMP/outputcount -i $LOCALTMP/outputref/transcriptome.idx -g $LOCALTMP/outputref/tr2g.txt -x $TECHNOLOGY -w $LOCALTMP/cellranger_v1_whitelist.txt --tmp $LOCALTMP/tempcount -t $SLURM_CPUS_PER_TASK -m 18000M --verbose $FASTAS"""
	fi

	if [[ $COUNTWORKFLOW == "lamanno" ]]
	then
		echo "Running a Lamanno count workflow"
		count_cmd="""srun kb count"""
	fi

	if [[ $COUNTWORKFLOW == "nucleus" ]]
	then
		echo "Running a Nucleus count workflow"
		count_cmd="""srun kb count -o $LOCALTMP/outputcount -i $LOCALTMP/outputref/transcriptome.idx -g $LOCALTMP/outputref/tr2g.txt -x $TECHNOLOGY -w $LOCALTMP/cellranger_whitelist_v1.txt --tmp $LOCALTMP/tempcount -t $SLURM_CPUS_PER_TASK -m 18000M --verbose -c1 $LOCALTMP/outputref/cdna_tr2c.txt -c2 $LOCALTMP/outputref/intron_tr2c.txt --workflow nucleus --h5ad $FASTAS"""
	fi

	if [[ $COUNTWORKFLOW == "kite" ]]
	then
		echo "Running a kite count workflow"
		count_cmd="""srun kb count"""
	fi
else
	echo "Using an included technology whitelist"

	if [[ $COUNTWORKFLOW == "default" ]]
	then
        	echo "Running default kb count workflow"
        	count_cmd="""srun kb count -o $LOCALTMP/outputcount -i $LOCALTMP/outputref/transcriptome.idx -g $LOCALTMP/outputref/tr2g.txt -x $TECHNOLOGY --tmp $LOCALTMP/tempcount -t $SLURM_CPUS_PER_TASK -m 18000M --verbose $FASTAS"""
	fi

	if [[ $COUNTWORKFLOW == "lamanno" ]]
	then
		echo "Running a Lamanno count workflow"
		count_cmd="""srun kb count"""
	fi

	if [[ $COUNTWORKFLOW == "nucleus" ]]
	then
		echo "Running a Nucleus count workflow"
		count_cmd="""srun kb count -o $LOCALTMP/outputcount -i $LOCALTMP/outputref/transcriptome.idx -g $LOCALTMP/outputref/tr2g.txt -x $TECHNOLOGY --tmp $LOCALTMP/tempcount -t $SLURM_CPUS_PER_TASK -m 18000M --verbose -c1 $LOCALTMP/outputref/cdna_tr2c.txt -c2 $LOCALTMP/outputref/intron_tr2c.txt --workflow nucleus --h5ad $FASTAS"""
	fi

	if [[ $COUNTWORKFLOW == "kite" ]]
	then
		echo "Running a kite count workflow"
		count_cmd="""srun kb count"""
	fi
fi

echo "Running kb count command: $count_cmd"
eval $count_cmd
echo "kb count has completed"
echo ""

#Copy back the results of count command
echo "Copying back the output data to $OUTPUTFOL"
srun rsync -r -q $LOCALTMP/outputref/ $OUTPUTFOL/
srun rsync -r -q $LOCALTMP/outputcount/ $OUTPUTFOL/
echo "copy back complete"


#Remove TMP folder
echo "Removing local temp folder"
rm -Rf $LOCALTMP/
echo "$SBATCH_JOB_NAME has finished"
