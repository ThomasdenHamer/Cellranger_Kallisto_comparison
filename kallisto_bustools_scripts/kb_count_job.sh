#!/bin/bash

#SBATCH --job-name=kbcount
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem=100000
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12


module purge
module add tools/miniconda/python3.8/4.9.2

#User variables
CONDAPATH=/exports/sascstudent/thomas/conda_envs/kbpython

#INPUTFASTA=/exports/sascstudent/thomas/inputdata/corrected_pbmc3k
#INPUTFASTA=/exports/sascstudent/thomas/inputdata/mesonephros_b_10x
#INPUTFASTA=/exports/sascstudent/thomas/inputdata/mouse_brain1k/nuclei_900_fastqs
#INPUTFASTA=/exports/sascstudent/thomas/hiPSC_susana
#INPUTFASTA=/exports/sascstudent/thomas/inputdata/chen_public_dataset


INDEXFOLDER=/exports/sascstudent/thomas/ref_indices/GRch38_kallisto_cdna
#INDEXFOLDER=/exports/sascstudent/thomas/ref_indices/GRch38_kallisto_intron
#INDEXFOLDER=/exports/sascstudent/thomas/ref_indices/GRcm38_kallisto_intron
#INDEXFOLDER=/exports/sascstudent/thomas/ref_indices/small_intronic

OUTPUTFOL="/exports/sascstudent/thomas/output/bustools/$1"

TECHNOLOGY="10xv3"
WHITELIST=/exports/sascstudent/thomas/cellranger_v1_whitelist.txt

COUNTWORKFLOW="default"

#Check user var
if [ -z "$1" ]
then
	echo "No Ouput folder given"
	exit 1
fi

if [ -z "$2" ] && [ -z "$3" ] && [ -z "$4" ] && [ -z "$5" ]
then
	echo "No input files given"
	exit 1
fi

#if [ -z "$2" ] && [ -z "$3" ]
#then
#	echo "No input files given"
#	exit 1
#fi

#Conda activate
conda activate $CONDAPATH


#Setup the tmp folders on the cluster
LOCALTMP=/tmp/$USER/$SLURM_JOB_ID
srun mkdir -p $LOCALTMP
srun chmod 700 /tmp/$USER
srun mkdir -p $LOCALTMP/outputcount
srun mkdir -p $LOCALTMP/infasta
srun mkdir -p $LOCALTMP/outputref


#Copy input files
echo "Copy input files"
echo "Copying the FASTA/Q files from $INPUTFASTA"

#srun cp $INPUTFASTA/*_R*.fastq.gz $LOCALTMP/infasta
#srun cp $INPUTFASTA/SRR10587$2/*_[12].fastq.gz $LOCALTMP/infasta/
#srun cp $INPUTFASTA/SRR10587$3/*_[12].fastq.gz $LOCALTMP/infasta/
srun cp $2 $LOCALTMP/infasta/
srun cp $3 $LOCALTMP/infasta/
srun cp $4 $LOCALTMP/infasta/
srun cp $5 $LOCALTMP/infasta/

#srun cp $2 $LOCALTMP/infasta/
#srun cp $3 $LOCALTMP/infasta/

echo "Copying index files from $INDEXFOLDER"
srun cp $INDEXFOLDER/*.txt $LOCALTMP/outputref
if [[ $TECHNOLOGY == "10xv1" ]]
then
	echo "Copying 10XV1 Whitelist"
	cp $WHITELIST $LOCALTMP/
	TECHNOLOGY=1,0,14:2,0,10:0,0,0
fi
echo "Done copying input files"


#FASTA files in order for 10xv1
#FASTAS=$(find $LOCALTMP/infasta/* -path '*.fastq.gz' | sort)
FASTAS=$(find $LOCALTMP/infasta/* -path '*R*.fastq.gz' | sort)
#echo $FASTAS
#Quantification of the dataset using the reference index

if [[ $TECHNOLOGY == *":"* ]]
then
	echo "Custom technology detected. Using provided whitelist"

	if [[ $COUNTWORKFLOW == "default" ]]
	then
		echo "Running default kb count workflow"
		count_cmd="""srun kb count -o $LOCALTMP/outputcount -i $INDEXFOLDER/transcriptome.idx -g $LOCALTMP/outputref/tr2g.txt -x $TECHNOLOGY -w $LOCALTMP/cellranger_v1_whitelist.txt --tmp $LOCALTMP/tempcount -t $SLURM_CPUS_PER_TASK -m 85G --verbose $FASTAS"""
	fi

	if [[ $COUNTWORKFLOW == "lamanno" ]]
	then
		echo "Running a Lamanno count workflow"
		count_cmd="""srun kb count"""
	fi

	if [[ $COUNTWORKFLOW == "nucleus" ]]
	then
		echo "Running a Nucleus count workflow"
		count_cmd="""srun kb count -o $LOCALTMP/outputcount -i $INDEXFOLDER/transcriptome.idx -g $LOCALTMP/outputref/tr2g.txt -x $TECHNOLOGY -w $LOCALTMP/cellranger_v1_whitelist.txt --tmp $LOCALTMP/tempcount -t $SLURM_CPUS_PER_TASK -m 85G --verbose -c1 $LOCALTMP/outputref/cdna_tr2c.txt -c2 $LOCALTMP/outputref/intron_tr2c.txt --workflow nucleus --h5ad $FASTAS"""
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
        	count_cmd="""srun kb count -o $LOCALTMP/outputcount -i $INDEXFOLDER/transcriptome.idx -g $LOCALTMP/outputref/tr2g.txt -x $TECHNOLOGY --tmp $LOCALTMP/tempcount -t $SLURM_CPUS_PER_TASK -m 85G --verbose $FASTAS"""
	fi

	if [[ $COUNTWORKFLOW == "lamanno" ]]
	then
		echo "Running a Lamanno count workflow"
		count_cmd="""srun kb count"""
	fi

	if [[ $COUNTWORKFLOW == "nucleus" ]]
	then
		echo "Running a Nucleus count workflow"
		count_cmd="""srun kb count -o $LOCALTMP/outputcount -i $INDEXFOLDER/transcriptome.idx -g $LOCALTMP/outputref/tr2g.txt -x $TECHNOLOGY --tmp $LOCALTMP/tempcount -t $SLURM_CPUS_PER_TASK -m 85G --verbose -c1 $LOCALTMP/outputref/cdna_tr2c.txt -c2 $LOCALTMP/outputref/intron_tr2c.txt --workflow nucleus --h5ad $FASTAS"""
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
srun rsync -r -q $LOCALTMP/outputref/*.txt $OUTPUTFOL/
srun rsync -r -q $LOCALTMP/outputcount/ $OUTPUTFOL/
echo "copy back complete"


#Remove TMP folder
echo "Removing local temp folder"
rm -Rf $LOCALTMP/
echo "$SBATCH_JOB_NAME has finished"
