#!/bin/bash

#SBATCH --job-name=STARsolo_count
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=8:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64000

module purge

#Input data
STARINDEXFOL=/exports/sascstudent/thomas/ref_indices/GRCh38_STAR

#Output folder
RESULTSFOLDER=/exports/sascstudent/thomas/output/STARsolo/$1

#User Variables
SOLOTYPE=CB_UMI_Simple
CBWHITELIST=/exports/sascstudent/thomas/cellranger_v1_whitelist.txt
CBWHITELISTNAME="cellranger_v1_whitelist.txt"
CBSTART=1
CBLEN=14
UMISTART=15
UMILEN=10
CELLFILTER=EmptyDrops_CR
SOLOFEATURES=Gene
SOLOMULTIMAPPING=Unique
SOLOUMIDEDUP=1MM_All


#Setup local tmp folders on cluster node
LOCALTMP=/tmp/$USER/$SLURM_JOB_ID
srun mkdir -p $LOCALTMP
srun chmod 700 /tmp/$USER
srun mkdir $LOCALTMP/fastinput
srun mkdir $LOCALTMP/refinput
srun mkdir $LOCALTMP/output


#Copy input files
echo "Copying cmdline input files"
#Takes all positional arguments except the first one to use for copying the input files.
#The order of the input files does not matter, they only need their matching pair.
FILES="${@:2}"
for (( i=0; i<${#FILES[@]}; i++ ))
do
	echo "Copying: ${FILES[$i]}"
	srun cp -fr ${FILES[$i]} "$LOCALTMP/fastinput"
done


#Split given input files into R1 and R2. R1 containing the BC+UMI. R2 containing the cDNA fragment. This convention is consistent with the V3 chemistry of 10xChromium
R1reads=$(find $LOCALTMP/fastinput/ -iname "*_R1_*.gz" | sort | tr "\n" "," | sed 's/,*$//g')
R2reads=$(find $LOCALTMP/fastinput/ -iname "*_R2_*.gz" | sort | tr "\n" "," | sed 's/,*$//g')
echo "R1 files: $R1reads"
echo "R2 files: $R2reads"


#Copy reference. Default human reference built using the FASTA and GTF files from the Cellranger STAR index.
echo "copying reference from $STARINDEXFOL"
srun cp $STARINDEXFOL/* $LOCALTMP/refinput/

#Copy whitelist
echo "Copying whitelist: $CBWHITELIST"
srun cp $CBWHITELIST $LOCALTMP/
echo "Done copying all files"


#Running STARsolo command
solocmd="""srun /exports/sascstudent/thomas/STAR/STAR-2.7.10b/source/STAR --genomeDir $LOCALTMP/refinput/ \
	--soloType $SOLOTYPE --soloCBwhitelist $LOCALTMP/$CBWHITELISTNAME \
	--soloCBstart $CBSTART --soloCBlen $CBLEN --soloUMIstart $UMISTART --soloUMIlen $UMILEN \
	--soloFeatures $SOLOFEATURES --soloUMIdedup $SOLOUMIDEDUP --soloMultiMappers $SOLOMULTIMAPPING --soloCellFilter $CELLFILTER \
	--readFilesCommand zcat --runThreadN $SLURM_CPUS_PER_TASK \
	--readFilesIn $R1reads $R2reads"""


echo "Running STARsolo command: $solocmd"
cd $LOCALTMP/output/
eval $solocmd
echo "STARsolo complete"

#Copy back results to user provided resultsfolder. Creates folder and parents if needed.
echo "Copying results to $RESULTSFOLDER"
srun mkdir -p $RESULTSFOLDER
srun cp -fr $LOCALTMP/output/* $RESULTSFOLDER
echo "Copy back done"

#Remove local tmp folder
echo "Removing local temp folder"
srun rm -Rf $LOCALTMP/
echo "Script complete"

