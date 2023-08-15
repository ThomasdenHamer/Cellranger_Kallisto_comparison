# Code for the comparison of Cellranger and Kallisto

This repository contains all scripts used in the comparison of Cellranger and Kallisto | Bustools done in the following paper, TBD

## Folders content

* The cellranger_scripts folder has one SLURM script that runs the default cellranger analysis and uses a pre-built STAR index
* The kallisto_bustools_scripts folder has three scripts. The kb_slurm_job.sh is not depreciated. The kb_ref_job.sh SLURM script generates a kallisto index from a genome and GTF file using the KallistoBUStools wrapper and can create both a default and intronic index. The kb_count_job.sh SLURM script uses the KallistoBUStools wrapper to run a user specified workflow using the given technology. The different workflows are described in the KallistoBUStools wrapper documentation found here: [https://www.kallistobus.tools/](https://www.kallistobus.tools/)
* The starsolo_scripts folder has two scripts. The star_index_slurm.sh SLURM script creates a updated STAR index of the Cellranger pre-built STAR index using the genome and GTF file from this pre-built index. The starsolo_slurm.sh SLURM script runs the STARsolo program with user-specified parameters.
* The r_analysis script contains all R script used to analyse this specific project. The different scripts use a variety of different packages but the Seurat package is neccesary for all the single cell analyses.


## The Usage of various BASH scripts made for the analysis of the corresponding paper.

### star_index_slurm.sh

The star_index_slurm.sh SLURM script creates a STAR index. This script was used to generate a STAR index to be used in STARsolo.
The pre-built index provided by Cellranger did not work because it was built using an older version of STAR then the STAR version used for STARsolo.

The uses STAR version 2.7.10b

All the variables here are hardcoded into the scripts since this isn't meant to be run more than once

### starsolo_slurm.sh

The starsolo_slurm.sh SLURM script performs a STARsolo mapping/alignment, counting and filtering on a given dataset.

Usage


> ```sbatch starsolo.sh <output_folder> <index_folder> <FASTQ file[s]>```


All runs to process must be made up of one R1 and one R2 FASTQ file. The R1 contains the Barcode and UMI. The R2 file contains the cDNA fragment.
This is consistent with the V2/V3 10xChromium chemistry.

The input FASTQ[s] do not have to be given in any specific order as the scripts automatically groups R1 and R2 reads together with their lane information.
Because of this the files do need a consistent naming convention. This means only the lane and the R of the file can change

Examples




