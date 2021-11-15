# ATACSEQ WORKFLOW FOR High Performance Computing Class

### DIOP Khadidiatou

### Universite Clermont Auvergne HPC 2021


This is a ATACSeq snakemake pipeline for High Performance Computing course. We first did the ATACSeq data analysis with slurm (https://github.com/dijashis/Projet_HPC) . We use ATAC-Seq data from *Gomez-Cabrero et al. (2019)* from a murine B3 cell line. 
One of the goals of the study is to identify new nuclear sites following translocation of the transcription factor Ikaros after exposure to the drug Tamoxifen The original data set has 50,000 cells collected per sample 3 replicates per sample and 2 cell stages: 0 and 24h (harvest time after drug treatment). 
Use (https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) for reproducibility.


## To launch the pipeline you we need :

*  Files that should be in `data/mydatalocal/atacseq/` (**subsets** and **bowtie2** index)

*  One config files detailling all the different options and data needed to launch the pipe. `config.yaml` in the config directory 

*  `Snakefiles` (entrypoint of the workflow contains rules and scripts)

*  `env.yaml` ( we need conda with bioconda and conda-forge and the following dependancies for all rules in snakefile )

**dependencies**:

 	- FastQC==0.11.9
	- Cutadapt==3.5                            
	- Bowtie2                                
	- samtools==1.14
	- r-base==4.1.1
	- openjdk==10.0.2
	- picard==2.26.5
	- deepTools 
	- MACS2
	- bedtools


## Launching the pipeline 

`snakemake  --use-conda` . Need to activate snakemake in conda before. 


## Directory ##
	├── .gitignore
	├── README.md
	├──Snakefile
    ├── env.yaml
	├── config
	│      ├── config.yaml
	│ 
	├── data
	│   └── mydatalocal
	│        └── atacseq 
    |            ├── bowtie2
	│                └──subset
	└── results
  
 ## WORKFLOW STEPS
  Same steps as in this repository (https://github.com/dijashis/Projet_HPC) 
