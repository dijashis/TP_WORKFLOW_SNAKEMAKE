# ATACSEQ WORKFLOW FOR High Performance Computing Class

<h1>DIOP Khadidiatou</h1>

_Université Clermont auvergne HPC 2021_


This is a ATACSeq snakemake pipeline for High Performance Computing course. We first did the ATACSeq data analysis with slurm (https://github.com/dijashis/Projet_HPC) . We use ATAC-Seq data from *Gomez-Cabrero et al. (2019)* from a murine B3 cell line. 
One of the goals of the study is to identify new nuclear sites following translocation of the transcription factor Ikaros after exposure to the drug Tamoxifen The original data set has 50,000 cells collected per sample 3 replicates per sample and 2 cell stages: 0 and 24h (harvest time after drug treatment). 
Use (https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html) for reproducibility.


## To launch the pipeline you we need :

*  files that should be in _data/mydatalocal/atacseq/_ (_subsets_ and _bowtie2_ index)

*  One config files detailling all the different options and data needed to launch the pipe. _config.yaml_ in the config directory 

*  _Snakefiles_ (entrypoint of the workflow contains rules and scripts)

*  _env.yaml_ (need conda with bioconda and conda-forge and dependancies for all rules in snakefile)

**_dependencies_:

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

snakemake  --use-conda --cores all 

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
  Same steps as in (https://github.com/dijashis/Projet_HPC)
