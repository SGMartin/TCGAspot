<h1>TCGASpot</h1>

[![Build Status](https://api.travis-ci.org/SGMartin/TCGAspot.svg?branch=master)](https://travis-ci.org/SGMartin/TCGAspot) [![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A55.4-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

A complete pipeline summarising my Master's degree thesis analysis

## Introduction
TCGASpot is a snakemake pipeline whose main purpose is to apply a pharmacogenomic tool called ***VulcanSpot*** to  large amounts of genomic data from ***The Cancer Genome Atlas*** (TCGA). Check VulcanSpot online at: http://vulcanspot.org/

To that end, it retrieves single point mutation, copy number variation and mRNA expression data and reconstruct the genomic landscape of all patients of input project. It then classifies every alteration in Gain of Function (GoF), Loss of Function (LoF) or Unknown effect (Unknown) at gene level.


### Input

- Single point mutation files: Only MuTecT2 workflow files are supported at the moment.
- Copy number variation files: Level 3, GISTIC gene level copy number score.
- RNAseq mRNA expression data: Raw HTSeq counts from xenabrowser
- Metadata from the whole download: needed to map alliquots to cases.
- (Optional): Annotations.txt files included with previous data

#### Input folder structure
TCGASpot expects a certain input folder structure to be followed,  otherwise
it will fail. Below is an example based on the test suite, located in _.test/samples_:

```bash
samples/
├── CNV
│   └── UVM
│       └── UVM.focal_score_by_genes.txt
├── MAF
│   └── UVM
│       └── TCGA.UVM.mutect.6c7b01bc-b068-4e01-8b4d-0362f5959f65.DR-10.0.somatic.maf
├── METADATA
│   └── metadata.json
└── MRNA
    └── UVM
        └── TCGA-UVM.htseq_counts.tsv
```


### Output

- A summary CSV table.
- Summary plots.

### Installation

| Python package | Version | Description |
| --- | --- | --- |
| DESeq2 | 1.26.0 | Bioconductor package for RNA-seq DE |
| doMC | 1.3.5 | R package for parallel processing 
| edgeR | 3.28.0 | Bioconductor package for RNA-seq DE |
| foreach | 1.4.4 | R package for parallel processing
| matplotlib | 3.1.1  | Python package: graphics library |
| mygene  | 3.1.0   | Python package: annotation library |
| numpy | 1.17.12 | Python package: vectoriced arrays library |
| pandas | 0.25.1 | Python package: dataframe library |
| python | 3.6.7 | Programming Language  |
| R | 3.6.1 | Programming Language |
| scipy | 1.3.1  | Python package: statistics library |
| seaborn | 0.9.0 | Python package: graphics library |
| snakemake | 5.4 | Python based workflow manager |

The above table is a glossary of all required packages. However, you do not have to manually install all of them by yourself! **Only Python 3 and snakemake need to be installed**. Remaining dependencies can be managed either by:

- (**RECOMMENDED**) Running the pipeline letting snakemake handle dependencies using a conda directive:

```
snakemake  --use-conda --jobs  <ncores> --resources mem_mb=<max_ram_to_use>
```

- (**Requires Conda**) Installing every dependency at once using conda and both environment files: (___tcgaspot.yaml___, ___deseq2.yaml___) which can be found in _TCGASpot/envs/_. They contains a list
of the complete set of additional packages needed to run the pipeline end to end. Once conda is installed, open a terminal
and type:
```
conda env create -f /path/to/tcgaspot/envs/tcgaspot.yaml
conda env update -f /path/to/tcgaspot/envs/deseq2.yaml
conda activate tcgaspot
conda install snakemake
snakemake --jobs <ncores> --resources mem_mb=<max_ram_to_use>
```




### Example DAG
<p align="center">
  <img width="600" height="480" src="https://gitlab.com/bu_cnio/TCGAspot/blob/master/example_dag.svg">
</p>
