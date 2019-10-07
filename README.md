<h1>TCGASpot</h1>

[![Build Status](https://api.travis-ci.org/SGMartin/TCGAspot.svg?branch=master)](https://travis-ci.org/SGMartin/TCGAspot) [![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A55.4-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

A complete pipeline summarising my Master's degree thesis analysis

## Introduction
[WIP]
TCGASpot is a snakemake pipeline whose main purpose is to apply a pharmacogenomic tool called ***VulcanSpot*** to  large amounts of genomic data from ***The Cancer Genome Atlas*** (TCGA). Check VulcanSpot online at: http://vulcanspot.org/

To that end, it retrieves single point mutation, copy number variation and mRNA expression data and reconstruct the genomic landscape of all patients of input project. It then classifies every alteration in Gain of Function (GoF), Loss of Function (LoF) or Unknown effect (Unknown) at gene level.


### Input

- Single point mutation files: Only MuTecT2 workflow files are supported at the moment.
- Copy number variation files: Level 3, GISTIC gene level copy number score.
- RNAseq mRNA expression data: RNAseq-HTSeq-FPKM from xenabrowser
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
        └── TCGA-UVM.htseq_fpkm.tsv
```


### Output

- A summary CSV table.
- Summary plots.

### Installation
The following table is a glossary of all required packages. However, you do not
have to manually install all of them by yourself!

| Python package | Version | Description |
| --- | --- | --- |
| matplotlib | 3.1.1  | Python package: graphics library |
| mygene  | 3.1.0   | Python package: annotation library |
| numpy | 1.17.12 | Python package: vectoriced arrays library |
| pandas | 0.25.1 | Python package: dataframe library |
| python | 3.7.3 | Programming Language  |
| scipy | 1.3.1  | Python package: statistics library |
| seaborn | 0.9.0 | Python package: graphics library |
| snakemake | 5.4 | Python based workflow manager |

**Only Python 3 and snakemake need to be installed**. Remaining dependencies can be managed
either by:

- Running the pipeline letting snakemake handle dependencies using a conda directive (**RECOMMENDED**):

```
snakemake  --use-conda --jobs  <ncores> --resources mem_mb=<max_ram_to_use>
```

- (**Requires Conda**) Installing every dependency at once using conda and an environment file: (___tcgaspot.yaml___) which can be found in _TCGASpot/envs/_. It contains a list
of the complete set of additional packages needed to run the pipeline end to end. Once conda is installed, open a terminal
and type:
```
conda env create -f /path/to/tcgaspot/envs/tcgaspot.yaml
conda activate tcgaspot
conda install snakemake
snakemake --jobs <ncores> --resources mem_mb=<max_ram_to_use>
```




### Example DAG
<p align="center">
  <img width="600" height="480" src="https://github.com/SGMartin/TCGAspot/blob/master/example_dag.svg">
</p>
