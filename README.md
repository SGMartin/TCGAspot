<h1 align="center">TCGASpot</h1>

<p align="center">
  <img width="600" height="480" src="https://github.com/SGMartin/TCGAspot/blob/master/example_dag.svg">
</p>

## Introduction

TCGASpot is a snakemake pipeline whose main purpose is to apply a pharmacogenomic tool called ***VulcanSpot*** to  large amounts of genomic data from ***The Cancer Genome Atlas*** (TCGA). Check VulcanSpot online at: http://vulcanspot.org/

### Input

- Single point mutation files: Only MuTecT2 workflow files are supported at the moment.
- Copy number variation files: Level 3, GISTIC gene level copy number score.
- RNAseq mRNA expression data: RNAseq-HTSeq-FPKM from xenabrowser

### Output

- A summary CSV table.
- Summary plots.

### Dependencies

- Snakemake
- Python's mygene package

