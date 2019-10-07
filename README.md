<h1 align="center">TCGASpot</h1>
<img src=https://api.travis-ci.org/SGMartin/TCGAspot.svg?branch=master>
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

### Output

- A summary CSV table.
- Summary plots.

### Dependencies

- Snakemake v.5.4 or higher.

A conda environment file (tcgaspot.yaml) can be found in TCGASpot/envs/. It contains a list
of all additional packages needed to run the pipeline end to end. You can either install it right away with:
```
conda env create -f /path/to/tcgaspot/envs/tcgaspot.yaml
conda activate tcgaspot
conda install snakemake
snakemake --jobs <ncores>
```
or you might run the pipeline letting snakemake handle it:

```
snakemake  --use-conda --jobs  <ncores> 
```

### Example DAG
<p align="center">
  <img width="600" height="480" src="https://github.com/SGMartin/TCGAspot/blob/master/example_dag.svg">
</p>
