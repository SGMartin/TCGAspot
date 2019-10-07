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

#### Input folder structure
TCGASpot expects a certain input folder structure to be followed,  otherwise
it will fail. Below is an example based on the test suite, located in .test/samples:

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

### Dependencies

- Snakemake v.5.4 or higher.

A conda environment file (__tcgaspot.yaml__) can be found in __TCGASpot/envs/__. It contains a list
of the complete set of additional packages needed to run the pipeline end to end. 
You can either install it right away with:
```
conda env create -f /path/to/tcgaspot/envs/tcgaspot.yaml
conda activate tcgaspot
conda install snakemake
snakemake --jobs <ncores>
```
or you might run the pipeline letting snakemake handle it (**RECOMMENDED**):

```
snakemake  --use-conda --jobs  <ncores> --resources mem_mb=<max_ram_to_use>
```

### Example DAG
<p align="center">
  <img width="600" height="480" src="https://github.com/SGMartin/TCGAspot/blob/master/example_dag.svg">
</p>
