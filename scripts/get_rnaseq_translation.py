#!/usr/bin/env python3
'''
RNAseq translation module: This scripts generates a new master table to translate
ENSEMBL id to Hugo Symbol. It shouldn't be run on an analysis basis and is 
included for completeness purposes.
'''
import mygene
import pandas as pd

#TODO: Implement an exhaustive filter to tackle multiple transcripts assigned
# to the same gene.

def main():
	"""
	Translates Ensembl IDs to HUGO using mygene queries.
	"""
	# SNAKEMAKE I/O #
	rnaseq_transcripts 				= snakemake.input[0]
	rnaseq_transcripts_translated	= snakemake.output[0]
	
	# Load data and trim version number from transcript ensembl id
	print("Annotating Ensembl IDs to HUGO...")

	ensembl_transcripts = pd.read_csv(rnaseq_transcripts, sep='\t')
	ensembl_transcripts['Ensembl_ID'] = ensembl_transcripts['Ensembl_ID'].str.split('.').str[0].str.strip()

	MyGene     = mygene.MyGeneInfo()
	gene_query = ensembl_transcripts['Ensembl_ID']

	annotated_names = MyGene.getgenes(gene_query,
									  fields='symbol',
									  as_dataframe=True
									 )
	
	# Get rid of unknown transcript-genes relationships
	annotated_names = annotated_names[~(annotated_names['notfound'] == True)]	

	# Drop duplicates... we are losing some info, but we can afford it. 7 genes
	# duplicated. Only one transcript will be used for these

	annotated_names = annotated_names.drop_duplicates(keep='first')
	annotated_names.drop(['_id', '_score', 'notfound'], axis=1, inplace=True)
	annotated_names.reset_index(inplace=True)
	annotated_names.columns = ['Ensembl_ID', 'Hugo_Symbol']

	annotated_names.to_csv(rnaseq_transcripts_translated, sep=',', index=False)

if __name__ == "__main__":
	main()

