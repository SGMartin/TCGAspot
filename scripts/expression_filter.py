#!/usr/bin/env python3
'''
mRNA expression filtering module: This script takes a matrix of RNA-seq mRNA counts
from XenaBrowser and prepares it for DESEq2 analysis. It reverses xena transformation
and filters out probes with low counts. Samples are mapped to cases ID and 
normals are removed as these are not used.

Additionally, only probes whose HUGO_SYMBOL is present in affy snp 6.0 array are
kept. The reason behind this is simple: only cnv_gains events are checked using
this data.
'''

#TODO: To me, for the future: Better naming for method vars.

import json
import pandas as pd
import numpy  as np

def main():

	pd.set_option('mode.chained_assignment', 'raise')

	# ------- SNAKEMAKE I/O ------ #

	raw_xena_data     = snakemake.input[0]
	rnaseq_annotation = snakemake.input[1]
	metadata		  = snakemake.input[2]

	where_to_save	  = snakemake.output[0]

	# Load dataframe	
	raw_xena_data = load_dataframe(raw_xena_data)

	# Annotate samples and remove those not matching to any case
	annotated_xena = annotate_samples_to_cases(raw_xena_data, metadata)

	collapsed_raw_xena = collapse_probes_to_cnv(raw_xena     = annotated_xena,
												rnaseq_annot = rnaseq_annotation
												)

	raw_count_xena 	   = raw_counts_from_xena(logtransformed=collapsed_raw_xena)

	# Drop duplicated probes, there are only 4-6 of over 19000 and their gscore
	# low anyway. TODO: Maybe a more thoughtful collapsing method
	raw_count_xena.drop_duplicates(subset='Hugo_Symbol', keep='first', inplace=True)
	raw_count_xena.to_csv(where_to_save, sep=',', index=False)


# TODO: sample filtering based on metadata or annotations could be done here
def load_dataframe(raw_xena_data:str) -> pd.DataFrame:
	'''
	This method loads RNA-seq expr. data. To make it efficient, only
	tumoral samples are loaded. Returns a pandas dataframe.
	'''

	# Get all columns and use them for usecols later
	file_header = ""
	with(open(raw_xena_data, 'r')) as buffer:
		file_header = buffer.readline()
	
	#get rid of trailing new line and get an array of tcga using \t
	file_header = file_header.rstrip()
	file_header = file_header.split('\t')

	samples = []
	for psample in file_header:

		if psample == 'Ensembl_ID':
			samples.append(psample)
		else:
			if int(psample.split('-')[3][0:2]) < 10: # normal samples are 11-12-13-14 + A/B/C
				samples.append(psample)

	data = pd.read_csv(raw_xena_data, sep='\t', usecols=samples)
	return data


def annotate_samples_to_cases(xena_data: pd.DataFrame, metadata:str) -> pd.DataFrame:
	'''
	Translates samples ID to case UUID using metadata. Returns a dataframe with
	all relevant columns renamed to cases.
	'''
	with open(metadata) as metadata_file:
		raw_metadata = json.load(metadata_file)
	
	patient_barcodes = dict()

	for project in raw_metadata:
		for patients in project['associated_entities']:

			trimmed_barcode = patients['entity_submitter_id'].split('-')[0:4]
			trimmed_barcode = '-'.join(trimmed_barcode)
			patient_barcodes[trimmed_barcode] = patients['case_id']
	
	# Unmapped cases and Ensembl_ID/Hugo_Symbol column remains unchanged
	xena_data.rename(columns=patient_barcodes, inplace=True)	
	
	# Filter samples not matching to any patient
	xena_data.drop(list(xena_data.filter(regex='TCGA')), axis=1, inplace=True)

	#TODO: improve this. For now use the first sample as case representative
	xena_data_annotated = xena_data.loc[:,~xena_data.columns.duplicated()].copy()
	
	return xena_data_annotated

# This method is not collapsing anymore, just annotating.
# TODO: Rename it
def collapse_probes_to_cnv(raw_xena: pd.DataFrame,
						   rnaseq_annot:str) -> pd.DataFrame:
	'''
	This method annotates expression transcripts and keeps those whose
	Hugo_Symbol matches affymetrix snp  6.0 array data including duplicates. 
	Returns a dataframe with the matched rows.
	'''

	# load rnaseq annotation
	rnaseq_annot = pd.read_csv(rnaseq_annot, sep=',', index_col='Ensembl_ID')
	rnaseq_annot = rnaseq_annot['Hugo_Symbol'].to_dict()

	#get rid of the version id
	raw_xena['Ensembl_ID']  = raw_xena['Ensembl_ID'].str.split('.').str[0].str.strip()
	raw_xena['Hugo_Symbol'] = raw_xena['Ensembl_ID'].map(rnaseq_annot)

	# attempt to free memory asap
	raw_xena.drop('Ensembl_ID', axis=1, inplace=True)
	del rnaseq_annot

	# Drop probes with no hugo gene associated
	raw_xena.dropna(axis='index', how='any', inplace=True)

	filtered_xena = raw_xena
	
	return filtered_xena	

def raw_counts_from_xena(logtransformed: pd.DataFrame) -> pd.DataFrame:
	'''
	This method attempts to reverse xena log2 tranform. to return raw counts
	for further analysis. Probes with <10 counts across all samples are removed 
	to speed up downstream analysis. Returns a dataframe of raw counts.
	'''

	# Xena performs a log2 +1 transformation to raw counts
	reverse_xena = np.exp2(logtransformed.loc[:, logtransformed.columns != 'Hugo_Symbol'])
	reverse_xena = reverse_xena - 1 
	reverse_xena = reverse_xena.round(decimals=0)
	reverse_xena['Hugo_Symbol'] = logtransformed['Hugo_Symbol']
	
	# Remove probes with very low counts (<10) across all samples.
#	reverse_xena['keep'] = reverse_xena.loc[:, reverse_xena.columns != 'Hugo_Symbol'].sum(axis='columns')
#	reverse_xena = reverse_xena[reverse_xena['keep'] >= 10]
#	reverse_xena.drop('keep', axis=1, inplace=True)

	return reverse_xena


if __name__ == "__main__":
	main()
