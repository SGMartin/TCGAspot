#!/usr/bin/env python3
"""
Section: Benchmarking VulcanSpot using TCGA data

Purpose: To assess VulcanSpot performance against a public dataset of cancer 
patients genomic data. 
"""

import json
import os
import sys

import pandas as pd


def main(rcnv, metadata, affy_db, save_as, metrics_save_as):

	### INPUT ###
	raw_cnv  = rcnv
	metadata = metadata
	affy_db  = affy_db

	###  OUTPUT ###
	where_to_save = save_as
	metrics 	  = metrics_save_as

	raw_cnv = pd.read_csv(raw_cnv, sep='\t')

	translated_cnv 		  = annotate_from_ensembl_to_hugo(raw_cnv, affy_db)
	filtered_cnv   		  = filter_cnv(translated_cnv)
	annotated_cnv  		  = match_cases_to_samples(filtered_cnv, metadata)
	annotated_cnv_samples = match_samples_to_tissue(annotated_cnv, metadata)

	# Check if annotations.txt is present
	# TODO: Log this
	input_cnv_directory = os.path.dirname(rcnv)
	annotation_file     = input_cnv_directory + '/annotations.txt'

	if(os.path.exists(annotation_file)):
		annotated_cnv= exclude_annotated_cases(annotated_cnv, annotation_file)


	report_metrics = generate_report_metrics(raw_cnv, annotated_cnv_samples)

	# SAVE RESULTS #
	annotated_cnv_samples.to_csv(where_to_save, sep=',', index=False)
	report_metrics.to_csv(metrics, sep=',', index=False)
	

def annotate_from_ensembl_to_hugo(cnv: pd.DataFrame, db: str) -> pd.DataFrame:
	'''
	Using a pregenerated, translated file, it matches SNP ENSEMBL values to
	MyGene HUGO values and gets rids of unmapped regions.
	'''

	# Load the database
	translated_affy_snp = pd.read_csv(db, sep=',', usecols=['Gene Symbol'])

	# these should match 1:1 because of how cnv data is arranged, so no need
	# for mapping
	cnv['Gene Symbol'] = translated_affy_snp['Gene Symbol']

	# delete unmapped regions
	cnv_alterations	= cnv[~cnv['Gene Symbol'].isnull()]

	# drop Gene ID and Cytoband  columns as these are not used
	cnv_alterations = cnv_alterations.drop(['Gene ID', 'Cytoband'], axis=1)

	del translated_affy_snp
	return cnv_alterations

def filter_cnv(cnv: pd.DataFrame) -> pd.DataFrame:
	'''
	Reshapes cnv files to a long format, using Gene symbol as id and then
	drop neutral cnv values.
	'''
	reshaped_cnv = pd.melt(frame=cnv,
	                       id_vars=['Gene Symbol'],
						   var_name='sample_id',
						   value_name='copy_number'
						  )

	# drop neutrals and rename
	filtered_cnv = reshaped_cnv[reshaped_cnv.loc[:,'copy_number'] != 0]
	gain_and_loss = {1: 'cnv_gain', -1: 'cnv_loss'}
	filtered_cnv = filtered_cnv.replace({'copy_number' : gain_and_loss})

	return filtered_cnv

#TODO: These two should be refactored
def match_samples_to_tissue(input_cnv: pd.DataFrame, metadata: str) -> pd.DataFrame:
	'''
	Translates sample ID from metadatada to TCGA barcode and then extracts
	original tissue using TCGA suffixes.
	https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
	'''

	with open(metadata) as metadata_file:
		raw_metadata = json.load(metadata_file)
	
	patient_barcodes = dict()

	for project in raw_metadata:
		for patients in project['associated_entities']:
			patient_barcodes[patients['entity_id']] = patients['entity_submitter_id']
	
	input_cnv['sample'] = input_cnv['sample_id'].map(patient_barcodes,
													 na_action='ignore'
													)
	input_cnv['sample'] = translate_barcode_to_sample_type(input_cnv['sample'])

	return input_cnv		

def match_cases_to_samples(input_cnv: pd.DataFrame,
						   metadata: str) -> pd.DataFrame:
	'''
	Translates samples ID to case UUID using metadata. Observations are
	named after alliquots. We have to rename them using the corresponding
	patient ID. This enables later matching with mutation data
	'''
	
	print("Reading metadata...")

	with open(metadata) as metadata_file:
		raw_metadata = json.load(metadata_file)
	
	patient_barcodes = dict()

	# Match cases and aliquots ids together. Case id matches maf files
	for project in raw_metadata:
		for patients in project['associated_entities']:
			patient_barcodes[patients['entity_id']] = patients['case_id']
	
	# map is faster than rename for large amounts of data, so we'll use that
	input_cnv['case_id'] = input_cnv['sample_id'].map(patient_barcodes,
													  na_action='ignore'
													  )	
	
	return input_cnv

def exclude_annotated_cases(rcnv: pd.DataFrame, annotations: str) -> pd.DataFrame:
	''' 
		Loads annotations.txt file for this project and returns a dataframe of
		excluded cases id. Cases/Items in annotations.txt are spared
		in certain cases based on notification category.
	'''
	cases_to_exclude = pd.read_csv(annotations,
								   sep='\t',
								   usecols=['entity_id', 'category']
								   )

	categories_to_keep = ['History of acceptable prior treatment related to a prior/other malignancy',
						  'Acceptable treatment for TCGA tumor',
						  'Alternate sample pipeline'
						 ]
	
	acceptable_alterations = cases_to_exclude['category'].isin(categories_to_keep)
	cases_to_exclude   	   = cases_to_exclude[~acceptable_alterations]

	excluded_cases   = rcnv['case_id'].isin(cases_to_exclude['entity_id'])
	filtered_cnv	 = rcnv[~excluded_cases]

	return filtered_cnv

def translate_barcode_to_sample_type(barcode:str) -> str:
	'''
	Translates tumor barcode to sample type according to
	https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
	'''
	CodesDictionary = {
        "01":"TP",
        "02":"TR",
        "03":"TB",
        "04":"TRBM",
        "05":"TAP",
        "06":"TM",
        "07":"TAM",
        "08":"THOC",
        "09":"TBM",
        "10":"NB",
        "11":"NT",
        "12":"NBC",
        "13":"NEBV",
        "14":"NBM",
        "15":"15SH",
        "16":"16SH",
        "20":"CELLC",
        "40":"TRB",
        "50":"CELL",
        "60":"XP",
        "61":"XCL",
        "99":"99SH"
    }

	sampleCode = str(barcode).split('-')
	doubleDigit = sampleCode[3][0:2]
	
	return CodesDictionary[doubleDigit]

def generate_report_metrics(rcnv: pd.DataFrame, acnv: pd.DataFrame) -> pd.DataFrame:
	'''
	Generates some metrics for each step performed. These are later saved
	to a file.
	'''
	#  Report metrics #
	raw_cases            = len(rcnv.columns) 
	raw_cases 			-= 3 # accounting for the gene id, cytoband and ensembl cols.
	raw_alterations 	 = len(rcnv.index)
	
	filtered_cases		 = acnv['case_id'].nunique()
	filtered_alterations = len(acnv.index)

	report = {'Item':['cases',
					  'cases_filtered',
					  'alterations',
					  'alterations_filtered'
					  ],
			  'Count':[raw_cases,
			  		   filtered_cases,
					   raw_alterations,
					   filtered_alterations
					   ]
			 }
	
	report = pd.DataFrame(report)

	return report

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])