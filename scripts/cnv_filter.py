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

	# Hugo and metadata Annot.
	ensembl_cnv     = annotate_from_ensembl_to_hugo(raw_cnv, affy_db)
	annotated_cnv   = annotate_cnv_from_metadata(ensembl_cnv, metadata)

	# Check if annotations.txt is present
	# TODO: Log this
	input_cnv_directory = os.path.dirname(rcnv)
	annotation_file     = input_cnv_directory + '/annotations.txt'

	if(os.path.exists(annotation_file)):
		annotated_cnv = exclude_annotated_items(annotated_cnv, annotation_file)

	# filter CNV neutrals and multiple aliquots
	filtered_cnv = filter_cnv(annotated_cnv)					

	# SAVE RESULTS #
	filtered_cnv.to_csv(where_to_save, sep=',', index=False)

	# Get metrics #
	report_metrics = generate_report_metrics(raw_cnv, filtered_cnv)
	report_metrics.to_csv(metrics, sep=',', index=False)
	

def annotate_from_ensembl_to_hugo(cnv: pd.DataFrame, db: str) -> pd.DataFrame:
	'''
	Using a pregenerated, translated file, it matches SNP ENSEMBL values to
	MyGene HUGO values and gets rids of unmapped and ambiguous regions.
	Returns a melted, reshaped dataframe
	'''

	# Load the database
	translated_affy_snp = pd.read_csv(db, sep=',', usecols=['Gene Symbol'])

	# these should match 1:1 because of how cnv data is arranged, so no need
	# for manual index mapping
	cnv['Gene Symbol'] = translated_affy_snp['Gene Symbol']

	# delete unmapped genomic regions
	cnv_alterations	= cnv[~cnv['Gene Symbol'].isnull()]

	# drop Gene ID and Cytoband  columns as these are not used
	cnv_alterations = cnv_alterations.drop(['Gene ID', 'Cytoband'], axis=1)

	# reshape data for easier manip.
	cnv_alterations = pd.melt(frame=cnv_alterations,
							  id_vars='Gene Symbol',
							  var_name='sample_id',
							  value_name='copy_number'
							 )
	
	# get genes which mapped to two regions
	is_ambiguous = translated_affy_snp['Gene Symbol'].duplicated(keep='first')
	ambiguous_genes = translated_affy_snp[is_ambiguous]

	# drop reference, unmapped regions
	ambiguous_genes = ambiguous_genes[~ambiguous_genes['Gene Symbol'].isnull()]
	ambiguous_genes = ambiguous_genes['Gene Symbol']

	# get ambiguous genes with shared alterations. These will be kept.
	is_cnv_ambiguous = cnv_alterations['Gene Symbol'].isin(ambiguous_genes)
	cnv_ambiguous    = cnv_alterations[is_cnv_ambiguous]
	
	is_cnv_ambiguous_consensus = cnv_ambiguous.duplicated(keep='first')
	cnv_ambiguous_consensus    = cnv_ambiguous[is_cnv_ambiguous_consensus]

	# drop ambiguous, non consensual from data
	cnv_alterations = cnv_alterations[~is_cnv_ambiguous]

	# rescue consensus data
	cnv_alterations = cnv_alterations.append(cnv_ambiguous_consensus)

	del translated_affy_snp
	del ambiguous_genes
	del cnv_ambiguous_consensus
	return cnv_alterations


def annotate_cnv_from_metadata(cnv: pd.DataFrame, metadata:str) -> pd.DataFrame:
	'''
	Annotates CNV using metadata. Returns a dataframe with aliquot and sample 
	columns.
	'''
	# get sample barcode and alliquot
	cnv_annotated 			 = get_barcode_from_sample(cnv, metadata)
	cnv_annotated['aliquot'] = cnv_annotated['barcode'].str.split('-').str[3].str[2]

	# get sample type from barcode
	cnv_annotated['sample']	= translate_barcode_to_tumor(cnv_annotated['barcode'])
	cnv_annotated.drop('barcode', axis=1, inplace=True)

	# get case id from sample
	cnv_annotated_samples    = 	get_case_from_sample(cnv_annotated, metadata)

	return cnv_annotated_samples

def exclude_annotated_items(cnv: pd.DataFrame, annotations:str) -> pd.DataFrame:
	'''
	Loads annotations.txt and filters out flagged items and cases. Cases/Items 
	are spared in certain cases based on notification category. Returns a
	filtered dataframe
	'''
	# Exclude both annotated aliquots and cases
	flagged_items = pd.read_csv(annotations,
								sep='\t',
								usecols=['entity_type', 
								   		 'entity_id',
										 'category'
										]
								)

	spared_categories  = ['History of acceptable prior treatment related to a prior/other malignancy',
						  'Acceptable treatment for TCGA tumor',
						  'Alternate sample pipeline'
						 ]
	
	flagged_items = flagged_items[~flagged_items['category'].isin(spared_categories)]

	# get cases and filter them
	excluded_cases = flagged_items[flagged_items['entity_type'] == 'case']['entity_id']
	cnv			   = cnv[~cnv['case_id'].isin(excluded_cases)]

	# get aliquots and filter them
	aliquots = flagged_items[flagged_items['entity_type'] == 'aliquot']['entity_id']
	cnv_clean = cnv[~cnv['sample_id'].isin(aliquots)]

	del flagged_items
	return cnv_clean

def filter_cnv(cnv: pd.DataFrame) -> pd.DataFrame:
	'''
	Drop entries when multiple aliquots from the same source do not match. Then
	it drops neutrals and rename gains and loss to cnv_gain and cnv_loss. 
	Returns a filtered dataframe.
	'''
	#TODO: maybe count aliquots and get a consensus

	# Get cases alterations from with multiple aliquots from the same sample
	mult_aliquots = cnv.duplicated(subset=['sample', 'Gene Symbol', 'case_id'])

	# Keep cases with multiple aliquots and drop all but consensus, of which
	# one is kept
	consensus_aliquots = cnv[mult_aliquots]
	consensus_aliquots = consensus_aliquots[consensus_aliquots.duplicated(keep='first')]

	# now drop all multiple aliquots instances from the mother matrix
	cnv = cnv[~mult_aliquots]
	cnv = cnv.append(consensus_aliquots) # rescue consensus

	# Now drop neutrals and rename
	#annotated_cnv = reshaped_cnv[reshaped_cnv.loc[:,'copy_number'] != 0]

	filtered_cnv = cnv[cnv.loc[:,'copy_number'] != 0]
	
	gof_lof = {-1:'cnv_loss', 1:'cnv_gain'}
	clean_cnv = filtered_cnv.replace({'copy_number': gof_lof})

	del filtered_cnv
	del consensus_aliquots
	return clean_cnv


def get_barcode_from_sample(input_cnv: pd.DataFrame, metadata: str) -> pd.DataFrame:
	'''
	Translates sample ID from metadatada to TCGA barcode and then extracts
	original tissue and aliquot number.
	'''

	with open(metadata) as metadata_file:
		raw_metadata = json.load(metadata_file)
	
	patient_barcodes = dict()

	for project in raw_metadata:
		for patients in project['associated_entities']:
			patient_barcodes[patients['entity_id']] = patients['entity_submitter_id']
	
	input_cnv['barcode'] = input_cnv['sample_id'].map(patient_barcodes,
													 na_action='ignore'
													)
	return input_cnv		


def get_case_from_sample(input_cnv: pd.DataFrame, metadata: str) -> pd.DataFrame:
	'''
	Translates samples ID to case UUID using metadata. Observations are
	named after alliquots. Returns a dataframe with a new column called
	case_id
	'''
	
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

def translate_barcode_to_tumor(barcode:str) -> str:
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
	raw_alterations 	 = len(rcnv.index) * raw_cases
	
	filtered_cases		 = acnv['case_id'].nunique()
	filtered_alterations = len(acnv.index)

	report = {'Item':['cases',
					  'cases_filtered',
					  'alterations_considered',
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