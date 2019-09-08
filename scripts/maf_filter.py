#!/usr/bin/env python3
"""
Section: Benchmarking VulcanSpot using TCGA data

Purpose: To assess VulcanSpot performance against a public dataset of cancer 
patients genomic data. 
"""

import os.path
import sys

import numpy  as np
import pandas as pd

def main(input_maf:str, alterations_save_as:str, metrics_save_as:str):

	raw_mutation_data =  input_maf
	where_to_save 	  =  alterations_save_as
	metrics			  =  metrics_save_as
	

	mutation_data = load_dataframe(raw_mutation_data)
	filtered_maf  = filter_quality_variants(mutation_data)

	# Check if there is an annotations.txt file present nearby
	# TODO: Log this
	input_maf_directory = os.path.dirname(input_maf)
	annotation_file     = input_maf_directory + '/annotations.txt'

	if(os.path.exists(annotation_file)):
		filtered_maf = exclude_annotated_cases(filtered_maf, annotation_file)

	# Translate barcode to sample type
	filtered_maf['sample'] = translate_barcode_to_sample_type(filtered_maf['Tumor_Sample_Barcode'])
	filtered_maf.drop('Tumor_Sample_Barcode', axis=1, inplace=True)

	# SAVE OUTPUT #
	filtered_maf.to_csv(where_to_save, sep=',', index=False)
	
	#TODO: This should be a method
	#  Report metrics #
	raw_cases            = mutation_data['case_id'].nunique()
	raw_alterations 	 = len(mutation_data.index)
	filtered_cases       = filtered_maf['case_id'].nunique()
	filtered_alterations = len(filtered_maf.index)

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
	report.to_csv(metrics, sep=',', index=False)


def load_dataframe(input_file: str) -> pd.DataFrame:
	"""
	Loads input TCGA mutect .maf file as pandas dataframe, appends a project 
	column and returns it.
	"""
	columns_of_interest = ["Hugo_Symbol", "Variant_Classification", 
						   "Tumor_Sample_Barcode", "t_depth", "n_depth",
						   "t_ref_count", "t_alt_count", "case_id",
						   "Consequence"
						  ]

	split_project = input_file.split('.') #get project name

	# search for project name based on relative location of "mutect"
	# note that this will fail depending on dir names... but it's the
	# best we can do
	# SCAN for MUTECT, VARSCAN etc...
	# 	
	p_index 	  = split_project.index('mutect')
	project_index = p_index - 1
	project_name   = split_project[project_index]

	raw_maf  = pd.read_csv(input_file,
						   sep='\t',
						   header=5,
						   usecols=columns_of_interest,
						   low_memory=False
						   )

	raw_maf['Project'] =  project_name

	return raw_maf

def filter_quality_variants(mutationData: pd.DataFrame) -> pd.DataFrame:
	"""
	Filters variants based on alteration consequence and sequencing depth.
	Returns a dataframe with coding only alterations with sufficient ( > 30) 
	depth.
	"""

	# drop non coding, intronic and silent mutations
	not_relevant = ["3'UTR", "5'UTR", "Intron", "RNA", "3'Flank", 
					"5'Flank", "Translation_Start_Site", "IGR", "Silent"
				   ]
	
	mutationData = mutationData[~mutationData['Variant_Classification'].isin(not_relevant)]

	# drop remaining non-exonic mutations
	# look https://github.com/mskcc/vcf2maf/issues/88 for more info about this
	# step
	is_intronic  = mutationData['Consequence'] == 'intron_variant'
	mutationData = mutationData[~is_intronic]

	# drop low depth mutations
	mutationData = mutationData[mutationData['t_depth'] >= 30]
	mutationData = mutationData[mutationData['n_depth'] >= 30]

	# Calculate variant allele frequency for each variant
	mutationData['VAF'] = mutationData['t_alt_count'] / mutationData['t_depth']

	# Drop entries with VAF < 0.2. Vulcan requires this
	mutationData = mutationData[mutationData['VAF'] >= 0.2]

	# Drop unused columns

	mutationData.drop(
					labels=["t_depth", "n_depth", "t_ref_count",
							"t_alt_count", 'Consequence'],
					axis='columns',
					inplace=True
					 )

	return mutationData

#TODO: make this customizable
def exclude_annotated_cases(maf: pd.DataFrame, annotations: str) -> pd.DataFrame:
	'''
		Loads annotations.txt file for this project and returns a dataframe
		of excluded cases id. Cases/Items in annotations.txt are spared	in 
		certain cases based on notification category.
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

	excluded_cases   = maf['case_id'].isin(cases_to_exclude['entity_id'])
	filtered_maf	 = maf[~excluded_cases]

	return filtered_maf

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


if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3])

