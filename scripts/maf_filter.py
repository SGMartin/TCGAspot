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
	columns_of_interest = ["Hugo_Symbol", "Variant_Classification", "t_depth",
						   "n_depth", "t_ref_count", "t_alt_count", "case_id",
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
    # https://github.com/mskcc/vcf2maf/issues/88
	not_relevant = ["3'UTR", "5'UTR", "Intron", "RNA", "3'Flank", 
					"5'Flank", "Translation_Start_Site", "IGR", "Silent",
					"Splice_Region"
				    ]
	
	mutationData = mutationData[~mutationData['Variant_Classification'].isin(not_relevant)]

	# drop remaining non-exonic mutations
	# TODO: look above link and "Splice_region" stuff
	is_intronic  = mutationData['Consequence'] == 'intron_variant'
	mutationData = mutationData[~is_intronic]

	# drop low depth mutations
	mutationData = mutationData[mutationData['t_depth'] >= 30]
	mutationData = mutationData[mutationData['n_depth'] >= 30]

	# Calculate variant allele frequency for each variant
	mutationData['VAF'] = mutationData['t_alt_count'] / mutationData['t_depth']

	# Drop entries with VAF < 0.2. Vulcan requires this
	mutationData = mutationData[mutationData['VAF'] >= 0.2]

	# Drop depth columns

	mutationData.drop(
					labels=["t_depth", "n_depth", "t_ref_count", "t_alt_count"],
					axis='columns',
					inplace=True
					 )

	return mutationData

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3])

