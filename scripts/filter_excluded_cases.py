#!/usr/bin/env python3

import sys

import pandas as pd

def main(rdata, annot, save_as, metrics_save_as):

	### INPUT ###
	data_to_filter    = rdata
	input_annotations = annot

	### OUTPUT ###
	where_to_save 	= save_as
	report			= metrics_save_as

	raw_alterations  = pd.read_csv(data_to_filter, sep=',')
	
	# Excluding cases
	cases_to_exclude = load_excluded_cases(input_annotations)
	excluded_cases   = raw_alterations['case_id'].isin(cases_to_exclude['entity_id'])

	filtered_cases = raw_alterations[~excluded_cases]

	# get metrics
	metrics = generate_report(raw_alterations, filtered_cases)

	# Save to file
	filtered_cases.to_csv(where_to_save, sep=',', index=False)
	metrics.to_csv(report, sep=',', index=False)


def load_excluded_cases(input_annotation: str) -> pd.DataFrame:
	'''Loads annotations.txt file for this project and returns a dataframe
		of excluded cases id. Cases/Items in annotations.txt are spared
		in certain cases based on notification category.
	'''

	categories_to_keep = ['History of acceptable prior treatment related to a prior/other malignancy',
						  'Acceptable treatment for TCGA tumor',
						  'Alternate sample pipeline'
						 ]

	raw_annotations = pd.read_csv(input_annotation,
								  sep='\t',
								  usecols=['entity_id', 'category']
								  )

	acceptable_alterations = raw_annotations['category'].isin(categories_to_keep)
	filtered_annotations   = raw_annotations[~acceptable_alterations]

	return filtered_annotations


def generate_report(r_alterations: pd.DataFrame, 
					c_alterations: pd.DataFrame) -> pd.DataFrame:
	'''
	Generates some metrics for logging.
	'''
	
	raw_cases 		= r_alterations['case_id'].nunique()
	raw_alterations = len(r_alterations.index)

	filtered_cases 		 = c_alterations['case_id'].nunique()
	filtered_alterations = len(c_alterations.index)

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
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


