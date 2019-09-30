#!/usr/bin/env python3
"""
Gain of function correction: This scripts attempts to check wether genes with
a GoF predicted by a cnv_gain event are indeed overexpressed downstream. To
do that it uses filtered expression data. In summary, a simple lookup is 
performed.
"""

import sys

import pandas as pd

#TODO: This module should grow once metrics and other checks are added. 

def main(cases_table: str, expression_data: str, where_to_save:str):
	'''

	'''

	# Load cases table for this project and its expression table
	cases_table 	 = pd.read_csv(cases_table, sep=',')
	expression_table = pd.read_csv(expression_data, sep=',')

	corrected_alterations = match_gof_to_expression(cases_table, expression_table)

	corrected_alterations.to_csv(where_to_save, sep=',', index=False)

def match_gof_to_expression(cases: pd.DataFrame, expression: pd.DataFrame) -> pd.DataFrame:
	'''
	This methods checks wether GoF alterations originated from cnv_gain events
	are indeed present downstream in expression data. Returns a dataframe of
	corrected events.
	'''
	# get rid of expression data not present in this project GoF set
	genes_of_interest = cases[cases['Consequence'] == 'GoF']['Hugo_Symbol'].unique()
	cases_of_interest = cases['case_id'].unique()

	expression = expression[expression['Hugo_Symbol'].isin(genes_of_interest)]
	expression = expression[expression['case_id'].isin(cases_of_interest)]	

	# filter by aliquot A... for now
	expression = expression[expression['aliquot'] == 'A']
	expression.drop('aliquot', axis=1, inplace=True)

	# select GoF of interest
	gofs	   = cases['Consequence'] == 'GoF'
	absent_snv = cases['Variant_Classification'] == 'None'
	has_cnv    = cases['copy_number'] == 'cnv_gain'

	alterations_to_check = cases[gofs & absent_snv & has_cnv]

	# perform a left join to check wether they are indeed overexpr.
	overexpressed = pd.merge(left=alterations_to_check,
							 right=expression,
							 how='left',
							 on=['Hugo_Symbol', 'case_id', 'sample'],
							 suffixes=('','_2'),
							 indicator=True
							)

	# keep alterations present ONLY in expr. data
	overexpressed = overexpressed[overexpressed['_merge'] == 'both']
	overexpressed.drop(labels=['_merge'],
					   axis=1,
					   inplace=True
					  )

	# remove checked alterations from original data
	cases = cases[~(gofs & absent_snv & has_cnv)]
	cases = cases.append(overexpressed, ignore_index=True)

	return cases


if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3])
