#!/usr/bin/env python3
"""
Gain of function correction: This scripts filters the table output by
cnv_de_expr.R to check which cnv-GoF events were indeed overexpressed. Those
not passing the threshold are reclassified as Unknown.
"""

import pandas as pd
pd.set_option('mode.chained_assignment', 'raise') #TODO: line get_consensus_from_dupl


#TODO: This module should grow once metrics and other checks are added. 

def main():

	## SNAKEMAKE I/O ## 
	cases_table		= snakemake.input[0]
	expression_data = snakemake.input[1]
	where_to_save 	= snakemake.output[0]

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
	
	# Get rid of those cnv events with no significance or little log fold
	# TODO: use cases and bias for more accurate class.

	is_significant 		   = expression['adj_pval'] <= 0.05
	is_expressed 		   = expression['log2fc'] >= 1

	genes_verified = expression.loc[(is_significant & is_expressed), 'Hugo_Symbol']

	# Select cases whose GoF class. comes from cnv only
	is_gof   = cases['Consequence'] == 'GoF'
	is_cnv   = cases['copy_number'] == 'cnv_gain'
	no_other = cases['Variant_Classification'] == 'None'
	is_valid = cases['Hugo_Symbol'].isin(genes_verified)

	#Set their values
	cases.loc[(is_gof & is_cnv & no_other), 'Consequence'] = 'Unknown'
	cases.loc[(is_gof & is_cnv & no_other & is_valid), 'Consequence'] = 'GoF'

	del expression
	return cases


if __name__ == "__main__":
	main()
