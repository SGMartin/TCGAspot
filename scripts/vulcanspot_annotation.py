#!/usr/bin/env python3

"""
Section: SL in TCGA data

Purpose: Annotate merged CNV and maf with data from Cancer Gene Census
and SL identified by VulcanSpot.
"""

import sys

import pandas as pd

def main(input_alterations: str, vulcan_db: str,
		 where_to_save:str):


	tcga_data   = pd.read_csv(input_alterations, sep=',')
	
	vulcan_annotated_tcga = annotate_vulcan_data(vulcan_db, tcga_data)

	vulcan_annotated_tcga.to_csv(where_to_save, sep=',', index=False)



def annotate_vulcan_data(vulcan_table: str,
						 tcga_data: pd.DataFrame) -> pd.DataFrame:
	'''
	Attempts to match TCGA alterations to VulcanSpot database. Returns a 
	dataframe with two appended columns, Vulcan_Local and Vulcan_Pancancer.
	These columns values are boolean and represent matches in a specific context
	or against Pancancer cohort.
	'''

	# Load available genes in vulcanSpot
	vulcan_genes = pd.read_csv(vulcan_table,
							   sep=',',
							   usecols=['geneA', 'tissue', 'alteration'],
							   low_memory=False
							   ) 

	vulcan_genes.drop_duplicates(keep='first', inplace=True)
	
	# Annotate local vulcan alterations #
	tcga_data = pd.merge(tcga_data, vulcan_genes,
						  how='left',
						  left_on=['Hugo_Symbol', 'Context', 'Consequence'], 
						  right_on=['geneA', 'tissue', 'alteration'], 
						  indicator=True
						)

	tcga_data.drop(['tissue', 'geneA', 'alteration'], inplace=True, axis=1)

	both_to_bool = { 'both': True, 'left_only' : False}	
	tcga_data['Vulcan_Local'] = tcga_data['_merge'].map(both_to_bool)
	tcga_data.drop('_merge', inplace=True, axis=1)
	
	# annotate global vulcan alterations
	#TODO: Refactor this
	vulcan_pancancer_genes = vulcan_genes[vulcan_genes['tissue'] == 'PANCANCER']

	tcga_data = pd.merge(tcga_data, vulcan_pancancer_genes,
						  how='left',
						  left_on=['Hugo_Symbol', 'Consequence'], 
						  right_on=['geneA', 'alteration'], 
						  indicator=True
						)

	tcga_data.drop(['tissue', 'geneA', 'alteration'], inplace=True, axis=1)

	both_to_bool = { 'both': True, 'left_only' : False}	
	tcga_data['Vulcan_Pancancer'] = tcga_data['_merge'].map(both_to_bool)
	tcga_data.drop('_merge', inplace=True, axis=1)

	del vulcan_genes

	return tcga_data


if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3])
