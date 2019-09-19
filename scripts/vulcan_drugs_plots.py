#!/usr/bin/env python3
'''
Vulcan drugs plots module: This scripts uses VulcanSpot drug data to draw 
different plots. TODO:
'''

import sys
import pandas as pd


def main(summary: str, vulcan_treatments_db: str):

	tcga_summary = pd.read_csv(summary, sep=',', usecols=['geneA', 'tissue',
							   							  'alteration', 'drug',
												          'Dscore', 'LINCS'
														  ])


def plot_lincs_and_pandrugs(summary: pd.DataFrame, vulcan_treatments_db:str):
	'''
	Plots the amount of drugs prescribed using PanDrugs vs those coming from
	LINCS at the project level.
	'''

	vulcan_treatments = pd.read_csv(vulcan_treatments_db, sep=',')

	# Count patients by TCGA project
	patient_count = summary.groupby('Project')['case_id'].nunique()

	# Get all druggable alterations
	summary_druggable = summary[summary['Vulcan_Local'] | summary['Vulcan_Pancancer']]

	summary_druggable.drop(['Variant_Classification', 'sample',
							'copy_number', 'Role',
							], 
							axis=1,
							inplace=True
							)
	
	# Select treatments for druggable, affected genes.
	genes_of_interest = summary_druggable['Hugo_Symbol'].unique()
	vulcan_treatments = vulcan_treatments[vulcan_treatments['geneA'].isin(genes_of_interest)]

	patients_matched_tissue_treatments = pd.merge(left=summary_druggable,
							    				  right=vulcan_treatments,
												  how='left',
												  left_on=['Hugo_Symbol', 'Consequence', 'Context'],
												  right_on=['geneA', 'alteration', 'tissue']
												  )
	



	





if __name__ == "__main__":
	main(sys.argv[1:])