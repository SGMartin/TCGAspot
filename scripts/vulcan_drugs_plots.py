#!/usr/bin/env python3
'''
Vulcan drugs plots module: This scripts uses VulcanSpot drug data to draw 
different plots. TODO:
'''

import sys
import pandas as pd


def main(summary: str, vulcan_treatments_db: str):

	# Build a drugs dataframe for the plots
	drug_summary = pd.read_csv(summary, sep=',', usecols=['case_id','Project',
														  'Hugo_Symbol',
														  'Consequence',
														  'Context',
														  'Vulcan_Local',
														  'Vulcan_Pancancer'
														  ])
	
	vulcan_drugs_db = pd.read_csv(vulcan_treatments_db, sep=',')

	# We have to split this operations in two parts to keep mem. usage low
	# Add matched tissue drugs first
	drugs_local_summary = pd.merge(left=drug_summary[drug_summary['Vulcan_Local']],
								   right=vulcan_drugs_db,
								   how='left',
								   left_on=['Hugo_Symbol', 'Consequence', 'Context'],
								   right_on=['geneA', 'alteration', 'tissue']
								   )

	drugs_local_summary['has_dscore'] = drugs_local_summary['Dscore'] > 0
	drugs_local_summary['has_lincs']  = drugs_local_summary['LINCS'] > 0

	tcga_drug_summary = drugs_local_summary.groupby('Project')['has_dscore',
															   'has_lincs'
															  ].sum()

	tcga_drug_summary 		  = pd.DataFrame(data=tcga_drug_summary,
											 columns=['pandrugs_based_local',
											 		  'lincs_based_local']
											)

	
	# delete this object to free allocated mem.
	del drugs_local_summary

	drugs_pancancer_summary = pd.merge(left=drug_summary[drug_summary['Vulcan_Pancancer']],
									   right=vulcan_drugs_db[vulcan_drugs_db['tissue'] == 'PANCANCER'],
									   how='left',
									   left_on=['Hugo_Symbol', 'Consequence'],
									   right_on=['geneA', 'alteration']
									   )
	

	drugs_pancancer_summary['has_dscore'] = drugs_pancanser_summary['Dscore'] > 0
	drugs_pancancer_summary['has_lincs']  = drugs_pancancer_summary['LINCS'] > 0

	tcga_drug_pancancer_summary = drugs_pancancer_summary.groupby('Project')['has_dscore', 'has_lincs'].sum()

	del drugs_pancancer_summary

	tcga_drug_pancancer_summary = pd.DataFrame(data=tcga_drug_pancancer_summary,
											   columns=['pandrugs_based_pancancer',
											   			'lincs_based_pancancer'
													   ]
											   )






if __name__ == "__main__":
	main(sys.argv[1:])