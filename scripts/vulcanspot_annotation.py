#!/usr/bin/env python3
"""
VulcanSpot annotation module: This script uses a local copy of VulcanSpot 
database (shipped with the pipeline) to annotate genetic alterations present
in VulcanSpot. Vulcan_Local and Vulcan_Pancancer refer to either SLs 
relationships found on cell lines representative of the case's tumor OR to 
SLs relationships found on the pancancer cohor.
"""


import sys

import pandas as pd

def main(input_alterations: str, context_translation: str, vulcan_db: str,
		 gscore_table:str, where_to_save:str):


	tcga_data   = pd.read_csv(input_alterations, sep=',')

	# Get the Gscore of all genes.
	tcga_gscore_annotated = annotate_gscores(gscore_table, tcga_data)

	# Assign TCGA projects to a suitable CCLE context and annotate alterations
	ccle_translated_tcga = translate_tcga_ccle_tissues(context_translation, tcga_gscore_annotated)
	vulcan_annotated_tcga = annotate_vulcan_data(vulcan_db, ccle_translated_tcga)

	del tcga_data
	del ccle_translated_tcga

	vulcan_annotated_tcga.to_csv(where_to_save, sep=',', index=False)



def translate_tcga_ccle_tissues(context_translation: str,
								tcga_data: pd.DataFrame) -> pd.DataFrame:
	'''
	Translates from TCGA project names to the most suitable VulcanSpot context
	based on an external equivalence table built manually.
	'''
	context_translation = pd.read_csv(context_translation, sep='\t')
	context_translation = context_translation.set_index('TCGA')['Vulcanspot'].to_dict()

	tcga_data['Context']      = tcga_data['Project'].map(context_translation)

	return tcga_data

def annotate_gscores(gscore_table:str, tcga_data:pd.DataFrame) -> pd.DataFrame:
	'''
	This method uses the table provided as input to annotate the pandrugs gscore
	of all genes present in input dataframe. Returns the same dataframe with a 
	'gscore' column. 
	'''

	gscore_table = pd.read_csv(gscore_table,
							   sep=',',
							   index_col='Hugo_Symbol',
							   usecols=['Hugo_Symbol', 'Gscore']
							   )

	# map is way faster for this
	gscore_table = gscore_table['Gscore'].to_dict()

	tcga_data['gscore'] = tcga_data['Hugo_Symbol'].map(gscore_table)

	return tcga_data


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
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
