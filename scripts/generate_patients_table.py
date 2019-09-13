#!/usr/bin/env python3
"""
Section: SL in TCGA data

Purpose: Generate a table of genetic alterations based on CNV and MAF data
from the TCGA. Patients genetic landscape is inferred from these sources.
"""
import sys

import pandas as pd

def main(input_maf:str, input_cnv:str, cancer_census:str, 
	     context_translation:str, where_to_save: str, where_to_save_metrics:str):

	# Loading input dataframes #
	maf_to_merge = pd.read_csv(input_maf,
							   sep=',',
							   usecols=['Hugo_Symbol', 'Chromosome',
							   		    'Start_Position', 'End_Position',
							    		'Variant_Classification',
							   			'VAF', 'case_id', 'sample', 'Project'
										],
							   low_memory=False
							   )

	cnv_to_merge = pd.read_csv(input_cnv,
							   sep=',',
							   usecols=['Gene Symbol', 'case_id',
							   		    'copy_number', 'sample'
									   ],
							   low_memory=False
							   )

	merged_alterations = merge_cnv_to_mutations(maf_to_merge, cnv_to_merge)

	# Annotate with CGC and translate projects to CCLE tissues
	annotated_tcga			 = annotate_cancer_gene_census(cancer_census,
												   	       merged_alterations
														   )
	context_translated_tcga  = translate_tcga_ccle_tissues(context_translation,
														   annotated_tcga
														  )
	# Classify alterations in GoF/LoF/Unknown
	context_translated_tcga['Consequence'] = context_translated_tcga.apply(func=annotate_gof_lof,
																		   axis='columns'
																		   )
	# clean inconsistent alterations
	tcga_func_annotated = delete_inconsistent_alterations(context_translated_tcga)

	# SAVING #
	tcga_func_annotated.to_csv(where_to_save, sep=',', index=False)
	
	# METRICS #
	metrics = generate_metrics(merged_alterations, tcga_func_annotated)
	metrics.to_csv(where_to_save_metrics, sep=',', index=False)

def merge_cnv_to_mutations(maf: pd.DataFrame, cnv: pd.DataFrame) -> pd.DataFrame:
	'''
	Merges maf and cnv together by outer join and then fills missing rows
	with None instead of NaN.
	'''
	# rename is faster than copying missing genes from one column to another
	# after merging.

	renamed_cnv = cnv.rename({'Gene Symbol' : 'Hugo_Symbol'}, axis='columns')

	merged_maf_cnv =  maf.merge(
							right=renamed_cnv,
					        how='outer',
						    left_on=['case_id', 'sample', 'Hugo_Symbol'],
							right_on=['case_id','sample', 'Hugo_Symbol']
							)
		
	# Fill NaN in Project column for cnv
	cases_projects = maf.loc[:,['case_id', 'Project']]
	cases_projects = cases_projects[~cases_projects.duplicated()] # keeps only first case entry

	fill_dict = cases_projects.set_index('case_id')['Project'].to_dict()

	merged_maf_cnv['Project'] = merged_maf_cnv['case_id'].map(fill_dict)

	# Fill nas in copy_number and variant_classification due to merging.
	merged_maf_cnv['Variant_Classification'].fillna(value='None', inplace=True)
	merged_maf_cnv['copy_number'].fillna(value='None', inplace=True)

	return merged_maf_cnv

def annotate_cancer_gene_census(cancer_census, merged_alterations) -> pd.DataFrame:
	"""
	Annotate roles based on cancer gene census.
	"""
	
	cancer_gene_census = pd.read_csv(cancer_census,
									 sep='\t',
									 usecols=['Gene Symbol', 'Role in Cancer']
								    )
	
	cancer_gene_census['Role in Cancer'].fillna(value='None', inplace=True)
	cancer_gene_census = cancer_gene_census.set_index('Gene Symbol')['Role in Cancer'].to_dict()

	merged_alterations['Role'] = merged_alterations['Hugo_Symbol'].map(cancer_gene_census)
	merged_alterations['Role'].fillna(value='None', inplace=True)
	
	del cancer_gene_census
	return merged_alterations


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


def annotate_gof_lof(tcga_data: pd.Series) -> str:
	'''
	Uses Variant_Classification columns, copy number column when available
	and VAF to predict alteration consequence based on VulcanSpot criteria.
	Returns a dataframe where all alterations are classified either in GoF,
	LoF or Unknown.
	'''

	truncated_protein = ['De_novo_Start_OutOfFrame', 'Frame_Shift_Del',
						 'Frame_Shift_Ins', 'In_Frame_Del','In_Frame_Ins',
						 'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site',
						 'Start_Codon_Del', 'Start_Codon_Ins',
						 'Stop_Codon_Del', 'Stop_Codon_Ins'
						]
	
	is_truncated = tcga_data['Variant_Classification'] in truncated_protein
	is_missense  = tcga_data['Variant_Classification'] == 'Missense_Mutation'
	has_cnv 	 = tcga_data['copy_number'] != 'None'
	is_oncogene  = 'oncogene' in tcga_data['Role']

	result = 'Unknown'

	if has_cnv:
		if tcga_data['copy_number'] == 'cnv_loss':
			result = 'LoF'

		else: # might be a GoF, check for somatic point mutations
			if (is_truncated == False):
				result = 'GoF'
			else:
				result = 'LoF'
	else:
		if is_truncated  & (tcga_data['VAF'] >= 0.7):
			result = 'LoF'
		
		if is_missense:
			if is_oncogene & (tcga_data['VAF'] >= 0.2):
				result = 'GoF'
			else:
				if tcga_data['VAF'] >= 0.7:
					result = 'LoF'

	return result


def delete_inconsistent_alterations(tcga_data: pd.DataFrame) -> pd.DataFrame:
	'''
	Delete rows where a GoF was predicted on a gene annotated as tumor suppressor
	gene. These although rare (<1%) might happen due to how TCGA GISTIC-2 data
	was generated.
	'''

	is_GoF = tcga_data['Consequence'] == 'GoF'
	is_TSG = tcga_data['Role'].str.contains('TSG')

	tcga_cleared = tcga_data[~(is_GoF & is_TSG)]

	return tcga_cleared

#TODO: This could be refactored taking into account the above method.

def generate_metrics(raw_tcga_data: pd.DataFrame,
			         clean_tcga_data: pd.DataFrame) -> pd.DataFrame:
	
	'''
	Generates metrics for this step. Returns a dataframe of cases and alterations
	count, both before and after filtering.
	'''
	
	# count inconsistent alterations
	is_GoF = raw_tcga_data['Consequence'] == 'GoF'
	is_TSG = raw_tcga_data['Role'].str.contains('TSG')

	raw_inconsistent_alterations = len(raw_tcga_data[(is_GoF & is_TSG)].index)
	raw_merged_alterations		 = len(raw_tcga_data.index)
	filtered_merged_alterations  = len(clean_tcga_data.index)

	raw_merged_cases       = raw_tcga_data['case_id'].nunique()
	raw_inconsistent_cases = raw_tcga_data[(is_GoF & is_TSG)]['case_id'].nunique()
	filtered_merged_cases  = clean_tcga_data['case_id'].nunique()

	report = {'Item':['merged_cases', 'merged_alterations',
					  'inconsistent_cases', 'inconsistent_alterations',
					  'filtered_cases', 'filtered_alterations'
					 ],
			  'Count':[raw_merged_cases, raw_merged_alterations,
			   		   raw_inconsistent_cases, raw_inconsistent_alterations,
			  		   filtered_merged_cases, filtered_merged_alterations
					  ]
			 }

	report = pd.DataFrame(report)
	
	return report


if __name__ == "__main__":
	main(
	     sys.argv[1], sys.argv[2], sys.argv[3],
		 sys.argv[4], sys.argv[5], sys.argv[6]
		 )
