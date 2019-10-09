#!/usr/bin/env python3
"""
Patientes table module: This script merges alterations from filtered MAF and CNV
files on a cases basis. It then annotates all alterations using Cosmic's Cancer
Gene Census and attempts to predict the functional impact of each one based on
VulcanSpot criteria.
"""

import pandas as pd

pd.set_option('mode.chained_assignment', None) #TODO: line get_consensus_from_dupl

def main():
	
	## SNAKEMAKE I/O ##
	input_maf 		= snakemake.input[0]
	input_cnv		= snakemake.input[1]
	cancer_census 	= snakemake.input[2]
	
	where_to_save 		  = snakemake.output[0]
	where_to_save_metrics = snakemake.output[1]

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
	# Classify alterations in GoF/LoF/Unknown
	annotated_tcga['Consequence'] = annotated_tcga.apply(func=annotate_gof_lof,
														 axis='columns'
														 )

	# Resolve conflicts arising from genes with multiple alterations
	tcga_conflicts_resolved = get_consensus_from_duplicates(annotated_tcga)

	# clean inconsistent alterations
	tcga_func_annotated = delete_inconsistent_alterations(tcga_conflicts_resolved)
	
	#TODO: these could be kept depending on user config
	# Drop rows which won't be used anymore
	tcga_func_annotated = tcga_func_annotated.drop(['Chromosome', 'Start_Position',
							      					'End_Position', 'VAF',
													'priority'
							     					],
							 						axis=1,
							  					   )
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
	
	merged_alterations['CGC'] = merged_alterations['Hugo_Symbol'].isin(cancer_gene_census['Gene Symbol'])

	cancer_gene_census['Role in Cancer'].fillna(value='None', inplace=True)
	cancer_gene_census = cancer_gene_census.set_index('Gene Symbol')['Role in Cancer'].to_dict()

	merged_alterations['Role'] = merged_alterations['Hugo_Symbol'].map(cancer_gene_census)
	merged_alterations['Role'].fillna(value='None', inplace=True)
	
	del cancer_gene_census
	return merged_alterations


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
	is_cgc		 = tcga_data['CGC']

	result = 'Unknown'

	if has_cnv:
		if tcga_data['copy_number'] == 'cnv_loss':
			result = 'LoF'
		else: # might be a GoF, check for somatic point mutations
			if not is_truncated and not is_missense:
				result = 'GoF'
			elif not is_truncated and is_missense and is_oncogene:
				result = 'GoF'
			elif is_truncated & (tcga_data['VAF'] >=0.7):
				result = 'LoF'
			else:
				result = 'Unknown'
	else:
		if is_truncated  & (tcga_data['VAF'] >= 0.7):
			result = 'LoF'
		
		if is_missense:
			if is_cgc:
				if is_oncogene & (tcga_data['VAF'] >= 0.2):
					result = 'GoF'
				if not is_oncogene & is_cgc & (tcga_data['VAF'] >= 0.7):
					result = 'LoF'
			else:
				result = 'Unknown'
			
	return result

def get_consensus_from_duplicates(tcga_data: pd.DataFrame) -> pd.DataFrame:
	'''
	Retrieves genes with multiple hits and therefore multiple consequences
	predictions and returns a consensus of all of them. Criteria is
	LoF > GoF > Unknown meaning that a single LoF prediction overrides other
	values and so on
	'''
	
	# Drop cases where all predictions for a gene are the same
	pre_filtered = tcga_data.drop_duplicates(subset=['Hugo_Symbol', 'case_id',
									  				'Project', 'sample',
									  				'Consequence'
									 				], 
							 				 keep='first'
							  				)

	priorities = {'LoF': 2, 'GoF': 1, 'Unknown': 0}

	pre_filtered['priority'] = pre_filtered['Consequence'].map(priorities)
	
	# Get row ID of those alterations with the highest priority after grouping
	alts_to_keep = pre_filtered.groupby(['Hugo_Symbol', 'case_id',
										 'Project', 'sample'
										]
									   )['priority'].transform(max) == pre_filtered['priority']
	
	filtered = pre_filtered[alts_to_keep]
	return filtered


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
	main()