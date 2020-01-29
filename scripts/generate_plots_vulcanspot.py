#!/usr/bin/env python3
"""
VulcanSpot plots module. This module generate plots related to VulcanSpot
drug prescription.
"""

import matplotlib.pyplot as plt
import pandas  as pd
import numpy   as np
import seaborn as sns

def main():

	## Snakemake I/O ##
	summary 			 = snakemake.input[0]
	vulcan_treatments_db = snakemake.input[1]

	cases_druggable		 		= snakemake.output[0]
	alterations_count_local		= snakemake.output[1]
	alterations_count_pancancer = snakemake.output[2]
	drug_sources				= snakemake.output[3]


	summary = pd.read_csv(summary, sep=',', low_memory=False)

	#TODO Improve this by plotting multiple gscores or generating dynamics graphs
	# Keep gscore >= 0.6 for these plots. Otherwise, there is little meaning...

	summary = summary[summary['gscore'] >= 0.6]

	report_lincs_pandrugs(summary, vulcan_treatments_db, drug_sources)
	report_patients_alterations_boxplot(summary,alterations_count_local, alterations_count_pancancer)
	report_patient_summary(summary, cases_druggable)

def report_patient_summary(summary: pd.DataFrame, where_to_save: str):
	'''
	Builds a report of cases with at least one alteration druggable  using
	synth. lethality as predicted by VulcanSpot.
	'''
	# Building a frequency table
	totals     = summary.groupby('Project')['case_id'].nunique().reset_index()
	localv     = summary[summary['Vulcan_Local']].groupby('Project')['case_id'].nunique()
	pancancer  = summary[summary['Vulcan_Pancancer']].groupby('Project')['case_id'].nunique()
	
	totals.rename({'case_id':'total'}, axis='columns', inplace=True)

	totals['pancancer'] = totals['Project'].map(pancancer)
	totals['local']		= totals['Project'].map(localv).fillna(0)

	#TODO: save this to a table before % and melting... useful for report?
	totals['pancancer'] = totals['pancancer'] / totals['total'] * 100
	totals['local']		= totals['local'] / totals['total'] * 100

	# totals removed as they are always 100
	totals = totals.melt(id_vars='Project',
						 var_name='type',
						 value_name='count',
						 value_vars=['pancancer', 'local']
						)

	# rounding to the units
	totals['count'] = totals['count'].round(0)
	
	# renaming for easier plotting... better than handling matplotlib labels manually
	label_dict = {
				  'pancancer': 'Pan-cancer',
				  'local': 'Matched tissue'
				 }	
	totals['type'] = totals['type'].map(label_dict)

	# Create a new figure
	fig, ax = plt.subplots(figsize=(7,13))
	sns.set_style('whitegrid')

	sns.set_color_codes('pastel')
	figure_pancancer = sns.barplot(x='count',
								   y='Project',
								   data=totals[totals['type'] == 'Pan-cancer'],
								   label='Pan-cancer',
								   color='b'
								   )

	sns.set_color_codes('muted')
	figure_local = sns.barplot(
						 x='count',
						 y='Project',
						 data=totals[totals['type'] == 'Matched tissue'],
						 label='Matched tissue',
						 color='b'
						 )

	# Figure details
	plt.title('Cases druggable by VulcanSpot')

	# Re-order artists and handles so that matched tissue legend bar comes 
	# before pancancer.
	handles, labels = ax.get_legend_handles_labels()
	handles = [handles[1], handles[0]]
	labels  = [labels[1], labels[0]]

	ax.legend(
			  handles=handles,
			  labels=labels,
			  ncol=2,
			  loc='lower center',
			  bbox_to_anchor=(0.5,-0.1),
			  frameon=True
			  )

	ax.set(xlim=(0,100),
	       ylabel='TCGA tumor type', 
	       xlabel='Fraction of cases'
		   )

	sns.despine(left=True)

	plt.savefig(where_to_save, format='svg')



def report_lincs_pandrugs(summary: pd.DataFrame, vulcan_treatments_db:str, where_to_save:str):
	'''
	Builds a report to show how many drugs are prescribed using LINCS vs
	using PanDrugs on a project level. 
	'''

	# Get total patients count
	patients_total  = summary.groupby('Project')['case_id'].nunique().reset_index()

	drug_data = pd.read_csv(vulcan_treatments_db,
							sep=',', 
							usecols=['geneA', 'tissue', 'drug',
									 'alteration','Dscore', 'LINCS']
							)
	
	# Matched tissue first
	matched_alterations = summary[summary['Vulcan_Local']]

	# Drug patients
	drugged_patients = pd.merge(matched_alterations,
							    right=drug_data,
								how='left',
								left_on=['Hugo_Symbol', 'Consequence', 'Context'],
								right_on=['geneA', 'alteration', 'tissue']
								)
	

	drugged_patients['has_dscore'] = drugged_patients['Dscore'] >= 0.6
	drugged_patients['has_lincs']  = drugged_patients['LINCS'] >= 0.9

	drugged_patients.drop(labels=['geneA', 'tissue', 'Consequence',
						   		  'Dscore', 'LINCS', 'Vulcan_Local',
						   		  'Vulcan_Pancancer', 'gscore'],
						   axis=1,
						   inplace=True)
	
	# Build a frequency table
	drugged_dscore = drugged_patients[drugged_patients['has_dscore']].groupby('Project')['case_id'].nunique()
	drugged_lincs  = drugged_patients[drugged_patients['has_lincs']].groupby('Project')['case_id'].nunique()

	patients_total.rename({'case_id':'total'}, axis='columns', inplace=True)

	patients_total['PanDrugs'] = patients_total['Project'].map(drugged_dscore).fillna(0)
	patients_total['Lincs']  = patients_total['Project'].map(drugged_lincs).fillna(0)

	# relative counts and round up
	patients_total['PanDrugs'] = patients_total['PanDrugs'] / patients_total['total'] * 100
	patients_total['Lincs']  = patients_total['Lincs'] / patients_total['total'] * 100 

	# totals removed as they are always 100
	druggable_counts = patients_total.melt(id_vars='Project',
						 				   var_name='Drug prescription',
						 				   value_name='count',
						 				   value_vars=['PanDrugs', 'Lincs']
										   )
	
	druggable_counts['count'] = druggable_counts['count'].round(1)

	# Create a new figure
	fig, ax = plt.subplots(figsize=(7,13))
	sns.set_style('whitegrid')

	sns.set_color_codes('pastel')
	figure_lincs = sns.barplot(x='count',
							   y='Project',
							   data=druggable_counts[druggable_counts['Drug prescription'] == 'Lincs'],
							   label='Lincs',
							   color='b'
							   )

	sns.set_color_codes('muted')
	figure_pandrugs = sns.barplot(
						 x='count',
						 y='Project',
						 data=druggable_counts[druggable_counts['Drug prescription'] == 'PanDrugs'],
						 label='PanDrugs',
						 color='b'
						 )

	# Figure details
	plt.title('Cases druggable by VulcanSpot: drug sources')

	# Re-order artists and handles so that matched tissue legend bar comes 
	# before pancancer.
	handles, labels = ax.get_legend_handles_labels()
	handles = [handles[1], handles[0]]
	labels  = [labels[1], labels[0]]

	ax.legend(
			  handles=handles,
			  labels=labels,
			  ncol=2,
			  loc='lower center',
			  bbox_to_anchor=(0.5,-0.1),
			  frameon=True
			  )

	ax.set(xlim=(0,100),
	       ylabel='TCGA tumor type', 
	       xlabel='Fraction of cases'
		   )

	sns.despine(left=True)
	plt.savefig(where_to_save, format='svg')
	

#TODO: refactor this in two methods?
def report_patients_alterations_boxplot(summary: pd.DataFrame, save_local:str, save_pancancer:str):
	'''
	This method generates two plots: one counts the average number of druggable
	alterations per patient in a matched tissue context and the other performs
	the same count on a pan-cancer context.
	'''

	patients_alterations = summary.drop(['Variant_Classification', 'Role',
										 'Context', 'copy_number','sample'
										],
										axis=1
										)
	
	in_local	 = patients_alterations['Vulcan_Local']
	in_pancancer = patients_alterations['Vulcan_Pancancer']


	patients_alterations = patients_alterations[in_local | in_pancancer]


	patients_alterations = pd.pivot_table(data=patients_alterations,
										  values=['Vulcan_Local','Vulcan_Pancancer'],
										  index=['case_id', 'Project'],
										  aggfunc=np.sum

									)
	
	patients_alterations.reset_index(inplace=True)

	patients_alterations.columns = ['case_id', 'project', 'local', 'pancancer']

	alterations	 = pd.melt(frame=patients_alterations,
						   id_vars=['case_id', 'project'],
						   value_vars=['local', 'pancancer'],
						   var_name='context',
						   value_name='count'
						  )
		
	#Let's plot... too big to hue? TODO: Axes too different in scale to hue.
	local_alt	  = alterations[alterations['context'] == 'local'].sort_values('project', ascending=True)
	pancancer_alt = alterations[alterations['context'] == 'pancancer'].sort_values('project', ascending=True)

	boxplot_alterations(data=local_alt,
						context='matched tissue',
						save_to=(save_local)
						)
	
	boxplot_alterations(data=pancancer_alt,
						context='pan-cancer',
						save_to=(save_pancancer)
						)


def boxplot_alterations(data: pd.DataFrame, context:str, save_to:str):

	fig, ax = plt.subplots(figsize=(15, 7))

	sns.set_style('whitegrid')
	sns.set_color_codes('pastel')

	boxplot = sns.barplot(x='project',
						  y='count',
						  data=data,
						  color='b'
						 )

	plt.title(f"Summary of patients druggable alterations: {context}")
	ax.set(xlabel='TCGA tumor type', ylabel='Patients druggable alterations')

	boxplot.set_xticklabels(boxplot.get_xticklabels(), rotation=30)
	plt.savefig(save_to, format='svg')


if __name__ == "__main__":
	main()