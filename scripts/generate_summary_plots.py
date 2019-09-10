#!/usr/bin/env python3

import sys

import matplotlib.pyplot as plt
import numpy    as np
import pandas  as pd 
import seaborn as sns

def main(summary: str, save_to:str):

	druggable_cases     = save_to + '/cases_druggable.svg'
	gof_lof_alterations = save_to + '/alterations_classified.svg'
	alterations_count 	= save_to + '/alterations_count'
	


	summary = pd.read_csv(summary, sep=',', low_memory=False)

	report_patients_alterations_boxplot(summary,alterations_count)
	report_patient_summary(summary, druggable_cases)
	report_gof_lof_alterations(summary, gof_lof_alterations)


def report_gof_lof_alterations(summary:pd.DataFrame ,where_to_save: str):

	# Create a new fig
	fig  = plt.figure(figsize=(15,6))
	sns.set_color_codes('muted')
	countplot = sns.countplot(x='Consequence',
							  data=summary,
							  order=summary['Consequence'].value_counts().index
							 )

	countplot.set(xlabel='Alterations classified', ylabel='Count')

	plt.savefig(where_to_save, format='svg')

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

#TODO: refactor this in two methods?
def report_patients_alterations_boxplot(summary: pd.DataFrame, where_to_save:str):

	patients_alterations = summary.drop(['Variant_Classification', 
										'VAF','Role', 'Context',
										'copy_number','sample'
										],
										axis=1)
	
	in_local	 = patients_alterations['Vulcan_Local']
	in_pancancer = patients_alterations['Vulcan_Pancancer']


	patients_alterations = patients_alterations[in_local | in_pancancer]

	# Drop duplicates: genes with multiple hits. We'll consider them a single
	# target
	patients_alterations.drop_duplicates(inplace=True)

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
						save_to=(where_to_save + '_local.svg')
						)
	
	boxplot_alterations(data=pancancer_alt,
						context='pan-cancer',
						save_to=(where_to_save + '_pancancer.svg')
						)


def boxplot_alterations(data: pd.DataFrame, context:str, save_to:str):

	fig, ax = plt.subplots(figsize=(15, 7))

	sns.set_style('whitegrid')
	sns.set_color_codes('pastel')

	boxplot = sns.boxplot(x='project',
						  y='count',
						  data=data
						 )

	plt.title(f"Summary of {context} patients druggable alterations")
	ax.set(xlabel='Project', ylabel='Patients druggable alterations')

	boxplot.set_xticklabels(boxplot.get_xticklabels(), rotation=30)
	plt.savefig(save_to, format='svg')

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])




