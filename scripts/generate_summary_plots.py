#!/usr/bin/env python3

import sys

import matplotlib.pyplot as plt
import pandas  as pd 
import seaborn as sns

def main(summary: str, save_to:str):

	druggable_cases     = save_to + '/cases_druggable.svg'
	gof_lof_alterations = save_to + '/alterations_classified.svg'


	summary = pd.read_csv(summary, sep=',', low_memory=False)

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
	fig, ax = plt.subplots(figsize=(6,15))
	

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

	ax.legend(ncol=2,
			  loc='lower center',
			  bbox_to_anchor=(0.5,-0.1),
			  frameon=True
			  )

	ax.set(xlim=(0,100),
	       ylabel='TCGA tumor type', 
	       xlabel='Fraction of cases'
		   )

	sns.despine(left=True)
	# Save it
	plt.savefig(where_to_save, format='svg')

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])




