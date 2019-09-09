#!/usr/bin/env python3

import sys

import matplotlib.pyplot as plt
import pandas  as pd 
import seaborn as sns

def main(summary: str, save_to:str):

	report_patient_summary(summary, save_to)



def report_patient_summary(summary: str, where_to_save: str):
	'''
	Builds a report of cases with at least one alteration druggable  using
	synth. lethality as predicted by VulcanSpot.
	'''

	summary = pd.read_csv(summary, sep=',', usecols=['case_id', 'Project', 
													 'Vulcan_Local', 'Vulcan_Pancancer'
													])
	
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
	totals['total']	= 100


	totals = totals.melt(id_vars='Project',
						 var_name='type',
						 value_name='count',
						 value_vars=['total', 'pancancer', 'local']
						)

	# rounding to the units
	totals['count'] = totals['count'].round(0)

	# renaming for easier plotting... better than handling matplotlib labels manually
	label_dict = {'total' : 'All cases',
				  'pancancer': 'Pan-cancer',
				  'local': 'Matched tissue'
				 }
	
	totals['type'] = totals['type'].map(label_dict)

	# Create a new figure
	fig, ax = plt.subplots(figsize=(15,6))

	sns.set_style('whitegrid')
	sns.set_color_codes('pastel')

	figure = sns.barplot(x='count',
						 y='Project',
						 hue='type',
						 data=totals
						 )
	
	# Figure details
	plt.title('Cases druggable by VulcanSpot')

	ax.legend(ncol=1,
			  loc='lower right',
			  frameon=True
			  )

	ax.set(xlim=(0,100),
	       ylabel='TCGA project', 
	       xlabel='Fraction of cases'
		   )

	# Save it
	plt.savefig(where_to_save)

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])




