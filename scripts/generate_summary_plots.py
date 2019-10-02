#!/usr/bin/env python3
'''
Summary plot module: This scripts uses the seaborn package to generate various
plots regarding alterations kept after all the pipeline has finished.
'''

import sys

import matplotlib.pyplot as plt
import numpy    as np
import pandas  as pd 
import seaborn as sns

def main(summary: str, save_to:str):

	summary = pd.read_csv(summary,
						  sep=',',
						  usecols=['Hugo_Symbol', 'Variant_Classification',
						  		   'copy_number', 'sample',
								   'gscore','Consequence']
						  )

	ranges = [0, 0.2, 0.4, 0.6, 0.8]

	for i in ranges:	# it is really the fastest way to do this... geez
		
		alt_plots = summary.loc[summary['gscore'] >= i]
		plot_name = f"{save_to}/alterations_classified_{i}.svg"
		report_gof_lof_alterations(alt_plots.copy(), plot_name)


def report_gof_lof_alterations(summary:pd.DataFrame, where_to_save: str):
	'''
	This methods reports the relative proportion of GoF, LoF and Unknown
	alterations present in the dataframe. Saves a barplot to where_to_save
	argument.
	'''
	alts = summary

	alts['snv']  = alts['Variant_Classification'] != 'None'
	alts['scna'] = alts['copy_number'] != 'None'
	alts['both'] = np.logical_and(alts['snv'], alts['scna']) 

	alts.drop(['Variant_Classification', 'copy_number'], axis=1, inplace=True)

	alts = alts.groupby('Consequence').agg(np.sum) # counts True as 1
	alts = alts.reset_index()
	alts.drop('gscore', axis=1, inplace=True)

	alts = pd.melt(frame=alts,
				   id_vars='Consequence',
				   value_vars=['scna', 'snv', 'both'],
				   var_name='source',
				   value_name='count'
				   )
	
	# renaming for a better plot legend
	readable_names = {'snv':'SNV', 'scna': 'CNV', 'both': 'Both'}
	alts['Source'] = alts['source'].map(readable_names)

	alts.drop('source', axis=1, inplace=True)
	
	# Get a % out of total count
	alts['count'] = alts['count'] / alts['count'].sum() * 100
	
	# Create a new fig
	sns.set(style='whitegrid')
	sns.set_color_codes('muted')

	fig = plt.figure(figsize=(15,6))

	plt.title('Alterations classified by source')
	barplot = sns.catplot(x='Consequence', y='count',
						  data=alts, hue='Source', 
						  kind='bar')
	
	barplot.set(xlabel='Gene alterations', ylabel="% of total alterations")

	plt.savefig(where_to_save, format='svg')

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])




