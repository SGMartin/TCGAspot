#!/usr/bin/env python3
"""
PanDrugs Annotation module:
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns


#TODO: This could be added to vulcanspot plots
def main():
	'''

	'''
	## SNAKEMAKE INPUT/OUTPUT ##
	summary	      = snakemake.input[0]
	pandrugs_db   = snakemake.input[1]
	where_to_save = snakemake.output[0]	

	pandrugs_db = pd.read_csv(pandrugs_db,
							  sep='\t',
							  usecols=['gene_symbol', 'show_drug_name', 'status',
							  		   'target_marker', 'resistance', 'score']
							  )
	

	# Filter gene-drug relationships of interest
	high_score    = pandrugs_db['score'] >= 0.6
	sensible      = pandrugs_db['resistance'] == 'sensitivity'
	direct_target = pandrugs_db['target_marker'] == 'target'
	
	pandrugs_db = pandrugs_db[high_score & sensible & direct_target]

	# Get drugs which are approved and in clinical trials
	CT_drug   	  = pandrugs_db['status'] == 'Clinical trials'
	FDA_drug	  = pandrugs_db['status'] == 'Approved'	

	# Load TCGA data
	tcga_data = pd.read_csv(summary,
							sep=',',
							usecols=['Hugo_Symbol', 'case_id', 'Consequence',
									 'Project', 'Context', 'Vulcan_Local',
									 'Vulcan_Pancancer', 'gscore'
									]
							)
	tcga_data = tcga_data[tcga_data['gscore'] >= 0.6]
	
	# Mark alterations druggable by PanDrugs both by CT drugs and APP drugs
	
	pandrugs_present_CT  = tcga_data['Hugo_Symbol'].isin(pandrugs_db[CT_drug]['gene_symbol'])
	pandrugs_present_FDA = tcga_data['Hugo_Symbol'].isin(pandrugs_db[FDA_drug]['gene_symbol'])

	gain_or_unknown  = tcga_data['Consequence'] != 'LoF'


	# Done as Pineiro et. al. 2018
	tcga_data['PanDrugs_CT']  = np.logical_and(pandrugs_present_CT, gain_or_unknown)
	tcga_data['PanDrugs_APP'] = np.logical_and(pandrugs_present_FDA, gain_or_unknown)

	del pandrugs_db # not needed anymore

	# Build a frequency table
	pandrugs_druggable_CT         = tcga_data['PanDrugs_CT']
	pandrugs_druggable_APP        = tcga_data['PanDrugs_APP']
	
	vulcan_druggable_matched      = tcga_data['Vulcan_Local']
	vulcan_druggable_pancancer    = tcga_data['Vulcan_Pancancer'] 

	cases_all = tcga_data.groupby('Project')['case_id'].nunique().reset_index()
	
	cases_pandrugs_CT 		= tcga_data[pandrugs_druggable_CT].groupby('Project')['case_id'].nunique()
	cases_pandrugs_APP      = tcga_data[pandrugs_druggable_APP].groupby('Project')['case_id'].nunique()
	cases_vulcan_matched 	= tcga_data[vulcan_druggable_matched].groupby('Project')['case_id'].nunique()
	cases_vulcan_pancancer  = tcga_data[vulcan_druggable_pancancer].groupby('Project')['case_id'].nunique()

	cases_all.rename({'case_id': 'total'}, axis='columns', inplace=True)

	cases_all['PanDrugs_CT']  = cases_all['Project'].map(cases_pandrugs_CT).fillna(0)
	cases_all['PanDrugs_APP'] = cases_all['Project'].map(cases_pandrugs_APP).fillna(0)
	cases_all['Vulcan_m'] 	  = cases_all['Project'].map(cases_vulcan_matched).fillna(0)
	cases_all['Vulcan_p']     = cases_all['Project'].map(cases_vulcan_pancancer).fillna(0)

	cases_all['PanDrugs_CT']  = cases_all['PanDrugs_CT'] / cases_all['total'] * 100
	cases_all['PanDrugs_APP'] = cases_all['PanDrugs_APP'] / cases_all['total'] * 100

	cases_all['Vulcan_m'] = cases_all['Vulcan_m'] / cases_all['total'] * 100
	cases_all['Vulcan_p'] = cases_all['Vulcan_p'] / cases_all['total'] * 100


	cases = cases_all.melt(id_vars='Project',
						   var_name='Drug prescription',
						   value_name='count',
						   value_vars=['PanDrugs_CT', 'PanDrugs_APP', 'Vulcan_m', 'Vulcan_p']
						  )

	cases['count'] = cases['count'].round(1)
	
	# Create data for overplotting
	cases_app   = cases[cases['Drug prescription'] == 'PanDrugs_APP']
	cases_ct	= cases[cases['Drug prescription'] == 'PanDrugs_CT']
	cases_vm    = cases[cases['Drug prescription'] == 'Vulcan_m']
	cases_vp    = cases[cases['Drug prescription'] == 'Vulcan_p']

	# Create a new figure
	fig, ax = plt.subplots(figsize=(7,13))
	sns.set_style('whitegrid')

	sns.set_color_codes('pastel')
	figure_pancancer = sns.barplot(x='count',
								   y='Project',
								   data=cases_vp,
								   label='VulcanSpot: Pan-cancer',
								   color='b'
								   )

	sns.set_color_codes('muted')
	figure_local = sns.barplot(
						 x='count',
						 y='Project',
						 data=cases_vm,
						 label='VulcanSpot: Matched tissue',
						 color='b'
						 )

	sns.set_color_codes('pastel')
	figure_pandrugs_ct = sns.barplot(x='count',
									 y='Project',
									 data=cases_ct,
									 label='PanDrugs: CT drugs',
									 color='g')

	sns.set_color_codes('muted')
	figure_pandrugs_app = sns.barplot(x='count',
								 	  y='Project',
								 	  data=cases_app,
								  	  label='PanDrugs: approved drugs',
								  	  color='g'
								  	 )
	


	# Figure details
	plt.title('Cases druggable by PanDrugs and VulcanSpot')

	# Re-order artists and handles so that matched tissue legend bar comes 
	# before pancancer.
	handles, labels = ax.get_legend_handles_labels()
	handles = [handles[3], handles[2], handles[1], handles[0]]
	labels  = [labels[3], labels[2], labels[1], labels[0]]

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

if __name__ == "__main__":
	main()