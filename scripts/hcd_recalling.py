
import requests
import matplotlib.pyplot as plt
import pandas  as pd
import seaborn as sns


def main():

	driver_list = pd.read_csv(
							  '~/Documents/TFM/RAW/Databases/driver_list.tsv',
							  sep='\t',
							  usecols=['Gene', 'Cancer']
							 )
	

	print('Found', driver_list['Gene'].nunique(), 'HCDs')
	vulcan_genes = get_vulcan_genes()

	plot_vulcan_drivers(driver_list, vulcan_genes)


def get_vulcan_genes() -> list:

	raw_list = requests.get('http://vulcanspot.org/api/genes?class=A').json()

	vulcan_gene_list = list()

	for genes in raw_list['data']:
		vulcan_gene_list.append(genes['key'])
	
	return vulcan_gene_list


def plot_vulcan_drivers(driver_list: pd.DataFrame, vulcan_list: list):

	unique_genes = driver_list['Gene'].unique()
	unique_genes = pd.DataFrame(unique_genes, columns=['Gene'])

	unique_genes['In Vulcan'] = unique_genes.isin(vulcan_list)
	print('Found',
	 	  unique_genes[unique_genes['In Vulcan'] == True]['Gene'].nunique(),
	      'in VulcanSpot'
	  	 )

	plt.figure(figsize=(25,15))
	fig = sns.catplot(x='In Vulcan', kind='count', data=unique_genes)
	fig.set(xlabel='H.C.D. in VulcanSpot', ylabel='')
	plt.savefig('/home/sagarcia/Desktop/test.png')


if __name__ == "__main__":
	main()