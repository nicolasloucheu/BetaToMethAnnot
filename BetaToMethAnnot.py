import pandas as pd
import os
from os import path
import pickle


#Choose input file
input_file = "GSE2366698_beta.csv"

# Change name of output folder
name_of_folder = "GSE2366698"




# Create the folder
if not os.path.exists(name_of_folder):
	os.mkdir(name_of_folder)

# List all chromosomes
chrom_lst = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y"]
top_z = pd.DataFrame(columns=[name_of_folder, 'CHR', 'MAPINFO'])
z_for_mean = pd.DataFrame(columns=[name_of_folder, 'CHR', 'MAPINFO'])

# Import beta-values
beta = pd.read_csv(input_file, index_col=0)
beta = beta.rename({beta.columns[0]: name_of_folder}, axis='columns')

# Import manifest
fields = ["IlmnID", "CHR", "MAPINFO"]
dtype_chr = {'CHR': 'str'}
manifest = pd.read_csv("hg38_manifest.csv", usecols=fields, index_col=0, dtype=dtype_chr)

# Merge beta and manifest
beta = pd.merge(manifest, beta, left_index=True, right_index=True)

# Generate a file and an index file for each chromosome 
for i in chrom_lst:
	beta_chrom = beta.loc[beta["CHR"] == i]
	beta_chrom = beta_chrom[[name_of_folder, "CHR", "MAPINFO"]]
	beta_chrom.to_csv(f"{name_of_folder}/{name_of_folder}_chrom_{i}.csv.gz", compression="gzip")
	index_chrom = beta_chrom.MAPINFO.to_list()
	with open(f"{name_of_folder}/{name_of_folder}_chrom_{i}_lst.txt", "wb") as fp:
		pickle.dump(index_chrom, fp)
	beta_chrom = beta_chrom.transpose()
	series_controls = pd.read_csv(f"chrom_means/chrom_{i}_means.csv.gz", compression="gzip", usecols=[0, 2, 3], index_col=0).transpose()
	cols_tokeep = series_controls.columns
	new_beta = beta_chrom[beta_chrom.columns & cols_tokeep].drop(['CHR', 'MAPINFO'], axis='index')
	z_scores = (abs(new_beta-series_controls.loc['MEAN'])/series_controls.loc['STD']).transpose().dropna(axis=0)
	z_scores = pd.merge(manifest, z_scores, left_index=True, right_index=True)
	z_scores = z_scores[[name_of_folder, "CHR", "MAPINFO"]]
	z_for_mean = z_for_mean.append(z_scores)
	z_scores.to_csv(f"{name_of_folder}/{name_of_folder}_z_{i}.csv.gz", compression="gzip")
	index_z = z_scores.MAPINFO.to_list()
	with open(f"{name_of_folder}/{name_of_folder}_index_{i}.txt", "wb") as fp:
		pickle.dump(index_z, fp)
	top_z = top_z.append(z_scores.loc[z_scores[name_of_folder] > 5]).sort_values(by=name_of_folder, axis=0, ascending=False)
	print(top_z)

top_z.to_csv(f"{name_of_folder}/top_z_scores.csv.gz", compression='gzip')
mean_z = z_for_mean[[name_of_folder]].mean().to_list()
with open(f"{name_of_folder}/z_score_mean.pickle", "wb") as fp:
		pickle.dump(mean_z, fp)