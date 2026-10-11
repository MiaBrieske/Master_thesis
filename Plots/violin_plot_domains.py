# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 10:06:40 2026

@author: miabr
"""

__usage__ = """ python3 violin_plot.py
					--mean_dir <PATH TO DIR WITH MEAN FILES>
					--out_dir <OUTPUT DIRECTORY>
					OPTIONAL:
					--cropped <DEFINE LOWER X-AXIS VALUE TO CROPP PLOT; DEFAULT: -1>
					"""

### IMPORTS ###
import matplotlib.pyplot as plt
import seaborn as sns
import os, sys
import pandas as pd

def mean_median_list(mean_dir):
	organism_list = []
	organism_list_2 = []
	mean_list = []
	mean_list_different_TF_TG_only = []
	median_list = []
	mean_abs_list = []
	median_abs_list = []
	for mean_files in os.listdir(mean_dir):
		organism = mean_files.split('_')[0]
		mean_file = os.path.join(mean_dir, mean_files)
		with open(mean_file, 'r') as f:
			
			f.readline()						#header
			line = f.readline()
			while line:
				line = line.strip().split('\t')
				if len(line) < 4:
					print(f"Skipping malformed line (fewer than 4 columns): {line}")
					line = f.readline()
					continue
				TF = line[0].split('-')[0]
				TG = line[0].split('-')[1]
				mean = float(line[1])
				median = float(line[2])
				m = int(line[4])
				mean_abs = float(line[5])
				median_abs = float(line[6])
				if m > 1:
					organism_list.append(organism)
					mean_list.append(mean)
					median_list.append(median)
					mean_abs_list.append(mean_abs)
					median_abs_list.append(median_abs)
					if TF != TG:
						mean_list_different_TF_TG_only.append(mean)
						organism_list_2.append(organism)
				
				line = f.readline()
	return organism_list, organism_list_2, mean_list, median_list, mean_abs_list, median_abs_list, mean_list_different_TF_TG_only

def violin_plot(df, output_dir, x, x_col, y_col):
	pairs = df[x_col].drop_duplicates().tolist()
	
	stats_df = df.groupby(x_col)[y_col].agg(n='count').reindex(pairs)
	
	plt.figure(figsize=(14, 7))
	
	color_map = {
		"plants" : "lightgreen",
		"bacteria" : "lightskyblue",
		"fungi" : "plum",
		"animals" : "orange"}
	
	ax = sns.violinplot(
		data = df,
		x=x_col,
		y=y_col,
		palette=color_map,
		order=pairs,
		inner='box',
		cut=0
		)
	if '_abs' in y_col:
		plt.ylim(0, 1)
	else:
		plt.ylim(x, 1)
	
	# n-Werte über jeder Violine einblenden
	for x_pos, pair in enumerate(pairs):
		n = stats_df.loc[pair, 'n']
		ax.text(
			x_pos, 1.02, 	# x = Position der Violine, y = knapp über ylim
			f"n={int(n)}",
			ha='center', va='bottom',
			fontsize=18, rotation=0
			)

	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)
	plt.xlabel(x_col, fontsize=18)
	plt.ylabel(f'{y_col} correlation', fontsize=18)
	plt.title(f'Spearman correlation comparison of domains ({y_col})', fontsize=18, pad=60)
	plt.tight_layout()
	
	
	if output_dir is not None:
		output_file = f"{output_dir}/violinplot_comparison_{y_col}.png"
		plt.savefig(output_file, dpi=300, bbox_inches="tight")
		print(f"Violin plot saved as: {output_file}")
	
	means = df.groupby(x_col)[y_col].mean().reindex(pairs)
	print(f"Means for {y_col}:")
	print(means)
	
	medians = df.groupby(x_col)[y_col].median().reindex(pairs)
	print(f"Medians for {y_col}:")
	print(medians)
	
	plt.show()
	plt.close()
	
def violin_plot_combined(df1, df2, output_dir, x, x_col, y_col, labels):
	color_map = {
		"plants" : "lightgreen",
		"bacteria" : "lightskyblue",
		"fungi" : "plum",
		"animals" : "orange"}
	
	fig, axes = plt.subplots(1,2, figsize=(14, 7), sharey=True)
	
	for ax, df, label in zip(axes, (df1, df2), labels):
		pairs = df[x_col].drop_duplicates().tolist()
		
		stats_df = df.groupby(x_col)[y_col].agg(n='count').reindex(pairs)
		
		sns.violinplot(
			data = df,
			x=x_col,
			y=y_col,
			hue=x_col,
			palette=color_map,
			order=pairs,
			inner='box',
			cut=0,
			ax=ax
				)
		if '_abs' in y_col:
			plt.ylim(0, 1)
		else:
			plt.ylim(x, 1)
	
		# n-Werte über jeder Violine einblenden
		for x_pos, pair in enumerate(pairs):
			n = stats_df.loc[pair, 'n']
			ax.text(
				x_pos, 1.02, 	# x = Position der Violine, y = knapp über ylim
				f"n={int(n)}",
				ha='center', va='bottom',
				fontsize=18, rotation=0
				)
		
		ax.tick_params(axis='x', labelsize=18)
		ax.tick_params(axis='y', labelsize=18)
		ax.set_xlabel(x_col, fontsize=18)
		ax.set_ylabel(f'{y_col} correlation', fontsize=18)
		ax.set_title(label, fontsize=18, pad=30)
		
		means = df.groupby(x_col)[y_col].mean().reindex(pairs)
		print(f"Means for {y_col}:")
		print(means)
	
	fig.suptitle(f"Spearman correlation comparison of domains ({y_col})", fontsize=12)
	fig.tight_layout()
	
	if output_dir is not None:
		output_file = f"{output_dir}/violinplot_comparison_{y_col}_combined.png"
		plt.savefig(output_file, dpi=300, bbox_inches="tight")
		print(f"Violin plot saved as: {output_file}")
	
	means = df.groupby(x_col)[y_col].mean().reindex(pairs)
	print(f"Means for {y_col}:")
	print(means)
	
	medians = df.groupby(x_col)[y_col].median().reindex(pairs)
	print(f"Medians for {y_col}:")
	print(medians)
	
	plt.show()
	plt.close()
	
def main(arguments):

	mean_dir = arguments[arguments.index('--mean_dir') + 1]
	out_dir = arguments[arguments.index('--out_dir') + 1]
	
	os.makedirs(out_dir, exist_ok=True)
	
	if '--cropped' in arguments:
		x = float(arguments[arguments.index('--cropped') + 1])
	else:
		x = -1
		
	organism_list, organism_list_2, mean_list, median_list, mean_abs_list, median_abs_list, mean_list_different_TF_TG_only = mean_median_list(mean_dir)
	
	df_mean = pd.DataFrame()
	df_median = pd.DataFrame()
	df_mean_abs = pd.DataFrame()
	df_median_abs = pd.DataFrame()
	df_mean_wsc = pd.DataFrame()
	
	df_mean = pd.DataFrame({
		'Domain': organism_list,
		'Mean': mean_list
	})
	df_median = pd.DataFrame({
		'Domain': organism_list,
		'Median': median_list
	})
	df_mean_abs = pd.DataFrame({
		'Domain': organism_list,
		'Mean_abs': mean_abs_list
	})
	df_median_abs = pd.DataFrame({
		'Domain': organism_list,
		'Median_abs': median_abs_list
	})
	df_mean_wsc = pd.DataFrame({
		'Domain': organism_list_2,
		'Mean_without_self_regulation': mean_list_different_TF_TG_only
	})

	
	violin_plot(df_mean, out_dir, x, x_col='Domain', y_col='Mean')
	violin_plot(df_median, out_dir, x, x_col='Domain', y_col='Median')
	violin_plot(df_mean_abs, out_dir, x, x_col='Domain', y_col='Mean_abs')
	violin_plot(df_median_abs, out_dir, x, x_col='Domain', y_col='Median_abs')
	violin_plot(df_mean_wsc, out_dir, x, x_col='Domain', y_col='Mean_without_self_regulation')
	
	df_mean_wsc = df_mean_wsc.rename(columns={'Mean_without_self_regulation': 'Mean'})
	
	violin_plot_combined(df_mean, df_mean_wsc, out_dir, x, x_col='Domain', y_col='Mean', labels=('Mean', 'Mean without self regulation'))

if '--mean_dir' in sys.argv and '--out_dir' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )