# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 14:33:58 2026

@author: miabr
"""

__version__ = """v0.1	"""

__usage__ = """
					python3 bHLH_ortholog_matrix.py (""" +  __version__  + """)
					--in <DIRECTORY_WITH_bHLH_RESULTS>
					--ref <bHLH_REFERENCE_FILE>
					--out <OUTPUT_DIRECTORY>
					
					OPTIONAL
					--heatmap <NO_INPUT_NEEDED>
			"""
			
###  Imports  ###
import os, sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

###  Functions  ###

def structure_out_file(bhlh_dir):
	structure_list = ['RefMember', 'Count/Total', 'Percentage [%]']
	for species in os.listdir(bhlh_dir):
		if species == 'log' or species == 'RESULTS' or species == '02_in_out_bHLH_analysis_trees':
			continue
		else:
			complete_path = os.path.join(bhlh_dir, species)
			if os.path.isdir(complete_path):
				name = species.split('.')
				structure_list.append(name[0])
			else:
				continue
	
	return structure_list
	
def refbHLHs_outfile(ref_file, outfile, structure_list):
	
	num_file = 1
	out = f"{outfile}_{num_file}.txt"
	
	with open (ref_file, 'r') as f:
		
		with open (out, 'w') as f_out:
			
			f_out.write('\t'.join(structure_list) + '\n')
			line = f.readline()
			while line:
				line = line.strip().split('\t')
				if len(line) < 3:
					print(line)
					line = f.readline()
					continue
				if line[2] == '':
					print(line)
					print(f"RefMember: {line[0]}")
					lines = [line[0], '0/0', '0']
					f_out.write('\t'.join(lines) + '\n')
				else:
					print(f"RefMember: {line[2]}")
					lines = [line[2], '0/0', '0']
					f_out.write('\t'.join(lines) + '\n')
				
				line = f.readline()
	return out, num_file

def plant_data(out, outfile, num_file, bhlh_dir):
	
	total = 0
	
	count_dict = {}
	with open(out, 'r') as f:
		
		header = f.readline()
		if header:
			line = f.readline()
			header = False
		while line:
			refbhlh = line.strip().split('\t')[0]
			count_dict[refbhlh] = 0  # Initialwert 0
			line = f.readline()
			
	for species in os.listdir(bhlh_dir):
		if species == 'log' or species == 'RESULTS' or species == '02_in_out_bHLH_analysis_tree':
			continue
		else:
			complete_path = os.path.join(bhlh_dir, species, "RESULTS")
			if not os.path.isdir(complete_path):
				print (f"Skip {species}: no RESULTS file")
				continue 
			else: 
				species = species.split('.')
				species = species[0]
				print (f"Processing {species}")
			
			for file in os.listdir(complete_path):
				if file == '04b_candidate_2_reference_mapping_file.txt':
					file_path = os.path.join(complete_path, file)
					
					if os.path.exists(file_path):
						print(file_path)
						num_file += 1
						total += 1
						old_out = out
						new_out = f"{outfile}_{num_file}.txt"
						
						with open(file_path, 'r') as f:
							
							species_hits = {}
							line = f.readline()
							while line:
								bhlh = line.strip().split('\t')[0]
								ref = line.strip().split('\t')[1]
								if ref not in species_hits:
									species_hits[ref] = []
								species_hits[ref].append(bhlh)
								line = f.readline()
								
							with open(old_out, 'r') as oldout:
								with open(new_out, 'w') as newout:
									print(f"Creating {new_out} adding data from {species}.")
									
									header_line = oldout.readline()
									newout.write(header_line)
									pos = header_line.strip().split('\t').index(species)
									
									line = oldout.readline()
									while line:
										line = line.strip().split('\t')
										ref = line[0]
										if ref in species_hits:
											bhlh_list = species_hits[ref]
											count_dict [ref] += 1
											bhlh = ','.join(bhlh_list)
										else: 
											bhlh = '-'
												
										if pos == len(line):
											line.append(bhlh)
										elif pos < len(line):
											line[pos] = bhlh
										else:
											line.extend([""] * (pos - len(line)))
											line.append(bhlh)
										
										line[1] = f"{count_dict[ref]}/{total}"
										line[2] = f"{count_dict[ref]/total*100:.2f}"
										newout.write('\t'.join(line) + '\n')
										
										line = oldout.readline()
										
						out = new_out
						os.remove(old_out)
						print(f"{old_out} removed.")
	###  sorting file  ###
	df = pd.read_csv(out, sep='\t')
	df_sorted = df.sort_values(by='Percentage [%]', ascending=False)
	result_file = f'{outfile}_RESULTS.txt'
	df_sorted.to_csv(result_file, sep='\t', index=False)
	
	print(f"{out} has been sorted and saved as {result_file}.")
	
	return result_file

def create_heatmap(result_file, outfile):
	
	df = pd.read_csv(result_file, sep='\t')
	
	###  only species columns  ###
	species_cols = df.columns[3:]
	
	###  Presence/Absence matrix  ###
	heatmap_data = (df[species_cols] != '-').astype(int)
	
	heatmap_data.index = df['RefMember']
	
	###  NEW: CSV export of the matirx  ###
	csv_out = f"{outfile}_heatmap.csv"
	heatmap_data.to_csv(csv_out)
	print(f"Presence/Absence-Matrix gespeichert als {csv_out}")
	
	plt.figure(figsize=(len(species_cols)*0.6 + 4, len(df)*0.25 + 4))
	
	sns.heatmap(
		heatmap_data,
		cmap=["white", "black"],
		cbar=False,
		linewidths=0.5,
		linecolor="gray"
		)
	
	plt.xlabel("Species")
	plt.ylabel("RefMemeber")
	plt.title("bHLH Ortholog Presence/Absence")
	
	out = f"{outfile}_heatmap.jpg"
	plt.savefig(out, dpi=300, bbox_inches='tight')
	plt.close()
	
	print(f"Heatmap saved {out}")
	
def main(arguments):
	
	ref_file = arguments[arguments.index('--ref') + 1]
	bhlh_dir = arguments[arguments.index('--in') + 1]
	output_prefix = arguments[arguments.index('--out') + 1]
	
	print(f"Creating structure of {output_prefix}-file.")
	structure_list = structure_out_file(bhlh_dir)
	print(f"Adding RefbHLHs to {output_prefix}-file.")
	out, num_file = refbHLHs_outfile(ref_file, output_prefix, structure_list)
	print(f"Adding plant data to {output_prefix}-file.")
	result_file = plant_data(out, output_prefix, num_file, bhlh_dir)
	
	###  optional heatmap  ###
	if '--heatmap' in sys.argv:
		create_heatmap(result_file, output_prefix)
		
if '--in' in sys.argv and '--out' in sys.argv and '--ref' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )