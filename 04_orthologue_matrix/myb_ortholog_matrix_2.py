# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 14:33:58 2026

@author: miabr
"""

"""
This script (def create_heatmap) was developed with the help of AI tools (Claude).
All sections of code were carefully reviewed, tested, and scientifically validated.
"""

__version__ = """v0.1	"""

__usage__ = """
					python3 myb_ortholog_matrix.py (""" +  __version__	+ """)
					--in <DIRECTORY_WITH_MYB_RESULTS>
					--ref <MYB_REFERENCE_FILE>
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

def structure_out_file(myb_dir):
	structure_list = ['RefMYB', 'Count/Total', 'Percentage [%]']
	for species in os.listdir(myb_dir):
		complete_path = os.path.join(myb_dir, species)
		if os.path.isdir(complete_path):
			name = species.split('_')
			structure_list.append('_'.join(name[1:]))
		else:
			continue
	
	return structure_list
	
def refmybs_outfile(refmyb_file, outfile, structure_list):
	
	num_file = 1
	out = f"{outfile}_{num_file}.txt"
	
	with open (refmyb_file, 'r') as f:
		
		with open (out, 'w') as f_out:
			
			f_out.write('\t'.join(structure_list) + '\n')
			line = f.readline()
			while line:
				line = line.strip().split('\t')
				if len(line) < 2:
					print(line)
					line = f.readline()
					continue
				lines = [line[1], '0/0', '0']
				f_out.write('\t'.join(lines) + '\n')
				
				line = f.readline()
	return out, num_file

def plant_data(out, outfile, num_file, myb_dir):
	
	total = 0
	
	count_dict = {}
	with open(out, 'r') as f:
		
		header = f.readline()
		if header:
			line = f.readline()
			header = False
		while line:
			refmyb = line.strip().split('\t')[0]
			count_dict[refmyb] = 0  # Initialwert 0
			line = f.readline()
			
	for species in os.listdir(myb_dir):
		complete_path = os.path.join(myb_dir, species, "RESULTS")
		if not os.path.isdir(complete_path):
			print (f"Skip {species}: no RESULTS file")
			continue 
		else: 
			species = species.split('_')
			species = '_'.join(species[1:])
			print (f"Processing {species}")
		
		for file in os.listdir(complete_path):
			if file.startswith('03b_'):
				file_path = os.path.join(complete_path, file)
				
				if os.path.exists(file_path):
					num_file += 1
					total += 1
					old_out = out
					new_out = f"{outfile}_{num_file}.txt"
					
					with open(file_path, 'r') as f:
						
						species_hits = {}
						line = f.readline()
						while line:
							myb = line.strip().split('\t')[0]
							ref = line.strip().split('\t')[1]
							if ref not in species_hits:
								species_hits[ref] = []
							species_hits[ref].append(myb)
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
										myb_list = species_hits[ref]
										count_dict [ref] += 1
										myb = ','.join(myb_list)
									else: 
										myb = '-'
											
									if pos == len(line):
										line.append(myb)
									elif pos < len(line):
										line[pos] = myb
									else:
										line.extend([""] * (pos - len(line)))
										line.append(myb)
									
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
	
	heatmap_data.index = df['RefMYB']
	
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
	plt.ylabel("RefMYB")
	plt.title("MYB Ortholog Presence/Absence")
	
	plt.tight_layout()
	
	out = f"{outfile}_heatmap.jpg"
	plt.savefig(out, dpi=300)
	plt.close()
	
	print(f"Heatmap saved {out}")
	
def main(arguments):
	
	ref_file = arguments[arguments.index('--ref') + 1]
	myb_dir = arguments[arguments.index('--in') + 1]
	output_prefix = arguments[arguments.index('--out') + 1]
	
	out_dir = output_prefix.rsplit('/', 1)[0]
	os.makedirs(out_dir, exist_ok=True)
	
	print(f"Creating structure of {output_prefix}-file.")
	structure_list = structure_out_file(myb_dir)
	print(f"Adding RefMYBs to {output_prefix}-file.")
	out, num_file = refmybs_outfile(ref_file, output_prefix, structure_list)
	print(f"Adding plant data to {output_prefix}-file.")
	result_file = plant_data(out, output_prefix, num_file, myb_dir)
	
	###  optional heatmap  ###
	if '--heatmap' in sys.argv:
		create_heatmap(result_file, output_prefix)
		
if '--in' in sys.argv and '--out' in sys.argv and '--ref' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )