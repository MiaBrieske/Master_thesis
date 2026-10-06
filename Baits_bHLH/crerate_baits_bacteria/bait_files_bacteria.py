# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 13:46:57 2026

@author: miabr
"""

__usage__ = """ python3 bait_files_bacteria.py
					--pep_file <PEP FILE FOR BAITS>
					--cds_file <CDS FILE FOR BAITS>
					--translation_file <TRANSLATION FILE>
					--out_path <PATH FOR RESULT FILES>
					--out_prefix <PREFIX FOR BAIT FILES>
					"""

import os, sys

def main(arguments):
	
	pep_file = arguments[arguments.index('--pep_file') + 1]
	cds_file = arguments[arguments.index('--cds_file') + 1]
	TG_list = arguments[arguments.index('--translation_file') + 1]
	out = arguments[arguments.index('--out_path') + 1]
	prefix = arguments[arguments.index('--out_prefix') + 1]
	
	os.makedirs(out, exist_ok=True)
	out = os.path.join(out, f"TFs_{prefix}.tsv")
	fasta = os.path.join(out, f"TFs_{prefix}.fasta")
	baits = os.path.join(out, f"{prefix}TFs_baits.fasta")
	
	with open(TG_list) as f:
		f.readline()
		TG_gene = {}
		TG_locus = {}
		line = f.readline()
		while line:
			line = line.strip().split('\t')
			if len(line) < 4:
				print(line)
			ID = line[0]
			gene = line[2]
			locus = line[3]
			if '|' in gene:
				gene = gene.split('|')
			else:
				gene = [gene]
			if '|' in locus:
				locus = locus.split('|')
			else:
				locus = [locus]
			for element in gene:
				TG_gene[ID] = element
			for element in locus:
				TG_locus[ID] = element
			
			line = f.readline()
		print(TG_list)
		print(len(TG_gene))
		print(len(TG_locus))
	
	with open(cds_file, 'r') as f:
			
		x = 0
		TG_gene_found = []
		TG_locus_found = []
		dictionary = {}
		line = f.readline()
		while line:
			line = line.strip()
			if line.startswith('>'):
				line = line.split(' ')
				ID = line[0].split('>')[1]
				if "gene" in line[1]:
					gene = line[1].split(']')
					gene = gene[0].split('=')[1]
					#print(gene.upper())
					locus = line[2].split(']')
					locus = locus[0].split('=')[1]
					
					found = False
					for key, value in TG_gene.items():
						if gene == value:
							dictionary[ID] = key
							x += 1
							TG_gene_found.append(key)
							found = True
							break
					if not found:
						for key, value in TG_locus.items():
							if locus == value:
								dictionary[ID] = key
								x += 1
								TG_locus_found.append(key)
				elif "locus_tag" in line[1]:
					locus = line[1].split(']')
					locus = locus[0].split('=')[1]
					for key, value in TG_locus.items():
						if locus == value:
							dictionary[ID] = key
							x += 1
							TG_locus_found.append(key)
						
			line = f.readline()
	#print(dictionary)
	with open(pep_file, 'r') as f:
		with open(out, 'w') as o:
			
			line = f.readline()
			while line:
				line = line.strip()
				if line.startswith('>'):
					gene = line.split('>')[1]
					print(gene)
					line_key = False
					for key, value in dictionary.items():
						#print(key)
						if key == gene:
							o.write(f">{value}\n")
							line_key = True
					'''
					if not line_key:
						o.write(f">{gene}\n")
					'''
				elif not line.startswith('>') and line_key:
					o.write(f"{line}\n")
				
				line = f.readline()
				
	for key, value in TG_locus.items():
		if key not in TG_locus_found:
			print(f"Locus not found: {key}")
			
	for key, value in TG_gene.items():
		if key not in TG_gene_found:
			print(f"Gene not found: {key}")
	
	with open(out, "r", encoding="utf-8") as f:
		content = f.read()

	with open(fasta, "w", encoding="utf-8") as out_txt:
		out_txt.write(content)
	
	with open(baits, "w", encoding="utf-8") as out_tsv:
		out_tsv.write(content)

if '--pep_file' in sys.argv and '--cds_file' in sys.argv and '--translation_file' in sys.argv and '--out_path' in sys.argv and '--out_prefix' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )
