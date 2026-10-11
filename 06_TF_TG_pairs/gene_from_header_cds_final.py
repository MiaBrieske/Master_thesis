# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 12:21:03 2026

@author: miabr
"""

__usage__ = """ python3 gene_from_header_cds_final.py
					--TG_list_unique <LIST WITH UNIQUE TG>
					--cds_model <CDS FILE FROM MODEL ORGANSIM BEFORE CLEAN>
					--out_file <PATH TO OUT FILE>
					"""

import sys, os

def main(arguments):
	
	TG_list_unique = arguments[arguments.index('--TG_list_unique') + 1]
	cds_file = arguments[arguments.index('--cds_model') + 1]
	out = arguments[arguments.index('--out_file') + 1]
	
	out_dir = out.rsplit('/', 1)[0]
	os.makedirs(out_dir, exist_ok=True)
	
	with open(TG_list_unique) as f:
		
		TG_gene = []
		TG_locus = []
		line = f.readline()
		while line:
			line = line.strip().split('\t')
			if len(line) < 2:
				gene = line[0]
				locus = None
			else:
				gene = line[1]
				locus = line[0]
				if ';' in locus:
					locus = locus.split(';')
					locus_list = True
				else:
					locus_list = False
			if locus_list:
				for element in locus:
					if gene not in TG_gene or element not in TG_locus:
						TG_gene.append(gene)
						TG_locus.append(element)
			else:
				if gene not in TG_gene or locus not in TG_locus:
					TG_gene.append(gene)
					TG_locus.append(locus)
			
			line = f.readline()
		print(len(TG_gene))
		print(len(TG_locus))
	
	with open(cds_file, 'r') as f:
		with open(out, 'w') as o:
			
			x = 0
			TG_gene_found = []
			TG_locus_found = []
			line = f.readline()
			while line:
				line = line.strip()
				if line.startswith('>'):
					line = line.split(' ')
					if line[0].startswith('>lcl_'):
						line = line[0].replace('_', '|', 1)
					ID = line[0].split('>')[1]
					if "gene" in line[1]:
						gene = line[1].split(']')
						gene = gene[0].split('=')[1]
						
						locus = line[2].split(']')
						locus = locus[0].split('=')[1]
						if ':' in locus:
							locus = locus.split(':')[1]
						
						if gene in TG_gene:
							o.write(f"{ID}\t{gene}\n")
							x += 1
							TG_gene_found.append(gene)
						elif locus in TG_locus:
							o.write(f"{ID}\t{locus}\n")
							x += 1
							TG_locus_found.append(locus)
					elif "locus_tag" in line[1]:
						locus = line[1].split(']')
						locus = locus[0].split('=')[1]
						if locus in TG_gene:
							o.write(f"{ID}\t{locus}\n")
							x += 1
							TG_locus_found.append(locus)
							
				line = f.readline()
	
	print(x)
	
	
	for element in TG_locus:
		if element not in TG_locus_found:
			print(f"Locus not found: {element}")
	'''
	for element in TG_gene:
		if element not in TG_gene_found:
			print(f"Gene not found: {element}")
	'''

if '--TG_list_unique' in sys.argv and '--cds_model' in sys.argv and '--out_file' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )