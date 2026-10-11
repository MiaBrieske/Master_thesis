# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 13:49:37 2026

@author: miabr
"""

__usage__ = """ python3 transcript_to_gene_final.py
					--TF_TG_list <LIST WITH HEADER Transcription factor	Transcript ID	Target gene -> Transcript ID COLUMN EMPTY>
					--gene_from_header <OUT FILE FROM gene_from_header_cds_final.py>
					--out <OUTPUT FILE>
					"""

import sys, os

def main(arguments):
	
	TF_TG_list = arguments[arguments.index('--TF_TG_list') + 1]
	transcript_to_gene = arguments[arguments.index('--gene_from_header') + 1]
	out = arguments[arguments.index('--out') + 1]
	
	out_dir = out.rsplit('/', 1)[0]
	os.makedirs(out_dir, exist_ok=True)
	
	with open(transcript_to_gene, 'r') as f:
		
		transcript_to_gene_dict = {}				#transcript_to_gene_dict = {gene : ID}
		line = f.readline()
		while line:
			line = line.strip().split('\t')
			ID = line[0]
			gene = line[1]
			transcript_to_gene_dict[gene] = ID
			
			line = f.readline()
	
	with open(TF_TG_list, 'r') as f:
		with open(out, 'w') as o:
			
			header = f.readline()
			header = header.strip()
			o.write(f"{header}\n")
			line = f.readline()
			while line:
				line = line.strip().split('\t')
				if len(line) >= 3:
					TF = line[0]
					gene = line[2]
				else:
					print(line[0])
				for key, value in transcript_to_gene_dict.items():
					if key == gene:
						o.write(f"{TF}\t{value}\t{gene}\n")
				line = f.readline()

if '--TF_TG_list' in sys.argv and '--gene_from_header' in sys.argv and '--out' in sys.argv:
	main( sys.argv )
else:
	sys.exit( __usage__ )