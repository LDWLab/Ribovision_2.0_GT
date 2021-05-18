#!/usr/bin/env python3
import re, sys, Bio.Align, getopt, subprocess, csv
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from io import StringIO
import time, ntpath

def usage():
	print (\
	"USAGE:\n./parse_accession.py -a [alignment_file_path]-h\n\
	-a: defines path to alignment file.\tREQUIRED\n\
	-o: defines path to output csv file \t\t\t\t\t\t\t\tREQUIRED\n\
	-h: prints this\
")

try:
	opts, args = getopt.getopt(sys.argv[1:], 'a:o:h', ['alignment=', 'help'])
except getopt.GetoptError:
	usage()
	sys.exit(2)

for opt, arg in opts:
	if opt in ('-h', '--help'):
		usage()
		sys.exit(2)
	elif opt in ('-a', '--alignment'):
		aln_path = arg
	elif opt in ('-o', '--output'):
		output_path = arg
	else:
		usage()
		sys.exit(2)

def parse_align(seqObj):
	''''
	Parses through the sequence object and generates a dictionary with
	keys the taxids and values the accession, name, sequence, and should add description for now it duplicates with name.
	Cuts off the last letter from the name to remove uL02e etc.
	'''
	outDict = dict()
	for seq in seqObj:
		name = '|'.join(seq.id.split('|')[1:])
		outDict[seq.id.split('_')[1]] = [name, seq.id.split('_')[0][:4], str(seq.seq), seq.id.split('_')[0][:4]]
	return outDict

def write_csv(list_to_write):
	with open(output_path, mode='a') as output_file:
		writer = csv.writer(output_file, delimiter=',')
		writer.writerow(list_to_write)

def output_csv(fastas_dict):
	for i in fastas_dict:
		alnName = ntpath.basename(aln_path).replace('.txt','').replace('e_new.fa','').replace('_new.fas','')
		prot_desc = f'Ribosomal protein {alnName}'
		if 'RecName: Full=' in prot_desc:
			prot_desc = re.sub(r'RecName: Full=', '', prot_desc)
		if '; AltName:' in prot_desc:
			prot_desc = re.sub(r'; AltName:.* ', '', prot_desc)
		if 'PREDICTED: ' in prot_desc:
			prot_desc = re.sub(r'PREDICTED: ', '', prot_desc)
		if '|' in fastas_dict[i][0]:
			write_csv([i,fastas_dict[i][0].split(" ")[0].split("|")[1], 'UNI',alnName,prot_desc,fastas_dict[i][2]])
			#print (i,fastas_dict[i][0].split(" ")[0].split("|")[1], 'UNI',alnName,fastas_dict[i][2], sep=',')
			pass
		elif 'MULTISPECIES:' in fastas_dict[i][0]:
			prot_desc = re.sub(r'MULTISPECIES: ', '', prot_desc)
			query_term = 'txid'+i+'[Organism] AND'+fastas_dict[i][0].split(":")[1]
			write_list = [i,fastas_dict[i][0].split(" ")[0], 'NCBI',alnName,prot_desc,fastas_dict[i][2]]
			write_csv(write_list)
			#fix_multispecie(query_term, fastas_dict[i][2],write_list)
		elif '(nucleomorph)' in prot_desc:
			print(i,fastas_dict[i][0].split(" ")[0],alnName,prot_desc,fastas_dict[i][2])
		elif '(macronuclear)' in prot_desc:
			print(i,fastas_dict[i][0].split(" ")[0],alnName,prot_desc,fastas_dict[i][2])
		else:
			write_csv([i,fastas_dict[i][0].split(" ")[0], 'NCBI',alnName,prot_desc,fastas_dict[i][2]])
			#print (i,fastas_dict[i][0].split(" ")[0], 'NCBI',alnName,fastas_dict[i][2], sep=',')
			pass
	return True

def main():
	sequences = SeqIO.parse(aln_path, "fasta")
	fastas_dict = parse_align(sequences)

	output_csv(fastas_dict)
	

if __name__ == "__main__":
	main()

#def check_GI(accession):
#	'''
#	Move out of here.
#	'''
#	print(accession)
#	geneinfo_output = subprocess.check_output("elink -db protein -id '"+accession+"' -target gene | efetch -format uid", shell=True)
#	time.sleep(1)
#	specieinfo_output = subprocess.check_output("elink -db protein -id '"+accession+"' -target taxonomy | efetch -format uid", shell=True)
#	time.sleep(1)
#	if geneinfo_output is None:
#		pass
#	else:
#		print (specieinfo_output,geneinfo_output,end=' ')