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
	-a: defines path to alignment file with txids (%%%). Works only on fasta type of alignments.\tREQUIRED\n\
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

def read_align_txt(alignment_path):
	'Reads the txt file'''
	with open(alignment_path, 'r') as myfile:
		data = myfile.read()
	return data

def parse_align(data):
	''''
	Parses through the text file and cuts it into a dictionary with
	keys the taxids and values the accession, name and sequence
	'''
	data = data.replace('OS=','[')
	data = re.sub(' OX=.+\n',']\n',data)
	data = data.replace('\n', '')
	data_list = re.split('%%%|>|\]|\[',data)
	i=0
	data_list=list(filter(None, data_list))
	fasta_data={}
	while i < len(data_list):
		try:
			fasta_data[data_list[i]]=[data_list[i+1],data_list[i+2],data_list[i+3]]
		except IndexError:
			print("Splitting failed on",aln_path,data_list)
			sys.exit(2)
		i+=4
	return fasta_data

def fix_taxids(fastas_dict):
	'Fix old, incorrect or unwanted taxids'
	if '1148' in fastas_dict.keys():
		fastas_dict['1080228'] = fastas_dict.pop('1148')
	if '269483' in fastas_dict.keys():
		fastas_dict['482957'] = fastas_dict.pop('269483')
	if '83333' in fastas_dict.keys():
		fastas_dict['511145'] = fastas_dict.pop('83333')
	if '45157' in fastas_dict.keys():
		fastas_dict['280699'] = fastas_dict.pop('45157')
	if '1936271' in fastas_dict.keys():
		fastas_dict['1841599'] = fastas_dict.pop('1936271')
	if '999953' in fastas_dict.keys():
		fastas_dict['185431'] = fastas_dict.pop('999953')
	if '5660' in fastas_dict.keys():
		fastas_dict['420245'] = fastas_dict.pop('5660')
	if '35128' in fastas_dict.keys():
		fastas_dict['296543'] = fastas_dict.pop('35128')
	if '44689' in fastas_dict.keys():
		fastas_dict['352472'] = fastas_dict.pop('44689')
	return fastas_dict

def fix_multispecie(query_term, orig_fasta_seq, write_list):
	'''
	In the case of entries with identical sequence but accession leading to a MULTISPECIE entry in NCBI.
	'''
	fetch_fastas = subprocess.check_output("esearch -db protein -query '"+query_term+"'| efetch -format fasta", shell=True)
	time.sleep(1)
	fetch_fastas=fetch_fastas.decode('utf-8')
	fasta_sequences = SeqIO.parse(StringIO(fetch_fastas), 'fasta')
	uniq_seq_acc = {}
	for record in fasta_sequences:
		if "MULTISPECIES" in record.description or "pdb|" in record.description:
			pass
		else:
			if orig_fasta_seq == record.seq:
				uniq_seq_acc[orig_fasta_seq]=record.id
	if len(uniq_seq_acc) == 0:
		write_list[4] = re.sub(r'MULTISPECIES: ', '', write_list[4])
		write_csv(write_list)
	else:
		if '|' in uniq_seq_acc[orig_fasta_seq]:
			write_list[1] = uniq_seq_acc[orig_fasta_seq].split("|")[1]
			write_list[2] = 'UNI'
			write_csv(write_list)
			#print (uniq_seq_acc[orig_fasta_seq].split("|")[1], orig_fasta_seq)
		else:
			write_list[1] = uniq_seq_acc[orig_fasta_seq]
			write_csv(write_list)
			#print (uniq_seq_acc[orig_fasta_seq], orig_fasta_seq)
	return True

def write_csv(list_to_write):
	with open(output_path, mode='a') as output_file:
		writer = csv.writer(output_file, delimiter=',')
		writer.writerow(list_to_write)

def output_csv(fastas_dict):
	for i in fastas_dict:
		alnName = ntpath.basename(aln_path).replace('.txt','')
		prot_desc = ' '.join(fastas_dict[i][0].split(" ")[1:]).rstrip()
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
	data = read_align_txt(aln_path)
	fastas_dict = parse_align(data)
	fixedtxid_fastas=fix_taxids(fastas_dict)
	output_csv(fixedtxid_fastas)
	

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