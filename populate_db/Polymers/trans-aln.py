#!/usr/bin/env python3
import re
import sys
import csv
import getopt
from Bio import SeqIO
from Bio import Entrez
import xml.dom.minidom as minidom
Entrez.email = 'ap211@gatech.edu'

#Run ./transl-aln.py alignment_file_path translation_table_path
#Add option to output number of fasta entries and number of processed entries
#Add option to switch between downloading and using downloaded translational table

def usage():
	print (\
	"USAGE:\n./transl-aln.py -a [alignment_file_path] -t [translation_table_path] -n -h\n\
	-a: defines path to alignment file. Works only on fasta type of alignments.\tREQUIRED\n\
	-t: defines path to translation table or txid list of your species.\t\tREQUIRED\n\
	-d: Using this and a txid list would connect to ENTREZ and download aliases from there.\n\
	-n: outputs numbers of fasta entries and number of processed entries\n\
	-h: prints this\
")

try:
	opts, args = getopt.getopt(sys.argv[1:], 'a:t:d;n;h', ['alignment=', 'translation=', 'download_aliases', 'number_missing', 'help'])
except getopt.GetoptError:
	usage()
	sys.exit(2)

num_miss=0
download_aliases=None
for opt, arg in opts:
	if opt in ('-h', '--help'):
		usage()
		sys.exit(2)
	elif opt in ('-a', '--alignment'):
		aln_path = arg
	elif opt in ('-t', '--translation'):
		transl_table_path = arg
	elif opt in ('-d', '--download_aliases'):
		download_aliases = 1
	elif opt in ('-n', '--number_missing'):
		num_miss=1
	else:
		usage()
		sys.exit(2)
		
try:
	aln_path
except NameError:
	usage()
	sys.exit(2)
try:
	transl_table_path
except NameError:
	usage()
	sys.exit(2)

#Parses through a translation table or txid list
def parse_transl_table(transl_table_path,alias_dict):
	with open(transl_table_path, 'r') as f:
		transl_table = f.read().splitlines()
	for line in transl_table:
		elements=line.split(',')
		if (";" in line):
			aliases=elements[1].split(';')
			alias_dict[elements[0]]=aliases
		else:
			alias_dict[elements[0]]=""
	return alias_dict

#Accesses Entrez and creates translation dictionary
def gen_al_dict(alias_dict):
	print("generating")
	for txid,other_names in alias_dict.items():
		print(txid)
		handle = Entrez.efetch(db="taxonomy", id=txid)
		entr_xml=minidom.parseString(handle.read())
		Sname=entr_xml.getElementsByTagName('ScientificName')
		other_names=list()
		parse_othernames(entr_xml, other_names)
		other_names.append(Sname[0].firstChild.nodeValue)
		alias_dict[txid]=other_names
	#Used for first time translation table generation
	for txid,names in alias_dict.items():
		print(txid,';'.join(names),sep=',')
	return alias_dict

#Further parsing of Entrez xml info for more aliases
#Other other names: EquivalentName,GenbankSynonym,DispName
def parse_othernames(xml_object,other_names):
	Onames=xml_object.getElementsByTagName('OtherNames')
	for node in Onames:
		synonyms=node.getElementsByTagName('Synonym')
		for x in synonyms:
			other_names.append(x.firstChild.nodeValue)
		name=node.getElementsByTagName('Name')
		for x in name:
			dispname=x.getElementsByTagName('DispName')
			for y in dispname:
				other_names.append(y.firstChild.nodeValue)
		eq_name = node.getElementsByTagName('EquivalentName')
		for x in eq_name:
			other_names.append(x.firstChild.nodeValue)
		gb_name = node.getElementsByTagName('GenbankSynonym')
		for x in gb_name:
			other_names.append(x.firstChild.nodeValue)
		return other_names


def translate(file_path,alias_dict,out_dict,seq_number):
	for record in SeqIO.parse(file_path, "fasta"):
		gi=record.id.split("_", 1)[0]
		gi = gi.replace("NR", "NR_")
		record.id = record.id.split("_", 1)[1]
		name_split=record.id.split("_",1)
		namettpf = record.id.split("_")
		name = re.sub("part_","",name_split[1])
		name = re.sub("_"," ",name)
		if re.search('OS=', name):
			name = re.sub(' OX=.+',']',name)
			name = name.replace('OS=','[')
		else:
			name = re.sub(' OX=.+','',name)
		if (num_miss == 1):
			print (">",record.id,sep='')
		check=0
		seq_number+=1
		if 'TT' == namettpf[1]:
			name = ">"+record.id+"_262724"
			out_dict[name]=record.seq
			print(">",namettpf[0],"_262724_",namettpf[1],"\n",record.seq,sep='')
			check=1
		elif 'PF' == namettpf[1]:
			name = ">"+record.id+"_186497"
			out_dict[name]=record.seq
			print(">",namettpf[0],"_186497_",namettpf[1],"\n",record.seq,sep='')
			check = 1
		for txid,aliases in alias_dict.items():
			if check == 1:
				break
			for alias in aliases:
				if re.search(name, alias):
					name = ">"+txid+"_"+record.id
					out_dict[name]=record.seq
					print (">",namettpf[0],"_",txid,"_",aliases[0].replace(" ","-"),"|",gi,"\n",record.seq,sep='')
					break
	return (out_dict,seq_number)

###EXECUTE###
alias_dict=dict()
out_dict=dict()
seq_number=0

parse_transl_table(transl_table_path,alias_dict)
#Checks if we are downloading or using provided aliases
if (download_aliases == 1):
	gen_al_dict(alias_dict)
	something,seqd=translate(aln_path,alias_dict,out_dict,seq_number)
	if (num_miss == 1):
		print (len(something),seqd)
else:
	something,seqd=translate(aln_path,alias_dict,out_dict,seq_number)
	if (num_miss == 1):
		print (len(something),seqd)
