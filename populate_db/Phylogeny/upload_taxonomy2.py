#!/usr/bin/env python3
import sys, re, os, csv, getpass, argparse, mysql.connector, time, json
from Bio import Entrez
import pandas as pd

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Update structural folds tables from ECOD and SCOPe, given a PDB ID')
	parser.add_argument('input_taxids', help='PDB identifier to query', type=str)
	parser.add_argument('user_name', help='Username for connecting to DESIRE', type=str)
	parser.add_argument('-host','--db_host', help='Defines database host (default: 130.207.36.75)', type=str, default='130.207.36.75')
	parser.add_argument('-schema','--db_schema', help='Defines schema to use (default: SEREB2)', type=str, default='SEREB2')
	parser.add_argument('-dl','--download_most_recent_phylogeny', help='Update latest phylogeny.', default=False, action="store_true")
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

def initiate_connection(uname, host, database):
	pw = getpass.getpass("Password: ")
	cnx = mysql.connector.connect(user=uname, password=pw, host=host, database=database)
	return cnx

def load_species_taxids(file_path):
	with open(file_path) as f:
		input_species_txids = f.read().splitlines()
	return input_species_txids

def download_phylogeny(input_taxids):
	Entrez.email = "ppenev@gatech.edu"  # Always tell NCBI who you are
	phylogeny_structure = dict()
	phylogeny_levels = ['superkingdom', 'phylum', 'class', 'order', 'family' , 'genus', 'species']
	phylogeny_levels.reverse()
	for taxid in input_taxids:
		time.sleep(1)
		oneTaxon = dict()
		handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
		records = Entrez.read(handle)
		oneTaxon[records[0]['TaxId']] = ['strain',records[0]['ScientificName']]
		lineages = records[0]["LineageEx"]
		for i in lineages:
			for j in sorted(phylogeny_levels):
				if j == i['Rank']:
					oneTaxon[i['TaxId']] = [i['Rank'],i['ScientificName']]
		phylogeny_structure[records[0]['TaxId']] = oneTaxon
		print(oneTaxon)
	json.dump(phylogeny_structure, open("./data/downloaded_phylogeny.json", 'w' ))
	return phylogeny_structure

def main(commandline_arguments):
	comm_args = create_and_parse_argument_options(commandline_arguments)
	input_taxids = load_species_taxids(comm_args.input_taxids)

	if comm_args.download_most_recent_phylogeny:
		phylogeny_structure = download_phylogeny(input_taxids)
	else:
		phylogeny_structure = json.load( open( "./data/downloaded_phylogeny.json" ) )

	print(phylogeny_structure)

	
	# cnx = initiate_connection(comm_args.user_name, comm_args.db_host, comm_args.db_schema)
	# cursor = cnx.cursor()
	# '''---Species Table---'''################################################################################
	# for strain_Id in input_taxids:
	# 	strain = phylogeny_structure[strain_Id][strain_Id][1]
	# 	query = ("INSERT INTO `Species`(`strain_id`,`strain`) VALUES('"+strain_Id+"','"+strain+"')")
	# 	#print(query)
	# 	cursor.execute(query)

	#Fix the recursive relationship
	# '''---TaxGroups Table---'''##############################################################################
	# empty = []
	# phylogeny_levels.insert(0,'strain')
	# cursor.execute("SET FOREIGN_KEY_CHECKS=0")
	# for strain_Id in input_taxids:
	# 	list_with_levels=['no_strain','no_species','no_genus','no_family','no_order','no_class','no_phylum','no_superkingdom']
	# 	for ids in phylogeny_structure[strain_Id].keys():
	# 		list_with_levels[phylogeny_levels.index(phylogeny_structure[strain_Id][ids][0])] = ids
	# 	for ids in phylogeny_structure[strain_Id].keys():
	# 		group_lvl = phylogeny_structure[strain_Id][ids][0] #should be group level
	# 		group_name = phylogeny_structure[strain_Id][ids][1] #should be the groupName
	# 		lineage_increment = 1
	# 		parent = None
	# 		for x in list_with_levels:
	# 			if type(x) != str:
	# 				try:
	# 					parent = list_with_levels[list_with_levels.index(ids)+lineage_increment]
	# 				except:
	# 					parent = ''
	# 				break
	# 			lineage_increment += 1
	# 		if ids not in empty:
	# 			empty.append(ids)
	# 			if parent:
	# 				if re.search (r'^no_*', parent):
	# 					parent = 0
	# 				query = ("INSERT INTO `TaxGroups`(`taxgroup_id`,`groupLevel`,`groupName`,`parent`) VALUES('"+ids+"','"+group_lvl+"','"+group_name+"','"+str(parent)+"')")
	# 			else:
	# 				query = ("INSERT INTO `TaxGroups`(`taxgroup_id`,`groupLevel`,`groupName`,`parent`) VALUES('"+ids+"','"+group_lvl+"','"+group_name+"','0')")
	# 			print(query, type(parent))
	# 			cursor.execute(query)

	# for strainID in onlyStrainIDs:
	# 	for ids in phylogeny_structure[strainID].keys():
	# 		if ids != strainID:
	# 			taxgroupID = ids
	# 			query = ("INSERT INTO `Species_TaxGroup`(`strain_id`,`taxgroup_id`) VALUES('"+strainID+"','"+taxgroupID+"')")
	# 			#print(query)
	# 			cursor.execute(query)

	# cursor.execute("SET FOREIGN_KEY_CHECKS=1")
	# cnx.commit()
	# cursor.close()
	# cnx.close()

if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))