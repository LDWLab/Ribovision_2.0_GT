#!/usr/bin/env python3
import sys, re, os, csv, getpass, argparse, mysql.connector, time, json
from Bio import Entrez
import pandas as pd
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Update species taxonomy given a list of taxids')
	parser.add_argument('input_taxids', help='Path to file with taxids, each taxid should be on a single line.', type=str)
	parser.add_argument('phylogeny_file', help='Use this file to store and read phylogeny.', type=str)
	parser.add_argument('user_name', help='Username for connecting to DESIRE', type=str)
	parser.add_argument('-pw','--password', help='Defines user password to use', type=str)
	parser.add_argument('-host','--db_host', help='Defines database host (default: 130.207.36.76)', type=str, default='130.207.36.76')
	parser.add_argument('-schema','--db_schema', help='Defines schema to use (default: DESIRE)', type=str, default='DESIRE')
	parser.add_argument('-dl','--download_most_recent_phylogeny', help='Update latest phylogeny.', default=False, action="store_true")
	parser.add_argument('-commit','--commit_changes', help='Commit the changes to the DB', action="store_true")
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

def download_phylogeny(input_taxids, newfileLoc, phylogeny_levels):
	Entrez.email = "ppenev@gatech.edu"  # Always tell NCBI who you are
	phylogeny_structure = dict()
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
	json.dump(phylogeny_structure, open(newfileLoc, 'w' ))
	return phylogeny_structure

def makeLevelTxid(txidDict):
	lvTotxid = dict()
	for k,v in txidDict.items():
		lvTotxid[v[0]] = k
	return lvTotxid

def getParentRecurse(phylogenyLevels, lvTotxid):
	if len(phylogenyLevels) == 0:
		return 0
	while phylogenyLevels[-1] not in lvTotxid:
		return getParentRecurse(phylogenyLevels[:phylogenyLevels.index(phylogenyLevels[-1])], lvTotxid)
	return lvTotxid[phylogenyLevels[-1]]

def constructParentsDict(txidDict, lvToTxid, phylogenyLevels):
	txidToParent = dict()
	for txid, levelName in txidDict.items():
		txidToParent[txid] = getParentRecurse(phylogenyLevels[:phylogenyLevels.index(levelName[0])], lvToTxid)
	return txidToParent

def main(commandline_arguments):
	comm_args = create_and_parse_argument_options(commandline_arguments)
	input_taxids = load_species_taxids(comm_args.input_taxids)
	phylogeny_levels = ['superkingdom', 'phylum', 'class', 'order', 'family' , 'genus', 'species', 'strain']
	
	if comm_args.download_most_recent_phylogeny:
		phylogeny_structure = download_phylogeny(input_taxids, comm_args.phylogeny_file, phylogeny_levels)
	else:
		phylogeny_structure = json.load( open( comm_args.phylogeny_file ) )

	for txid in input_taxids:
		lvTotxid = makeLevelTxid(phylogeny_structure[txid])
		parentDict = constructParentsDict(phylogeny_structure[txid], lvTotxid, phylogeny_levels)
		for k in phylogeny_structure[txid].keys():
			phylogeny_structure[txid][k].append(parentDict[k])

	cnx = initiate_connection(comm_args.user_name, comm_args.db_host, comm_args.db_schema)
	cursor = cnx.cursor()
	# '''---Species Table---'''################################################################################
	for strain_Id in input_taxids:
		strain = phylogeny_structure[strain_Id][strain_Id][1]
		query = ("INSERT INTO `Species`(`strain_id`,`strain`) VALUES('"+strain_Id+"','"+strain+"')")
		print(query)
		cursor.execute(query)

	#Fix the recursive relationship
	'''---TaxGroups Table---'''##############################################################################
	phylogeny_levels.insert(0,'strain')
	cursor.execute("SET FOREIGN_KEY_CHECKS=0")
	for strain_Id in input_taxids:
		list_with_levels=['no_strain','no_species','no_genus','no_family','no_order','no_class','no_phylum','no_superkingdom']
		for ids in phylogeny_structure[strain_Id].keys():
			list_with_levels[phylogeny_levels.index(phylogeny_structure[strain_Id][ids][0])] = ids
		for ids in phylogeny_structure[strain_Id].keys():
			group_lvl = phylogeny_structure[strain_Id][ids][0] #should be group level
			group_name = phylogeny_structure[strain_Id][ids][1] #should be the groupName
			parent = phylogeny_structure[strain_Id][ids][2] #should be the parent
			query = ("INSERT INTO `TaxGroups`(`taxgroup_id`,`groupLevel`,`groupName`,`parent`) VALUES('"+ids+"','"+group_lvl+"','"+group_name+"','"+str(parent)+"')")
			print(query, type(parent))
			try:
				cursor.execute(query)
			except:
				print("Skipped")

	for strainID in input_taxids:
		for ids in phylogeny_structure[strainID].keys():
			if ids != strainID:
				taxgroupID = ids
				query = ("INSERT INTO `Species_TaxGroup`(`strain_id`,`taxgroup_id`) VALUES('"+strainID+"','"+taxgroupID+"')")
				print(query)
				try:
					cursor.execute(query)
				except:
					print("Skipped")

	cursor.execute("SET FOREIGN_KEY_CHECKS=1")
	if comm_args.commit_changes:
		cnx.commit()
	cursor.close()
	cnx.close()

if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))