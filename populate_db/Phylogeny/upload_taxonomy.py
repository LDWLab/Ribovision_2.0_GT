#!/usr/bin/env python3
import getopt
import sys, re
import os
import csv
from Bio import Entrez
import pandas as pd
import getpass
import mysql.connector
import time

#def usage():
#	print (\
#	"USAGE:\n./taxonomy.py -i [input_file_path] -o [output_file_path] -h\n\
#	-i: defines path to input file. Works only on xlsx type of files.\tREQUIRED\n\
#	-o: defines output file\tREQUIRED\n\
#	-h: prints this\
#")

filename = './data/input_taxids.txt'
Entrez.email = "ppenev@gatech.edu"  # Always tell NCBI who you are

with open(filename) as f:
	input_taxes = f.read().splitlines()

'''-----------------------------------------------------------------------'''
'''-----------------------------------------------------------------------'''
'''-----------------------------------------------------------------------'''
uname = input('User name: ')
pw = getpass.getpass('Password: ')
cnx = mysql.connector.connect(user=uname,password=pw,host='130.207.36.76', database='DESIRE')
cursor = cnx.cursor()
'''-----------------------------------------------------------------------'''
'''-----------------------------------------------------------------------'''
'''-----------------------------------------------------------------------'''

taxon = []
verynew_taxon = []
taxonIDS = []
TAXON = {}
onlyStrainIDs = []
aList = ['superkingdom', 'phylum', 'class', 'order', 'family' , 'genus', 'species']
aList.reverse()
for i in input_taxes:
	time.sleep(1)
	oneTaxon = {}
	y = i
	handle = Entrez.efetch(db="Taxonomy", id=i, retmode="xml")
	records = Entrez.read(handle) #records contains the strain
	oneTaxon[records[0]['TaxId']] = ['strain',records[0]['ScientificName']]
	if records[0]['TaxId'] not in taxonIDS:
		taxonIDS.append(records[0]['TaxId'])
		onlyStrainIDs.append(records[0]['TaxId'])
	lineages = records[0]["LineageEx"] #want more than just lineages, want strain
	for i in lineages:
		dict_Tax={}
		for j in sorted(aList):
			if j == i['Rank']:
				oneTaxon[i['TaxId']] = [i['Rank'],i['ScientificName']]
				if i['TaxId'] not in taxonIDS:
					taxonIDS.append(i['TaxId'])
				taxon.append([i['TaxId'],i['Rank'],i['ScientificName']])
				verynew_taxon.append([(y), i['TaxId']])
	TAXON[records[0]['TaxId']] = oneTaxon
	print(oneTaxon)
new_taxon =  list(taxon)
noRepeats_newtaxon = []
for alist in new_taxon:
	if alist not in noRepeats_newtaxon:
		noRepeats_newtaxon.append(alist)
verynew_taxon = list(verynew_taxon)


'''---Species Table---'''################################################################################
for strain_Id in onlyStrainIDs:
	strain = TAXON[strain_Id][strain_Id][1]
	query = ("INSERT INTO `DESIRE`.`Species`(`strain_id`,`strain`) VALUES('"+strain_Id+"','"+strain+"')")
	#print(query)
	cursor.execute(query)

'''---TaxGroups Table---'''##############################################################################
empty = []
aList.insert(0,'strain')
cursor.execute("SET FOREIGN_KEY_CHECKS=0")
for strain_Id in onlyStrainIDs:
	list_with_levels=['no_strain','no_species','no_genus','no_family','no_order','no_class','no_phylum','no_superkingdom']
	for ids in TAXON[strain_Id].keys():
		list_with_levels[aList.index(TAXON[strain_Id][ids][0])] = ids
	for ids in TAXON[strain_Id].keys():
		group_lvl = TAXON[strain_Id][ids][0] #should be group level
		group_name = TAXON[strain_Id][ids][1] #should be the groupName
		lineage_increment = 1
		parent = None
		for x in list_with_levels:
			if type(x) != str:
				try:
					parent = list_with_levels[list_with_levels.index(ids)+lineage_increment]
				except:
					parent = ''
				break
			lineage_increment += 1
		if ids not in empty:
			empty.append(ids)
			if parent:
				if re.search (r'^no_*', parent):
					parent = 0
				query = ("INSERT INTO `DESIRE`.`TaxGroups`(`taxgroup_id`,`groupLevel`,`groupName`,`parent`) VALUES('"+ids+"','"+group_lvl+"','"+group_name+"','"+str(parent)+"')")
			else:
				query = ("INSERT INTO `DESIRE`.`TaxGroups`(`taxgroup_id`,`groupLevel`,`groupName`,`parent`) VALUES('"+ids+"','"+group_lvl+"','"+group_name+"','0')")
			print(query, type(parent))
			cursor.execute(query)

for strainID in onlyStrainIDs:
	for ids in TAXON[strainID].keys():
		if ids != strainID:
			taxgroupID = ids
			query = ("INSERT INTO `DESIRE`.`Species_TaxGroup`(`strain_id`,`taxgroup_id`) VALUES('"+strainID+"','"+taxgroupID+"')")
			#print(query)
			cursor.execute(query)

cursor.execute("SET FOREIGN_KEY_CHECKS=1")
cnx.commit()
cursor.close()
cnx.close()