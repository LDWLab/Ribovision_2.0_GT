#!/usr/bin/env python3
import getopt
import sys
import os
import csv
from Bio import Entrez
import pandas as pd
import getpass
import mysql.connector
	
def usage():
	print (\
	"USAGE:\n./taxonomy.py -i [input_file_path] -o [output_file_path] -h\n\
	-i: defines path to input file. Works only on xlsx type of files.\tREQUIRED\n\
	-o: defines output file\tREQUIRED\n\
	-h: prints this\
")

filename = './input_taxids.txt'
Entrez.email = "A.N.Other@example.com"  # Always tell NCBI who you are

with open('./input_taxids.txt') as f:
	input_taxes = f.read().splitlines()

taxon = []
verynew_taxon = []
taxonIDS = []
TAXON = {}
onlyStrainIDs = []
for i in input_taxes:
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
		aList = ['superkingdom', 'phylum', 'class', 'order', 'family' , 'genus', 'species']
		for j in sorted(aList):
			if j == i['Rank']:
				oneTaxon[i['TaxId']] = [i['Rank'],i['ScientificName']]
				if i['TaxId'] not in taxonIDS:
					taxonIDS.append(i['TaxId'])
				taxon.append([i['TaxId'],i['Rank'],i['ScientificName']])
				verynew_taxon.append([(y), i['TaxId']])
	TAXON[records[0]['TaxId']] = oneTaxon
new_taxon =  list(taxon)
noRepeats_newtaxon = []
for alist in new_taxon:
	if alist not in noRepeats_newtaxon:
		noRepeats_newtaxon.append(alist)
verynew_taxon = list(verynew_taxon)

'''-----------------------------------------------------------------------'''
'''-----------------------------------------------------------------------'''
'''-----------------------------------------------------------------------'''
uname = input('User name: ')
pw = getpass.getpass('Password: ')
cnx = mysql.connector.connect(user=uname,password=pw,host='130.207.36.75', database='SEREB2')
cursor = cnx.cursor()
'''-----------------------------------------------------------------------'''
'''-----------------------------------------------------------------------'''
'''-----------------------------------------------------------------------'''

'''---Species Table---'''################################################################################
for strain_Id in onlyStrainIDs:
	strain = TAXON[strain_Id][strain_Id][1]
	query = ("INSERT INTO `SEREB2`.`Species`(`strain_id`,`strain`) VALUES('"+strain_Id+"','"+strain+"')")
	#print(query)
	cursor.execute(query)

'''---TaxGroups Table---'''##############################################################################
empty = []
for strain_Id in onlyStrainIDs:
	for ids in TAXON[strain_Id].keys():
		group_lvl = TAXON[strain_Id][ids][0] #should be group level
		group_name = TAXON[strain_Id][ids][1] #should be the groupName
		if ids not in empty:
			empty.append(ids)
			query = ("INSERT INTO `SEREB2`.`TaxGroups`(`taxgroup_id`,`groupLevel`,`groupName`,`taxid`) VALUES('"+ids+"','"+group_lvl+"','"+group_name+"','"+ids+"')")
			#print(query)
			cursor.execute(query)

'''---Species TaxGroups Table---'''#####################################################################
for strainID in onlyStrainIDs:
	for ids in TAXON[strainID].keys():
		if ids != strainID:
			taxgroupID = ids
			query = ("INSERT INTO `SEREB2`.`Species_TaxGroup`(`strain_id`,`taxgroup_id`) VALUES('"+strainID+"','"+taxgroupID+"')")
			#print(query)
			cursor.execute(query)
cnx.commit()
cursor.close()
cnx.close()		