#!/usr/bin/env python3
import csv
from collections import Counter
import numpy as np
import itertools 
import mysql.connector
import getpass

'''
Generating a Nomenclature table and a Old_name table. The data is new_name of the protein 
and the occurrence of the protein either in bacteria, archaea, or eukarya. 
If it is eukarya, whether it is available in mitochondria or chloroplast.
The table can be directly exported to the database. The write_csv function checks for the occurrence and inputs
the corresponding old name for each of the occurrence. Please pay special attention 
to the code used to upload the output to the database. 
''' 

uname = input("User name: ")
pw = getpass.getpass("Password: ")
cnx = mysql.connector.connect(user=uname, password=pw, host='130.207.36.76', database='DESIRE')
cursor = cnx.cursor()

def read_csv (csvfile):
	''' 
	Reads csv file with new names, occurences and old names,
	then returns a list of lists for the rows.
	'''
	csvReader = csv.reader(csvfile, delimiter= ',')
	next(csvReader)
	occurList = []
	csvList=[]
	for row in csvReader:
		csvList.append(row)
	return csvList

def execute_old_name_query(row, occurrence, lastrow_id):
	'''
	Checks the current occurence and executes queries for the Old_name table.
	'''
	if occurrence == 'B':
		if row[3]:
			query = ("INSERT INTO `DESIRE`.`Old_name`(`nomo_id`,`old_name`,`N_B_Y_H_A`) VALUES('"+str(lastrow_id)+"','"+row[3]+"','"+occurrence+"')")
			cursor.execute(query)
	elif occurrence == 'A':
		if row[6]:
			query = ("INSERT INTO `DESIRE`.`Old_name`(`nomo_id`,`old_name`,`N_B_Y_H_A`) VALUES('"+str(lastrow_id)+"','"+row[6]+"','"+occurrence+"')")
			cursor.execute(query)
	elif occurrence == 'E':
		if row[4]:
			query = ("INSERT INTO `DESIRE`.`Old_name`(`nomo_id`,`old_name`,`N_B_Y_H_A`) VALUES('"+str(lastrow_id)+"','"+row[4]+"','Y')")
			cursor.execute(query)
		if row[5]:
			query = ("INSERT INTO `DESIRE`.`Old_name`(`nomo_id`,`old_name`,`N_B_Y_H_A`) VALUES('"+str(lastrow_id)+"','"+row[5]+"','H')")
			cursor.execute(query)
	elif occurrence == 'm':
		if row[3]:
			query = ("INSERT INTO `DESIRE`.`Old_name`(`nomo_id`,`old_name`,`N_B_Y_H_A`) VALUES('"+str(lastrow_id)+"','"+row[3]+"','B')")
			cursor.execute(query)
		if row[4]:
			query = ("INSERT INTO `DESIRE`.`Old_name`(`nomo_id`,`old_name`,`N_B_Y_H_A`) VALUES('"+str(lastrow_id)+"','"+row[4]+"','mY')")
			cursor.execute(query)
		if row[5]:
			query = ("INSERT INTO `DESIRE`.`Old_name`(`nomo_id`,`old_name`,`N_B_Y_H_A`) VALUES('"+str(lastrow_id)+"','"+row[5]+"','mH')")
			cursor.execute(query)
	elif occurrence == 'c':
		if row[3]:
			query = ("INSERT INTO `DESIRE`.`Old_name`(`nomo_id`,`old_name`,`N_B_Y_H_A`) VALUES('"+str(lastrow_id)+"','"+row[3]+"','B')")
			cursor.execute(query)
		if row[4]:
			query = ("INSERT INTO `DESIRE`.`Old_name`(`nomo_id`,`old_name`,`N_B_Y_H_A`) VALUES('"+str(lastrow_id)+"','"+row[4]+"','cY')")
			cursor.execute(query)
	else:
		raise ValueError('Occurrence can only be B, A, E, m, or c!')
	return True

def upload_csv(csvList, cursor):
	lastrow_id=0
	for row in csvList:
		for occurrence in list(row[2]):
			if row[0][1] == 'S':
				query = ("INSERT INTO `DESIRE`.`Nomenclature`(`new_name`,`occurrence`,`MoleculeGroup`) VALUES('"+row[0]+"','"+occurrence+"','SSU')")
			elif row[0][1] == 'L':
				query = ("INSERT INTO `DESIRE`.`Nomenclature`(`new_name`,`occurrence`,`MoleculeGroup`) VALUES('"+row[0]+"','"+occurrence+"','LSU')")
			else:
				raise ValueError('Incorrect Molecule group! Can be only L or S: '+row[0])
			cursor.execute(query)
			lastrow_id = cursor.lastrowid
			execute_old_name_query(row, occurrence, lastrow_id)
			query = ("INSERT INTO `DESIRE`.`Old_name`(`nomo_id`,`old_name`,`N_B_Y_H_A`) VALUES('"+str(lastrow_id)+"','"+row[1]+"','BAN')")
			cursor.execute(query)
	cnx.commit()
	return True

with open("./ALL-BAN-NICK.csv") as csvfile:
	csvList = read_csv(csvfile)
	upload_csv(csvList, cursor)
	
cursor.close()
cnx.close()