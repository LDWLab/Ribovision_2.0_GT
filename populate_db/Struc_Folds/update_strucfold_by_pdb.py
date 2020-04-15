#!/usr/bin/env python3
import re, csv, sys, getopt, getpass, mysql.connector, argparse

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Update structural folds tables from ECOD and SCOPe, given a PDB ID')
	parser.add_argument('-pdb','--pdb_id', help='PDB identifier to query', type=str)
	parser.add_argument('-uname','--user_name', help='Username for connecting to DESIRE', type=str)
	parser.add_argument('-dl','--download_most_recent_fold_definitions', help='Update latest fold definitions.', required=False, default=False, action="store_true")
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

def download_latest_fold_defs(url, file_name):
	'''Download latest structural phylogeny file'''
	import urllib.request
	try:
		urllib.request.urlretrieve(url, file_name)
	except:
		raise ConnectionError("Failed downloading "+ url)
	return True

def initiate_connection(uname, host, database):
	pw = getpass.getpass("Password: ")
	cnx = mysql.connector.connect(user=uname, password=pw, host=host, database=database)
	return cnx

def parse_ecod_definitions(pdbid):
	ecod_defs = list()
	with open("ecod.latest.domains.txt") as file:
		for line in file:
			if re.match(r'^#', line):
				continue
			if pdbid == line.split('\t')[4]:
				ecod_defs.append([line.split('\t')[3],line.split('\t')[5],line.split('\t')[7],line.split('\t')[8],
				line.split('\t')[9],line.split('\t')[10],line.split('\t')[11],line.split('\t')[12],line.split('\t')[13]])
	if len(ecod_defs) == 0:
		raise ValueError("Entered pdb_id: "+ pdbid+" was not found in the latest ecod file!")
	return ecod_defs

def upload_struc_fold(level, name, class_sys, parent, external_id, cursor, cnx):
	if level != 'Architecture':
		query = "INSERT INTO `SEREB`.`Structural_Folds`(`Level`,`Name`,`classification_system`,`parent`,`external_id`) VALUES('"+level+"','"+name+"','"+class_sys+"','"+parent+"','"+external_id+"')"
	else:
		query = "INSERT INTO `SEREB`.`Structural_Folds`(`Level`,`Name`,`classification_system`,`external_id`) VALUES('"+level+"','"+name+"','"+class_sys+"','"+external_id+"')"
	print(query)
	cursor.execute(query)
	lastrow_id = str(cursor.lastrowid)
	cnx.commit()
	return lastrow_id

def check_then_upload_struc_fold(ecod_definitions, cursor, cnx):
	levels = ['Architecture', 'X', 'H', 'T', 'F']
	for entry_for_upload in ecod_definitions:
		parent = "NULL"
		last_statement = "IS "+parent
		external_ecod_id = ['0'] + entry_for_upload[0].split('.')
		for fold_name,fold_level,fold_id in zip(entry_for_upload[-5:], levels, external_ecod_id):
			fold_name = fold_name.replace("\"","")
			if fold_level != 'Architecture':
				last_statement = "= '"+parent+"'"
			cursor.execute("SELECT SEREB.Structural_Folds.struc_fold_id FROM SEREB.Structural_Folds WHERE\
							SEREB.Structural_Folds.Level = '"+fold_level+"' AND\
							SEREB.Structural_Folds.Name = '"+fold_name+"' AND\
							SEREB.Structural_Folds.external_id = '"+fold_id+"' AND\
							SEREB.Structural_Folds.parent "+last_statement)
			result = cursor.fetchall()
			if len(result) == 0:
				parent = upload_struc_fold(fold_level, fold_name, 'ECOD', parent, fold_id, cursor, cnx)
			if len(result) == 1:
				parent = str(result[0][0])
			

def main(commandline_arguments):
	comm_args = create_and_parse_argument_options(commandline_arguments)
	if comm_args.download_most_recent_fold_definitions:
		download_latest_fold_defs('http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt', 'ecod.latest.domains.txt')
	ecod_defs = parse_ecod_definitions(str(comm_args.pdb_id).lower())

	cnx = initiate_connection(comm_args.user_name, '130.207.36.75', 'SEREB')
	cursor = cnx.cursor()

	check_then_upload_struc_fold(ecod_defs, cursor, cnx)

	cursor.close()
	cnx.close()

if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))