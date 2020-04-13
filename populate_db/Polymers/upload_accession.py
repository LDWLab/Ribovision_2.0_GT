#!/usr/bin/env python3
import csv, sys, getopt, getpass, mysql.connector

def usage():
	print (\
	"USAGE:\n./upload_accession.py -c [csv_file_path]-h\n\
	-c: defines path to csv file with txids, accessions, database, protein name, description, and sequence.\tREQUIRED\n\
	-h: prints this\
")

try:
	opts, args = getopt.getopt(sys.argv[1:], 'c:h', ['alignment=', 'help'])
except getopt.GetoptError:
	usage()
	sys.exit(2)

for opt, arg in opts:
	if opt in ('-h', '--help'):
		usage()
		sys.exit(2)
	elif opt in ('-c', '--alignment'):
		csv_path = arg
	else:
		usage()
		sys.exit(2)

#uname = input("User name: ")
pw = getpass.getpass("Password: ")
cnx = mysql.connector.connect(user='ppenev', password=pw, host='130.207.36.76', database='SEREB2')
cursor = cnx.cursor()

def read_csv(csv_path):
	with open(csv_path, 'r') as csv_file:
		reader = csv.reader(csv_file)
		csv_list = list(reader)
	return csv_list

def superkingdom_info(ID):
	'''
	Gets the superkingdom for a strain ID
	'''
	#print(ID)
	cursor.execute("SELECT SEREB2.TaxGroups.groupName FROM SEREB2.Species_TaxGroup\
		INNER JOIN SEREB2.TaxGroups ON SEREB2.Species_TaxGroup.taxgroup_id=SEREB2.TaxGroups.taxgroup_id\
		INNER JOIN SEREB2.Species ON SEREB2.Species_TaxGroup.strain_id=SEREB2.Species.strain_id\
		WHERE SEREB2.TaxGroups.groupLevel = 'superkingdom' AND SEREB2.Species.strain_id = '"+ID+"'")
	results = cursor.fetchall()
	#print(ID,results)
	try:
		superkingdom=(results[0][0])
	except:
		raise ValueError ("No result for specie "+str(ID)+" in the MYSQL query")
	return superkingdom

def check_nomo_id(occur, name):
	'''
	Gets nom_id for new name and superkingdom
	'''
	#cursor.execute("SELECT SEREB2.Nomenclature.nom_id FROM SEREB2.Nomenclature\
	#	INNER JOIN SEREB2.Old_name ON SEREB2.Nomenclature.nom_id=SEREB2.Old_name.nomo_id\
	#	WHERE SEREB2.Old_name.old_name = '"+name+"' AND SEREB2.Old_name.N_B_Y_H_A = 'BAN' AND SEREB2.Nomenclature.occurrence = '"+occur+"'")
	cursor.execute("SELECT SEREB2.Nomenclature.nom_id FROM SEREB2.Nomenclature\
		WHERE SEREB2.Nomenclature.new_name = '"+name+"' AND SEREB2.Nomenclature.occurrence = '"+occur+"'")
	result = cursor.fetchall()
	#nom_id=result[0][0]
	try:
		nom_id=result[0][0]
	except:
		raise ValueError ("No result for nom_id "+name+" and occurrence "+occur+" in the MYSQL query")
	return nom_id

def upload_resi(poldata_id, fullseq):
	i = 1
	for resi in fullseq:
		query = "INSERT INTO `SEREB2`.`Residues`(`PolData_id`,`resNum`,`unModResName`) VALUES('"+poldata_id+"','"+str(i)+"','"+resi+"')"
		cursor.execute(query)
		#print(query)
		i+=1
	return True

def main():
	csv_list = read_csv(csv_path)
	for entry in csv_list:
		superK = superkingdom_info(entry[0])
		nom_id = check_nomo_id(superK[0], entry[3])	
		query = "INSERT INTO `SEREB2`.`Polymer_Data`(`GI`,`strain_ID`,`nomgd_id`, `GeneDescription`) VALUES('"+entry[1]+"','"+str(entry[0])+"','"+str(nom_id)+"','"+entry[4]+"')"
		print(query)
		cursor.execute(query)
		lastrow_id = str(cursor.lastrowid)
		query = "INSERT INTO `SEREB2`.`Polymer_metadata`(`polymer_id`,`accession_type`,`polymer_type`, `Fullseq`) VALUES('"+str(lastrow_id)+"','LDW-prot','protein','"+entry[5]+"')"
		cursor.execute(query)
		#print(query)
		upload_resi(str(lastrow_id), entry[5])
	

if __name__ == "__main__":
	main()

cnx.commit()
cursor.close()
cnx.close()
print("Success!")