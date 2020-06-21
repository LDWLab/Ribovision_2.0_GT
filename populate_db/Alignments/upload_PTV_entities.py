#!/usr/bin/env python3
import csv, sys, getopt, getpass, mysql.connector
from Bio import AlignIO


def usage():
	print (\
	"USAGE:\nupload_PTV_entities.py -t [tsv_file_path]-h\n\
	-t: defines path to tab separated file.\tREQUIRED\n\
	-h: prints this\
")

try:
	opts, args = getopt.getopt(sys.argv[1:], 't:h', ['tsv_file=', 'help'])
except getopt.GetoptError:
	usage()
	sys.exit(2)

for opt, arg in opts:
	if opt in ('-h', '--help'):
		usage()
		sys.exit(2)
	elif opt in ('-t', '--tsv_file'):
		tsv_path = arg
	else:
		usage()
		sys.exit(2)

#uname = input("User name: ")
pw = getpass.getpass("Password: ")
cnx = mysql.connector.connect(user='ppenev', password=pw, host='130.207.36.76', database='SEREB2')
cursor = cnx.cursor()

def read_tsv(tsv_path):
	'''
	Reads the tsv file and returns the data.
	'''
	with open(tsv_path, 'r') as tsv_file:
		reader = csv.reader(tsv_file, delimiter='\t')
		tsv_list = list(reader)
	return tsv_list

def strain_info(pdb_id):
	'''Gets the strain id for a pdbid'''
	cursor.execute("")
	results = cursor.fetchall()
	try:
		strain_id = (results[0][0])
	except:
		raise ValueError("No results for PDB ID "+str(pdb_id)+" in the MYSQL query!")
	return strain_id

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
		raise ValueError ("No result for specie "+str(ID)+" in the MYSQL query!")
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
		raise ValueError ("No result for nom_id "+name+" and occurrence "+occur+" in the MYSQL query!")
	return nom_id

def check_polymer(taxid, nomid):
	'''
	Gets polymer id for a given taxid and nomid (LDW-prot requirement)
	'''
	cursor.execute("SELECT SEREB2.Polymer_Data.PData_id FROM SEREB2.Polymer_Data\
		INNER JOIN SEREB2.Polymer_metadata ON SEREB2.Polymer_Data.PData_id = SEREB2.Polymer_metadata.polymer_id\
		WHERE SEREB2.Polymer_metadata.accession_type = 'LDW-prot' AND SEREB2.Polymer_Data.nomgd_id = "+nomid+" AND SEREB2.Polymer_Data.strain_id = "+taxid)
	result = cursor.fetchall()
	try:
		pol_id=result[0][0]
	except:
		pol_id = 'NOVAL'
		print("No result for nom_id "+nomid+" and taxid "+taxid+" in the MYSQL query!")
		#raise ValueError ("No result for nom_id "+nomid+" and taxid "+taxid+" in the MYSQL query!")
	return pol_id


def create_aln_mapping(sequence):
	'''
	Creates list of pairs mapping sequence to alignment positions.
	'''
	resNum = 0
	aln_pos = 0
	resi_mapping=[]
	for residue in sequence:
		aln_pos+=1
		if residue == '-':
			pass
		else:
			resNum+=1
			resi_mapping.append((resNum,aln_pos))
	return resi_mapping

def check_resi_id(seqnum, polid, resname):
	'''
	Gets residue id for a given polymer and a sequence number.
	Also checks if the data (resname) is correct.
	'''
	query = "SELECT SEREB2.Residues.resi_id, SEREB2.Residues.unModResName FROM SEREB2.Residues WHERE SEREB2.Residues.PolData_id = "+polid+" AND SEREB2.Residues.resNum = "+seqnum
	cursor.execute(query)
	result = cursor.fetchall()
	try:
		res_id=result[0][0]
	except:
		raise ValueError ("No result for pol_id "+polid+" and sequence number "+seqnum+" in the MYSQL query!")
	if result[0][1] == resname:
		return res_id
	else:
		raise ValueError ("Wrong residue name for polymer id "+polid+" and residue "+seqnum+". Is there something wrong with the numbering?")
	


def main():
	entities = read_tsv(tsv_path)

	for entry in entities:
		print(entry.id.split('_')[2])
		superK = superkingdom_info(entry.id.split('_')[1])
		nom_id = check_nomo_id(superK[0], entry.id.split('_')[0][:-1])
		polymer_id = check_polymer(str(entry.id.split('_')[1]),str(nom_id))
		if polymer_id == 'NOVAL':
			continue
		mapped_resis = create_aln_mapping(entry.seq)
		for seq_aln_pos in mapped_resis:
			print (entry.seq[seq_aln_pos[1]-1],end='')
			resi_id = check_resi_id(str(seq_aln_pos[0]), str(polymer_id), str(entry.seq[seq_aln_pos[1]-1]))
			#print(resi_id)
			query = "INSERT INTO `SEREB2`.`Aln_Data`(`aln_id`,`res_id`,`aln_pos`) VALUES('"+str(aln_id)+"','"+str(resi_id)+"','"+str(seq_aln_pos[1])+"')"
			#print(query)
			cursor.execute(query)
		print()

	

if __name__ == "__main__":
	main()

cnx.commit()
cursor.close()
cnx.close()
print("Success!")