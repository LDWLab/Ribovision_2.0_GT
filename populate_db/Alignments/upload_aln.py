#!/usr/bin/env python3
import re, csv, sys, getopt, getpass, mysql.connector, argparse
from Bio import AlignIO

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Upload data to alignment related tables in DESIRE from a fasta.\nMake sure upload_accession has updated the tables with Residue data.')
	parser.add_argument('alignment_file', help='Path to alignment file', type=str)
	parser.add_argument('source', help='Defines superkingdom source (e.g. abe)', type=str)
	parser.add_argument('-aln_method','--alignment_method', help='Alignment method used (default: PROMALS3D)', type=str, default='PROMALS3D')
	parser.add_argument('-host','--db_host', help='Defines database host (default: 130.207.36.75)', type=str, default='130.207.36.75')
	parser.add_argument('-schema','--db_schema', help='Defines schema to use (default: SEREB)', type=str, default='SEREB')
	parser.add_argument('-user_name','--uname', help='Defines user name to use (default: ppenev)', type=str, default='ppenev')
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

def read_align(aln_name):
	'''
	Reads the fasta file and gets the sequences.
	'''
	alignments = AlignIO.read(open(aln_name), "fasta")
	return alignments

def superkingdom_info(cursor, ID):
	'''
	Gets the superkingdom for a strain ID
	'''
	#print(ID)
	cursor.execute("SELECT SEREB.TaxGroups.groupName FROM SEREB.Species_TaxGroup\
		INNER JOIN SEREB.TaxGroups ON SEREB.Species_TaxGroup.taxgroup_id=SEREB.TaxGroups.taxgroup_id\
		INNER JOIN SEREB.Species ON SEREB.Species_TaxGroup.strain_id=SEREB.Species.strain_id\
		WHERE SEREB.TaxGroups.groupLevel = 'superkingdom' AND SEREB.Species.strain_id = '"+ID+"'")
	results = cursor.fetchall()
	#print(ID,results)
	try:
		superkingdom=(results[0][0])
	except:
		raise ValueError ("No result for specie "+str(ID)+" in the MYSQL query!")
	return superkingdom

def check_nomo_id(cursor, occur, name):
	'''
	Gets nom_id for new name and superkingdom
	'''
	#cursor.execute("SELECT SEREB.Nomenclature.nom_id FROM SEREB.Nomenclature\
	#	INNER JOIN SEREB.Old_name ON SEREB.Nomenclature.nom_id=SEREB.Old_name.nomo_id\
	#	WHERE SEREB.Old_name.old_name = '"+name+"' AND SEREB.Old_name.N_B_Y_H_A = 'BAN' AND SEREB.Nomenclature.occurrence = '"+occur+"'")
	cursor.execute("SELECT SEREB.Nomenclature.nom_id FROM SEREB.Nomenclature\
		WHERE SEREB.Nomenclature.new_name = '"+name+"' AND SEREB.Nomenclature.occurrence = '"+occur+"'")
	result = cursor.fetchall()
	#nom_id=result[0][0]
	try:
		nom_id=result[0][0]
	except:
		raise ValueError ("No result for nom_id "+name+" and occurrence "+occur+" in the MYSQL query!")
	return nom_id

def check_polymer(cursor, taxid, nomid):
	'''
	Gets polymer id for a given taxid and nomid (LDW-prot requirement)
	'''
	cursor.execute("SELECT SEREB.Polymer_Data.PData_id FROM SEREB.Polymer_Data\
					INNER JOIN SEREB.Polymer_metadata ON SEREB.Polymer_Data.PData_id = SEREB.Polymer_metadata.polymer_id WHERE \
					SEREB.Polymer_metadata.accession_type = 'LDW-prot' AND \
					SEREB.Polymer_Data.nomgd_id = "+nomid+" AND \
					SEREB.Polymer_Data.strain_id = "+taxid)
	result = cursor.fetchall()
	try:
		pol_id=result[0][0]
	except:
		pol_id = 'NOVAL'
		print("No result for nom_id "+nomid+" and taxid "+taxid+" in the MYSQL query!")
		#raise ValueError ("No result for nom_id "+nomid+" and taxid "+taxid+" in the MYSQL query!")
	return pol_id

def upaln_getid(cursor, aln_name, source_string, method_name):
	'''
	Uploads alignment name and method, then returns its primary key from the DB.
	'''
	query = "INSERT INTO `SEREB`.`Alignment`(`Name`,`Method`,`Source`) VALUES('"+aln_name+"','"+method_name+"','"+source_string+"')"
	print(query)
	cursor.execute(query)
	lastrow_id = str(cursor.lastrowid)
	return lastrow_id

def upload_pol_aln(cursor, pol_id, aln_id):
	'''Populate table Polymer_Alignments'''
	cursor.execute("SELECT PData_id,Aln_id FROM SEREB.Polymer_Alignments WHERE\
					PData_id = '"+pol_id+"' AND\
					Aln_id = '"+aln_id+"'")
	result = cursor.fetchall()
	if len(result) == 0:
		query = "INSERT INTO `SEREB`.`Polymer_Alignments`(`PData_id`, `Aln_id`) VALUES('"+pol_id+"','"+aln_id+"')"
		cursor.execute(query)
		return True
	if len(result) == 1:
		return True
	if len(result) > 1:
		raise ValueError("Bad primary key combination in Polymer_alignments! PData_id:"+pol_id+" Aln_id: "+aln_id)

def upload_aln_data(cursor, entry, seq_aln_pos, aln_id, polymer_id):
	'''Uploads aln_data table'''
	print (entry.seq[seq_aln_pos[1]-1],end='')
	resi_id = check_resi_id(cursor, str(seq_aln_pos[0]), str(polymer_id), str(entry.seq[seq_aln_pos[1]-1]))
	#print(resi_id)
	cursor.execute("SELECT aln_id,res_id FROM SEREB.Aln_Data WHERE\
			res_id = '"+str(resi_id)+"' AND\
			aln_id = '"+str(aln_id)+"'")
	result = cursor.fetchall()
	if len(result) == 0:
		query = "INSERT INTO `SEREB`.`Aln_Data`(`aln_id`,`res_id`,`aln_pos`) VALUES('"+str(aln_id)+"','"+str(resi_id)+"','"+str(seq_aln_pos[1])+"')"
		print(query)
		cursor.execute(query)
		return True
	if len(result) == 1:
		return False

def increment_by_range(sequence, increment_range):
	'''Cretes the mapping for a set of sequence ranges'''
	aln_pos = 0
	seq_pos = 0
	seq_order = list()
	resi_mapping=[]
	for one_range in increment_range:
		for ix in range(one_range[0], one_range[1]+1):
			seq_order.append(ix)
	for seq in sequence:
		if seq == '-':
			aln_pos+=1
			continue
		resi_mapping.append((seq_order[seq_pos],aln_pos+1))
		seq_pos+=1
		aln_pos+=1
	return resi_mapping

def create_aln_mapping(entry):
	'''
	Creates list of pairs mapping sequence to alignment positions.
	'''
	sequence = entry.seq
	increment_ranges = list()
	if len(entry.id.split("_")) == 3:
		increment_ranges.append((1,len(str(sequence).replace("-",""))))
	if len(entry.id.split("_")) > 3:
		string_ranges = re.sub(r'\/.*','', entry.id.split("_")[3])
		for one_range in string_ranges.split(","):
			increment_ranges.append((int(one_range.split('-')[0]), int(one_range.split('-')[1])))

	resi_mapping = increment_by_range(sequence, increment_ranges)

	return resi_mapping

def check_resi_id(cursor, seqnum, polid, resname):
	'''
	Gets residue id for a given polymer and a sequence number.
	Also checks if the data (resname) is correct.
	'''
	query = "SELECT SEREB.Residues.resi_id, SEREB.Residues.unModResName FROM SEREB.Residues WHERE \
			SEREB.Residues.PolData_id = "+polid+" AND \
			SEREB.Residues.resNum = "+seqnum
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

def main(commandline_arguments):
	comm_args = create_and_parse_argument_options(commandline_arguments)
	aln_path = comm_args.alignment_file
	source_string = comm_args.source
	
	pw = getpass.getpass("Password: ")
	cnx = mysql.connector.connect(user=comm_args.uname, password=pw, host=comm_args.db_host, database=comm_args.db_schema)
	cursor = cnx.cursor()

	alns = read_align(aln_path)
	aln_name = aln_path.split("/")[-1].replace('_txid_tagged.fas', '').replace('.fas', '').replace('.fa', '')	#Fix that with re
	aln_id = upaln_getid(cursor, aln_name, source_string, comm_args.alignment_method)
	for entry in alns:
		print(entry.id.split('_')[2])
		superK = superkingdom_info(cursor, entry.id.split('_')[1])
		nom_id = check_nomo_id(cursor, superK[0], entry.id.split('_')[0][:-1])
		polymer_id = check_polymer(cursor, str(entry.id.split('_')[1]),str(nom_id))
		if polymer_id == 'NOVAL':
			continue
		upload_pol_aln(cursor, str(polymer_id), str(aln_id))
		mapped_resis = create_aln_mapping(entry)
		for seq_aln_pos in mapped_resis:
			upload_aln_data(cursor, entry, seq_aln_pos, aln_id, polymer_id)
		print()
	
	cnx.commit()
	cursor.close()
	cnx.close()
	print("Success!")

if __name__ == "__main__":
	sys.exit(main(sys.argv[1:]))

