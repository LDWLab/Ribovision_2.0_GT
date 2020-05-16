from django.http import JsonResponse
from alignments.models import *

def dictfetchall(cursor):
	"Return all rows from a cursor as a dict"
	columns = [col[0] for col in cursor.description]
	return [
		dict(zip(columns, row))
		for row in cursor.fetchall()
	]

def sql_filtered_aln_query(aln_id, parent_id):
	from django.db import connection
	SQLStatement = 'SELECT CONCAT(Aln_Data.aln_id,"_",Aln_Data.res_id) AS id,strain,unModResName,aln_pos,Species.strain_id FROM SEREB.Aln_Data\
		INNER JOIN SEREB.Alignment ON SEREB.Aln_Data.aln_id = SEREB.Alignment.Aln_id\
		INNER JOIN SEREB.Residues ON SEREB.Aln_Data.res_id = SEREB.Residues.resi_id\
		INNER JOIN (\
			SELECT * from SEREB.Polymer_Data WHERE SEREB.Polymer_Data.PData_id IN \
				(SELECT PData_id from SEREB.Polymer_Alignments WHERE SEREB.Polymer_Alignments.Aln_id = '+str(aln_id)+')\
			AND SEREB.Polymer_Data.strain_id IN \
				(SELECT strain_id FROM SEREB.Species_TaxGroup WHERE taxgroup_id = '+str(parent_id)+')) as filtered_polymers\
		ON SEREB.Residues.PolData_id = filtered_polymers.PData_id\
		INNER JOIN SEREB.Species ON filtered_polymers.strain_id = SEREB.Species.strain_id\
		WHERE SEREB.Alignment.aln_id = '+str(aln_id)
	with connection.cursor() as cursor:
		cursor.execute(SQLStatement)
		raw_result = dictfetchall(cursor)
	if len(raw_result) == 0:
		raise Http404("We do not have this combination of arguments in our database.")
	return raw_result

def query_to_dict_structure(rawMYSQLresult, filter_element, nogap_tupaln=dict(), max_alnposition=0):
	for row in rawMYSQLresult:
		if int(row['aln_pos']) > max_alnposition:
			max_alnposition = int(row['aln_pos'])
		if (row['strain'], filter_element) in nogap_tupaln:
			nogap_tupaln[(row['strain'], filter_element)].append((row['unModResName'], row['aln_pos'], row['id'].split("_")[1]))
		else:
			nogap_tupaln[(row['strain'], filter_element)]=[]
			nogap_tupaln[(row['strain'], filter_element)].append((row['unModResName'], row['aln_pos'], row['id'].split("_")[1]))
	return nogap_tupaln, max_alnposition

def build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_alnposition):
	import re
	fasta_string=''
	for strain in nogap_tupaln:
		alignment_sequence_name = str(strain[1])+"_"+str(re.sub(' ','_',strain[0]))
		fasta_string+='\n>'+alignment_sequence_name+'\n'
		mem = 1
		for index, resi_pos in enumerate(nogap_tupaln[strain], start=1):
			if mem == resi_pos[1]:
				mem = mem+1
			elif mem < resi_pos[1]:
				diff = resi_pos[1]-mem
				for i in range(0,diff):
					mem = mem+1
					fasta_string+='-'
				mem = mem+1
			#else:
			#	raise ValueError("You are likely looking at cross-domain alignment with sequences from repeated species. For now this is not supported!")
			fasta_string+=resi_pos[0]
			if index == len(nogap_tupaln[strain]):
				if resi_pos[1] < max_alnposition:
					diff = max_alnposition-resi_pos[1]
					for index2,i in enumerate(range(0,diff), start=1):
						fasta_string+='-'
						if index2 == diff:
							fasta_string+='\n'
				else:
					fasta_string+='\n'
	literal_string = re.sub(r'\n\n','\n',fasta_string,flags=re.M)
	return literal_string.lstrip().encode('unicode-escape').decode('ascii')
