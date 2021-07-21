from django.http import JsonResponse, Http404
from Bio.SeqUtils import IUPACData
from alignments.models import *
from alignments.residue_api import *

def dictfetchall(cursor):
	"Return all rows from a cursor as a dict"
	columns = [col[0] for col in cursor.description]
	return [
		dict(zip(columns, row))
		for row in cursor.fetchall()
	]

def construct_query(aln_id, parent_id):
	return f'SELECT CONCAT(Aln_Data.aln_id,"_",Aln_Data.res_id) AS id,strain,unModResName,aln_pos,Species.strain_id, filtered_polymers.GI FROM Aln_Data\
		INNER JOIN Alignment ON Aln_Data.aln_id = Alignment.Aln_id\
		INNER JOIN Residues ON Aln_Data.res_id = Residues.resi_id\
		INNER JOIN (\
			SELECT * from Polymer_Data WHERE Polymer_Data.PData_id IN \
				(SELECT PData_id from Polymer_Alignments WHERE Polymer_Alignments.Aln_id = {str(aln_id)})\
			AND Polymer_Data.strain_id IN \
				(with recursive cte (taxgroup_id, groupName, parent, groupLevel) as \
		(\
		select taxgroup_id, groupName, parent, groupLevel\
			from TaxGroups\
			where parent = {str(parent_id)} or taxgroup_id = {str(parent_id)}\
			union all\
			select p.taxgroup_id, p.groupName, p.parent, p.groupLevel\
			from TaxGroups p\
			inner join cte\
				on p.parent = cte.taxgroup_id\
		)\
		select taxgroup_id from cte where (groupLevel REGEXP "strain"))) as filtered_polymers\
		ON Residues.PolData_id = filtered_polymers.PData_id\
		INNER JOIN Species ON filtered_polymers.strain_id = Species.strain_id\
		WHERE Alignment.aln_id = {str(aln_id)}'

def sql_filtered_aln_query(aln_id, parent_id):
	from django.db import connection
	SQLStatement = construct_query(aln_id, parent_id)
	with connection.cursor() as cursor:
		cursor.execute(SQLStatement)
		raw_result = dictfetchall(cursor)
	if len(raw_result) == 0:
		raise Http404("We do not have this combination of arguments in our database.")
	return raw_result

def get_fold_for_raw_result_range(raw_result_range, raw_result):
	for pos in range(raw_result_range[0],raw_result_range[1]):
		resi_data = resi_info(None, raw_result[pos]['resi_id'])
		if len(resi_data['Structural fold']) != 0:
			return resi_data['Structural fold'][0]

def para_aln(request, aln_id):
	from django.db import connection
	from alignments.views import extract_gap_only_cols, extract_species_list, construct_dict_for_json_response
	try:
		alignment = Alignment.objects.get(pk=aln_id)
	except Alignment.DoesNotExist:
		raise Http404("Alignment id "+str(aln_id)+" is not present in the database!")
	if alignment.method != 'structure_based':
		raise Http404("Alignment id "+str(aln_id)+" is not paralogous!")
	
	SQLStatement = 'SELECT CONCAT(Aln_Data.aln_id,"_",Aln_Data.res_id) AS id,resi_id,strain,unModResName,aln_pos,Species.strain_id,nomgd_id FROM Aln_Data\
		INNER JOIN Alignment ON Aln_Data.aln_id = Alignment.Aln_id\
		INNER JOIN Residues ON Aln_Data.res_id = Residues.resi_id\
		INNER JOIN Polymer_Data ON Residues.PolData_id = Polymer_Data.PData_id\
		INNER JOIN Species ON Polymer_Data.strain_id = Species.strain_id\
		WHERE Alignment.aln_id = '+str(aln_id)

	with connection.cursor() as cursor:
		cursor.execute(SQLStatement)
		raw_result = dictfetchall(cursor)
	if len(raw_result) == 0:
		raise Http404("We do not have alignment id "+str(aln_id)+" in our database.")

	currpos = 0
	startpos = 0
	ranges_of_permutation = list()
	for row in raw_result:
		if currpos == 0:
			currpos += 1
			continue
		if row['aln_pos'] < raw_result[currpos-1]['aln_pos']:
			ranges_of_permutation.append((startpos, currpos))
			startpos = currpos+1
		if currpos+1 == len(raw_result):
			ranges_of_permutation.append((startpos, currpos))
			break
		currpos+=1
	
	last_strain = ''
	fold_pattern_over_ranges, folds_with_ranges = list(), list()
	for single_range in ranges_of_permutation:
		if last_strain != '':
			if raw_result[single_range[0]]['strain'] != last_strain:
				if len(set(fold_pattern_over_ranges)) == 2:
					for i, j in enumerate(ranges_of_permutation): 
						folds_with_ranges.append((j, fold_pattern_over_ranges[i % len(fold_pattern_over_ranges)]))
				break
			fold_pattern_over_ranges.append(get_fold_for_raw_result_range(single_range, raw_result))
		else:
			fold_pattern_over_ranges.append(get_fold_for_raw_result_range(single_range, raw_result))
			last_strain = raw_result[single_range[0]]['strain']

	if len(set(fold_pattern_over_ranges)) == 1:
		last_nom = ''
		for i, single_range in enumerate(ranges_of_permutation):
			if last_nom != '':
				if raw_result[single_range[0]]['nomgd_id'] != last_nom:
					last_fold = get_fold_for_raw_result_range(single_range, raw_result)
					last_nom = raw_result[single_range[0]]['nomgd_id']
				folds_with_ranges.append((single_range, last_fold))
			else:
				last_fold = get_fold_for_raw_result_range(single_range, raw_result)
				folds_with_ranges.append((single_range, last_fold))
				last_nom = raw_result[single_range[0]]['nomgd_id']

	if len(set(fold_pattern_over_ranges)) > 2:
		raise Http404("Alignment with id "+str(aln_id)+" has more than 2 structural folds!")
	
	if len(ranges_of_permutation) % len(fold_pattern_over_ranges) != 0:
		#This will happen when the alignment has some species with permutation and others with a different fold.
		#Needs different handling
		raise Http404("Alignment with id "+str(aln_id)+" has strains with unequal fold assignments!")
	number_of_pattern_repeats = len(ranges_of_permutation)/len(fold_pattern_over_ranges)


	rawsqls = dict()
	for single_range, fold in folds_with_ranges:
		if fold not in rawsqls.keys():
			rawsqls[fold] = []
		for pos in range(single_range[0], single_range[1]):
			rawsqls[fold].append(raw_result[pos])
	
	import re
	nogap_tupaln = dict()
	max_alnposition = 0
	for fold, rawsql in rawsqls.items():
		sorted_sql = sorted(rawsql, key = lambda i: (i['strain'], i['aln_pos']))
		nogap_tupaln, max_alnposition= query_to_dict_structure(sorted_sql, fold[1].replace('_','-'), nogap_tupaln, max_alnposition)

	fastastring, frequency_list = build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_alnposition)
	
	gap_only_cols = extract_gap_only_cols(fastastring)
	filtered_spec_list = extract_species_list(fastastring)

	concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)

	response_dict = construct_dict_for_json_response([concat_fasta,filtered_spec_list,gap_only_cols,frequency_list])
	return JsonResponse(response_dict, safe = False)

def query_to_dict_structure(rawMYSQLresult, filter_element, nogap_tupaln=dict(), max_alnposition=0):
	for row in rawMYSQLresult:
		if int(row['aln_pos']) > max_alnposition:
			max_alnposition = int(row['aln_pos'])
		if (row['strain'], row['GI'], filter_element) in nogap_tupaln:
			nogap_tupaln[(row['strain'], row['GI'], filter_element)].append((row['unModResName'], row['aln_pos'], row['id'].split("_")[1]))
		else:
			nogap_tupaln[(row['strain'], row['GI'], filter_element)]=[]
			nogap_tupaln[(row['strain'], row['GI'], filter_element)].append((row['unModResName'], row['aln_pos'], row['id'].split("_")[1]))
	return nogap_tupaln, max_alnposition

def build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_alnposition, all_residues=IUPACData.protein_letters):
	import re
	from alignments.Shannon import gap_adjusted_frequency
	fasta_string=''
	alignment_rows = list()
	for strain in nogap_tupaln:
		row_residue_list = list()
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
					row_residue_list.append('-')
				mem = mem+1
			#else:
			#	raise ValueError("You are likely looking at cross-domain alignment with sequences from repeated species. For now this is not supported!")
			fasta_string+=resi_pos[0]
			row_residue_list.append(resi_pos[0])
			if index == len(nogap_tupaln[strain]):
				if resi_pos[1] < max_alnposition:
					diff = max_alnposition-resi_pos[1]
					for index2,i in enumerate(range(0,diff), start=1):
						fasta_string+='-'
						row_residue_list.append('-')
						if index2 == diff:
							alignment_rows.append(row_residue_list)
							fasta_string+='\n'
				else:
					alignment_rows.append(row_residue_list)
					fasta_string+='\n'
	alignment_columns = [x for x in map(list, zip(*alignment_rows))]
	frequency_list = list()
	for resi_column in alignment_columns:
		frequency_list.append(gap_adjusted_frequency(resi_column, all_residues))
	literal_string = re.sub(r'\n\n','\n',fasta_string,flags=re.M)
	return literal_string.lstrip().encode('unicode-escape').decode('ascii'), frequency_list
