import re

from django.shortcuts import render
from django.http import HttpResponse, Http404, JsonResponse
from django.urls import reverse_lazy
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.views.generic import ListView, CreateView, UpdateView

from alignments.models import *
from alignments.taxonomy_views import *


def sql_alignment_query(aln_id):
	alnposition = AlnData.objects.raw('SELECT * FROM SEREB.Aln_Data\
		INNER JOIN SEREB.Alignment ON SEREB.Aln_Data.aln_id = SEREB.Alignment.Aln_id\
		INNER JOIN SEREB.Residues ON SEREB.Aln_Data.res_id = SEREB.Residues.resi_id\
		INNER JOIN SEREB.Polymer_Data ON SEREB.Residues.PolData_id = SEREB.Polymer_Data.PData_id\
		INNER JOIN SEREB.Species ON SEREB.Polymer_Data.strain_id = SEREB.Species.strain_id\
		WHERE SEREB.Alignment.aln_id = '+str(aln_id)+'')
	alnpos=[]
	fastastring,max_aln_length = build_alignment(alnposition)

	return fastastring,max_aln_length

def sql_filtered_aln_query(aln_id, parent_id):
	SQLStatement = 'SELECT * FROM SEREB.Aln_Data\
					INNER JOIN SEREB.Alignment ON SEREB.Aln_Data.aln_id = SEREB.Alignment.Aln_id\
					INNER JOIN SEREB.Residues ON SEREB.Aln_Data.res_id = SEREB.Residues.resi_id\
					INNER JOIN (SELECT * from SEREB.Polymer_Data WHERE \
								SEREB.Polymer_Data.PData_id IN (SELECT PData_id from SEREB.Polymer_Alignments WHERE SEREB.Polymer_Alignments.Aln_id = '+str(aln_id)+')\
								AND \
								SEREB.Polymer_Data.strain_id IN (SELECT strain_id FROM SEREB.Species_TaxGroup WHERE taxgroup_id = '+str(parent_id)+')) as filtered_polymers\
								ON SEREB.Residues.PolData_id = filtered_polymers.PData_id\
					INNER JOIN SEREB.Species ON filtered_polymers.strain_id = SEREB.Species.strain_id\
					WHERE SEREB.Alignment.aln_id = '+str(aln_id)
	alnposition = AlnData.objects.raw(SQLStatement)
	if len(alnposition) == 0:
		raise Http404("We do not have this combination of arguments in our database.")
	fastastring,max_aln_length = build_alignment(alnposition)
	return fastastring,max_aln_length

def dictfetchall(cursor):
	"Return all rows from a cursor as a dict"
	columns = [col[0] for col in cursor.description]
	return [
		dict(zip(columns, row))
		for row in cursor.fetchall()
	]

def sql_filtered_aln_query_two_parents(aln_id, parent1_id, parent2_id):
	'''Queries DB for alignment from 2 parent taxids AND
	Constructs the alignment with parent group names at the front of sequence names.
	'''
	parent1_level = Taxgroups.objects.filter(taxgroup_id=parent1_id)[0].grouplevel
	parent2_level = Taxgroups.objects.filter(taxgroup_id=parent2_id)[0].grouplevel
	if parent1_level != parent2_level:
		raise Http404("For now we do not support comparisons between different taxonomic levels. Offending levels are: "+parent1_level+" and "+parent2_level)
	from django.db import connection
	SQLStatement = 'SELECT * FROM SEREB.Aln_Data\
					INNER JOIN SEREB.Alignment ON SEREB.Aln_Data.aln_id = SEREB.Alignment.Aln_id\
					INNER JOIN SEREB.Residues ON SEREB.Aln_Data.res_id = SEREB.Residues.resi_id\
					INNER JOIN (SELECT * from SEREB.Polymer_Data WHERE \
								SEREB.Polymer_Data.PData_id IN (SELECT PData_id from SEREB.Polymer_Alignments WHERE SEREB.Polymer_Alignments.Aln_id = '+str(aln_id)+')\
								AND \
								SEREB.Polymer_Data.strain_id IN (SELECT strain_id FROM SEREB.Species_TaxGroup WHERE taxgroup_id = '+str(parent1_id)+' OR taxgroup_id = '+str(parent2_id)+')) as filtered_polymers\
								ON SEREB.Residues.PolData_id = filtered_polymers.PData_id\
					INNER JOIN SEREB.Species ON filtered_polymers.strain_id = SEREB.Species.strain_id\
					WHERE SEREB.Alignment.aln_id = '+str(aln_id)
	with connection.cursor() as cursor:
		cursor.execute(SQLStatement)
		#print(dictfetchall(cursor)[0]['strain_id'])
		raw_result = dictfetchall(cursor)

	if len(raw_result) == 0:
		raise Http404("We do not have this combination of arguments in our database.")

	nogap_tupaln={}
	all_alnpositions=[]
	topgroup_name = ''
	for row in raw_result:
		all_alnpositions.append(row['aln_pos'])
		if (topgroup_name,row['strain']) in nogap_tupaln:
			nogap_tupaln[(topgroup_name,row['strain'])].append((row['unModResName'], row['aln_pos']))
		else:
			try:
				topgroup_query = Taxgroups.objects.filter(grouplevel=parent1_level, speciestaxgroup__strain=row['strain_id'])[0]
				topgroup_name = topgroup_query.groupname
			except:
				raise Http404("No superkingdom result for taxid"+row['strain_id']+"!")
			nogap_tupaln[(topgroup_name,row['strain'])]=[]
			nogap_tupaln[(topgroup_name,row['strain'])].append((row['unModResName'], row['aln_pos']))
	fasta_string=''

	for kingdom_strain in nogap_tupaln:
		strain = re.sub(' ','_',kingdom_strain[1])
		fasta_string+='\n>'+kingdom_strain[0]+'_'+strain+'\n'
		mem = 1
		for index, resi_pos in enumerate(nogap_tupaln[kingdom_strain], start=1):
			if mem == resi_pos[1]:
				mem = mem+1
			elif mem < resi_pos[1]:
				diff = resi_pos[1]-mem
				for i in range(0,diff):
					mem = mem+1
					fasta_string+='-'
				mem = mem+1
			else:
				raise ValueError("This shouldn't be possible!")
			fasta_string+=resi_pos[0]
			if index == len(nogap_tupaln[kingdom_strain]):
				if resi_pos[1] < max(all_alnpositions):
					diff = max(all_alnpositions)-resi_pos[1]
					for index2,i in enumerate(range(0,diff), start=1):
						fasta_string+='-'
						if index2 == diff:
							fasta_string+='\n'
				else:
					fasta_string+='\n'
	literal_string = re.sub(r'\n\n','\n',fasta_string,flags=re.M)
	return literal_string.lstrip().encode('unicode-escape').decode('ascii'),max(all_alnpositions)

def build_alignment(rawMYSQLresult):
	'''
	In here add a way to do the phase; we iterate over resis, so all we need to check is what phase the ones from PYRFU are.
	'''
	nogap_tupaln={}
	all_alnpositions=[]
	for row in rawMYSQLresult:
		all_alnpositions.append(row.aln_pos)
		if row.strain in nogap_tupaln:
			nogap_tupaln[row.strain].append((row.unModResName, row.aln_pos))
		else:
			nogap_tupaln[row.strain]=[]
			nogap_tupaln[row.strain].append((row.unModResName, row.aln_pos))
	fasta_string=''
	for strain in nogap_tupaln:
		strain1 = re.sub(' ','_',strain)
		fasta_string+='\n>'+strain1+'\n'
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
			else:
				raise ValueError("This shouldn't be possible!")
			fasta_string+=resi_pos[0]
			if index == len(nogap_tupaln[strain]):
				if resi_pos[1] < max(all_alnpositions):
					diff = max(all_alnpositions)-resi_pos[1]
					for index2,i in enumerate(range(0,diff), start=1):
						fasta_string+='-'
						if index2 == diff:
							fasta_string+='\n'
				else:
					fasta_string+='\n'
	literal_string = re.sub(r'\n\n','\n',fasta_string,flags=re.M)
	return literal_string.lstrip().encode('unicode-escape').decode('ascii'),max(all_alnpositions)

def pdbid_to_strainid(pdbid):
	'''Transforms PDBID to taxid'''
	structure_id = Threedstructures.objects.filter(structurename = pdbid)[0]
	secondary_id = SecondaryTertiary.objects.values("secondary_structure").filter(number_3d_structure = structure_id)[0]["secondary_structure"]
	anchor_taxid = Secondarystructures.objects.values("strain_fk").filter(secstr_id = secondary_id)[0]["strain_fk"]
	return anchor_taxid

def api_twc(request, align_name, tax_group1, tax_group2, anchor_structure):

	#### _____________Transform PDBID to taxid______________ ####
	anchor_taxid = pdbid_to_strainid(anchor_structure)

	#### This should be separate view with its own URL for serving multi-group alignments ####
	filter_strain = str(Species.objects.filter(strain_id = anchor_taxid)[0].strain).replace(" ", "_")

	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	fastastring,max_aln_length1 = sql_filtered_aln_query_two_parents(align_id,tax_group1,tax_group2)
	concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)
	
	#### _____________Trim down the alignment______________ ####
	from alignments.Shannon import species_index_to_aln_index, truncate_aln
	from Bio import AlignIO
	from io import StringIO
	alignment = list(AlignIO.parse(StringIO(concat_fasta), 'fasta'))[0]
	aln_anchor_map = species_index_to_aln_index(alignment, filter_strain)
	alignment = truncate_aln(alignment, list(aln_anchor_map[0].keys()), aln_anchor_map=aln_anchor_map[0])
	#### _______________Calculate TwinCons_________________ ####
	
	from TwinCons.bin import PhyMeas
	list_for_phymeas = ['-as',alignment.format("fasta"), '-r', '-bl']
	alnindex_score,sliced_alns,number_of_aligned_positions=PhyMeas.main(list_for_phymeas)
	list_for_topology_viewer = []
	for alnindex in alnindex_score:
		list_for_topology_viewer.append([alnindex,alnindex_score[alnindex][0]])
	return JsonResponse(list_for_topology_viewer, safe = False)

def twincons(request, align_name, tax_group1, tax_group2, anchor_structure):
	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	polymerid = PolymerData.objects.values("pdata_id").filter(polymeralignments__aln = align_id, strain = taxid)[0]["pdata_id"]
	chainid = Chainlist.objects.values("chainname").filter(polymer = polymerid)[0]["chainname"]
	context = {'pdbid': anchor_structure, 'chainid': chainid, 'entropy_address':align_name+"/"+str(tax_group)+"/"+str(anchor_structure)}
	return render(request, 'alignments/entropy_detail.html', context)

def entropy(request, align_name, tax_group, anchor_structure):
	from alignments import Shannon
	taxid = pdbid_to_strainid(anchor_structure)
	filter_strain = Species.objects.filter(strain_id = taxid)[0].strain
	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	polymerid = PolymerData.objects.values("pdata_id").filter(polymeralignments__aln = align_id, strain = taxid)[0]["pdata_id"]
	chainid = Chainlist.objects.values("chainname").filter(polymer = polymerid)[0]["chainname"]
	fastastring,max_aln_length = sql_filtered_aln_query(align_id,tax_group)
	aln_shannon_list = Shannon.main(['-a',fastastring,'-f','fastastring','--return_within','-s',filter_strain])
	#print(aln_shannon_list)
	context = {'pdbid': anchor_structure, 'chainid': chainid, 'shannon_dictionary': aln_shannon_list, 'entropy_address':align_name+"/"+str(tax_group)+"/"+str(anchor_structure)}
	return render(request, 'alignments/entropy_detail.html', context)

def api_entropy(request, align_name, tax_group, anchor_structure):
	from alignments import Shannon
	import os
	from django.conf import settings
	taxid = pdbid_to_strainid(anchor_structure)
	filter_strain = Species.objects.filter(strain_id = taxid)[0].strain
	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	fastastring,max_aln_length = sql_filtered_aln_query(align_id,tax_group)
	aln_shannon_list = Shannon.main(['-a',fastastring,'-f','fastastring','--return_within','-s',filter_strain])
	return JsonResponse(aln_shannon_list, safe = False)

def index(request):
	some_Alignments = Alignment.objects.all()
	superKingdoms = Taxgroups.objects.raw('SELECT * FROM SEREB.TaxGroups WHERE\
		 SEREB.TaxGroups.groupLevel = "superkingdom";')
	context = {
		'some_Alignments': some_Alignments,
		'superKingdoms': superKingdoms
	}
	return render(request, 'alignments/index.html', context)

def index_orthologs(request):
	three_d_structures = Threedstructures.objects.all()
	some_Alignments = Alignment.objects.all()
	superKingdoms = Taxgroups.objects.raw('SELECT * FROM SEREB.TaxGroups WHERE\
		 SEREB.TaxGroups.groupLevel = "superkingdom";')
	context = {
		'some_Alignments': some_Alignments,
		'superKingdoms': superKingdoms,
		'threeDstructures': three_d_structures
	}
	return render(request, 'alignments/index_orthologs.html', context)

def jalview(request):
	return render(request, "alignments/jalview.html")

def rProtein(request, align_name, tax_group):
	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	fastastring,max_aln_length = sql_filtered_aln_query(align_id,tax_group)
	context = {'fastastring': fastastring, 'aln_name':str(Alignment.objects.filter(aln_id = align_id)[0].name)}
	return render(request, 'alignments/detail.html', context)

def rRNA(request, align_name, tax_group):
	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	fastastring,max_aln_length = sql_filtered_aln_query(align_id,tax_group)
	context = {'fastastring': fastastring, 'aln_name':str(Alignment.objects.filter(aln_id = align_id)[0].name)}
	return render(request, 'alignments/rRNA.html', context)

def upload(request):
	context = {

	}
	return render(request, 'alignments/upload.html', context)

def submit(request):
	data_pairs = []

	context = {
		'data_pairs' : data_pairs
	}
	if request.method == 'POST' and 'filename' in request.FILES:
		file = request.FILES['filename']
		file_iterator = iter(file)
		# Skip the title line
		print(file_iterator.__next__().decode())
		while True:
			try:
				entry = file_iterator.__next__().decode().strip().split(',')
				data_pairs.append((int(entry[0]), float(entry[1])))
			except StopIteration:
				break
	return render(request, 'alignments/csvDisplay.html', context)

def submitAlignment(request):
	data_pairs = []
	context = {
		'data_pairs' : data_pairs
	}
	if request.method == 'POST' and 'filename' in request.FILES:
		file = request.FILES['filename']
		file_iterator = iter(file)
		while True:
			try:
				data_pairs.append((file_iterator.__next__().decode().strip(), file_iterator.__next__().decode().strip()))
			except StopIteration:
				break
	return render(request, 'alignments/alignmentsDisplay.html', context)
