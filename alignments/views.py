import re

from django.shortcuts import render
from django.http import HttpResponse, Http404, JsonResponse
from django.urls import reverse_lazy
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
	alnpos=[]
	if len(alnposition) == 0:
		raise Http404("We do not have this combination of arguments in our database.")
	fastastring,max_aln_length = build_alignment(alnposition)
	return fastastring,max_aln_length

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
	literal_string = re.sub(r'\n\n','\n',fasta_string)
	return literal_string.lstrip().encode('unicode-escape').decode('ascii'),max(all_alnpositions)

def api_twc(request, align_name, tax_group1, tax_group2, anchor_taxid):
	filter_strain = Species.objects.filter(strain_id = anchor_taxid)[0].strain
	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	fastastring1,max_aln_length1 = sql_filtered_aln_query(align_id,tax_group1)
	fastastring2,max_aln_length2 = sql_filtered_aln_query(align_id,tax_group2)
	fastastring1 = re.sub('>','>a_', fastastring1)
	fastastring2 = re.sub('>','>b_', fastastring2)
	concat_fasta = re.sub(r'\n>', '\\n' ,fastastring1+fastastring2)
	#concat_fasta = re.sub(r'\n\n', '\n' ,concat_fasta)
	return HttpResponse(concat_fasta, content_type="text/plain")

def entropy(request, align_name, tax_group, taxid):
	from alignments import Shannon
	filter_strain = Species.objects.filter(strain_id = taxid)[0].strain
	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	fastastring,max_aln_length = sql_filtered_aln_query(align_id,tax_group)
	aln_shannon_list = Shannon.main(['-a',fastastring,'-f','fastastring','--return_within','-s',filter_strain])
	#print(aln_shannon_list)
	context = {'shannon_dictionary': aln_shannon_list, 'entropy_address':align_name+"/"+str(tax_group)+"/"+str(taxid)}
	return render(request, 'alignments/entropy_detail.html', context)

def api_entropy(request, align_name, tax_group, taxid):
	from alignments import Shannon
	import os
	from django.conf import settings
	#file_ = open(os.path.join(settings.BASE_DIR, '1cbs_outliers.txt'))
	#string = file_.read()
	#return JsonResponse(string, safe = False)
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
	some_Alignments = Alignment.objects.all()
	superKingdoms = Taxgroups.objects.raw('SELECT * FROM SEREB.TaxGroups WHERE\
		 SEREB.TaxGroups.groupLevel = "superkingdom";')
	context = {
		'some_Alignments': some_Alignments,
		'superKingdoms': superKingdoms
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

