from django.shortcuts import render
from django.http import HttpResponse
from alignments.models import *
from django.urls import reverse_lazy
from django.views.generic import ListView, CreateView, UpdateView
import re

from django.http import JsonResponse

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

def buildTaxonomy(request):
	# nodeSet = []
	# for taxgroup in taxgroups:
	# 	nodeSet.append({'label' : taxgroup.groupname})
	# tree = {
	# 	'label' : 'Select a taxgroup:',
	# 	'nodes' : nodeSet
	# }

	# tree = {
	# 'label': 'root',
	# 'nodes': [
	# {
	# 	'label': 'item1',
	# 	'nodes': [
	# 	{
	# 		'label': 'item1.1'
	# 	},{
	# 		'label': 'item1.2',
	# 		'nodes': [
	# 		{
	# 			'label': 'item1.2.1'
	# 		}]
	# 	}]
	# },{
	# 	'label': 'item2'
	# }]}

	taxgroups = Taxgroups.objects.raw('SELECT * FROM SEREB.TaxGroups WHERE\
		 SEREB.TaxGroups.groupLevel = "superkingdom";')
	taxonomy = []
	for taxgroup in taxgroups:
		subtaxonomy = {
			'label' : taxgroup.groupname,
			'nodes' : buildTaxonomyRecurse(taxgroup.parent)
		}
		print(str(subtaxonomy['nodes']))
		taxonomy.append(subtaxonomy)
	tree = {
		'label' : 'Select a taxgroup:',
		'nodes' : taxonomy
	}
	return JsonResponse(tree, safe = False)

def buildTaxonomyRecurse(parentIndex):
	mySQLStr = 'SELECT * FROM SEREB.TaxGroups WHERE\
		 SEREB.TaxGroups.groupLevel = "' + str(parentIndex) + '";'
	taxgroups = Taxgroups.objects.raw(mySQLStr)
	print(mySQLStr)
	taxonomy = []
	for taxgroup in taxgroups:
		subtaxonomy = {
			'label' : taxgroup.groupname,
			'nodes' : {}#buildTaxonomyRecurse()
		}
		taxonomy.append(subtaxonomy)
	return taxonomy

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

def index(request):
	some_Alignments = Alignment.objects.all()
	superKingdoms = Taxgroups.objects.raw('SELECT * FROM SEREB.TaxGroups WHERE\
		 SEREB.TaxGroups.groupLevel = "superkingdom";')
	context = {
		'some_Alignments': some_Alignments,
		'superKingdoms': superKingdoms
	}
	return render(request, 'alignments/index.html', context)

def detail(request, align_name):
	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	fastastring,max_aln_length = sql_alignment_query(align_id)
	#print(fastastring)
	context = {'fastastring': fastastring, 'aln_name':str(Alignment.objects.filter(aln_id = align_id)[0].name)}
	return render(request, 'alignments/detail.html', context)

def rRNA(request, name):
	align_id = Alignment.objects.filter(name = name)[0].aln_id
	fastastring,max_aln_length = sql_alignment_query(align_id)
	context = {'fastastring': fastastring, 'aln_name':str(Alignment.objects.filter(aln_id = align_id)[0].name)}
	return render(request, 'alignments/rRNA.html', context)

class TaxgroupListView(ListView):
	model = Taxgroups
	context_object_name = 'taxgroups'

class TaxgroupCreateView(CreateView):
	model = Taxgroups
	fields = ('superkingdom', 'phyla', 'alignment')
	success_url = reverse_lazy('taxgroup_changelist')

class TaxgroupUpdateView(UpdateView):
	model = Taxgroups
	fields = ('superkingdom', 'phyla', 'alignment')
	success_url = reverse_lazy('taxgroup_changelist')
