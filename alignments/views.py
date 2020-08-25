import re

import urllib.request
import json

import datetime

import subprocess
from subprocess import Popen, PIPE

import os

from django.shortcuts import render
from django.http import HttpResponse, Http404, JsonResponse
from django.urls import reverse_lazy
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.views.generic import ListView, CreateView, UpdateView

from alignments.models import *
from alignments.taxonomy_views import *
from alignments.residue_api import *
from alignments.structure_api import *
from alignments.fold_api import *
from alignments.alignment_query_and_build import para_aln
import alignments.alignment_query_and_build as aqab


def trim_alignment(concat_fasta, filter_strain):
	'''Reads a fasta string into alignment and trims it down by filter sequence'''
	from alignments.Shannon import species_index_to_aln_index, truncate_aln
	from Bio import AlignIO
	from io import StringIO
	alignment = list(AlignIO.parse(StringIO(concat_fasta), 'fasta'))[0]
	aln_anchor_map, anchor_ix_in_alignment = species_index_to_aln_index(alignment, filter_strain)
	alignment = truncate_aln(alignment, list(aln_anchor_map.keys()), aln_anchor_map=aln_anchor_map)
	return alignment

def calculate_twincons(alignment):
	'''Calculates twincons score given an alignment object.
	Returns data in a list format for the topology viewer'''
	from TwinCons.bin import PhyMeas
	list_for_phymeas = ['-as',alignment.format("fasta"), '-r', '-bl']
	alnindex_score,sliced_alns,number_of_aligned_positions=PhyMeas.main(list_for_phymeas)
	list_for_topology_viewer = []
	for alnindex in alnindex_score:
		list_for_topology_viewer.append([alnindex,alnindex_score[alnindex][0]])
	return list_for_topology_viewer

def upload_custom_data_for_mapping(request):
	if request.method == 'POST' and 'filename' in request.FILES:
		data_pairs = []
		file = request.FILES['filename']
		file_iterator = iter(file)
		while True:
			try:
				entry = file_iterator.__next__().decode().strip().split(',')
				data_pairs.append((int(entry[0]), float(entry[1])))
			except StopIteration:
				break
		request.session['csv'] = data_pairs
	if request.method == 'GET':
		data_pairs = request.session.get('csv')
		return JsonResponse(data_pairs, safe = False)

def api_twc_with_upload(request, anchor_structure):
	#### _____________Transform PDBID to taxid______________ ####
	anchor_taxid = pdbid_to_strainid(anchor_structure)

	fastastring = request.session.get('fasta')
	#print('fastastring:\n' + fastastring)

	concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)
	#### _____________Trim down the alignment______________ ####
	alignment = trim_alignment(concat_fasta, str(anchor_taxid))

	#### _______________Calculate TwinCons_________________ ####
	list_for_topology_viewer = calculate_twincons(alignment)

	return JsonResponse(list_for_topology_viewer, safe = False)

def api_twc_parameterless(request):
	now = datetime.datetime.now()
	fileNameSuffix = "_" + str(now.year) + "_" + str(now.month) + "_" + str(now.day) + "_" + str(now.hour) + "_" + str(now.minute) + "_" + str(now.second) + "_" + str(now.microsecond)
	alignmentFileName = "./static/alignment" + fileNameSuffix + ".txt"
	ebiFileName = "./static/ebi_sequence" + fileNameSuffix + ".txt"
	mappingFileName = ebiFileName + ".map"

	fh = open(alignmentFileName, "w")
	fh.write(request.POST["fasta"])
	fh.close()

	fh = open(ebiFileName, "w")
	fh.write(">ebi_sequence\n")
	fh.write(request.POST["ebi_sequence"])
	fh.close()

	pipe = Popen("mafft --addfull " + ebiFileName + " --mapout " + alignmentFileName + "; cat " + mappingFileName, stdout=PIPE, shell=True)
	output = pipe.communicate()[0]
	text = output.decode("ascii")
	print(text)

	os.remove(alignmentFileName)
	os.remove(ebiFileName)
	os.remove(mappingFileName)

	context = {}
	return render(request, "alignments/dummyPage.html", context)

def api_twc(request, align_name, tax_group1, tax_group2, anchor_structure=''):

	#### _____________Transform PDBID to taxid______________ ####
	if anchor_structure != '':
		anchor_taxid = pdbid_to_strainid(anchor_structure)
		filter_strain = str(Species.objects.filter(strain_id = anchor_taxid)[0].strain).replace(" ", "_")

	fastastring = request.POST.get('fasta')
	if fastastring is None:
		#### _________Query database for the alignment__________ ####
		align_id = Alignment.objects.filter(name = align_name)[0].aln_id

		rawsqls = []
		for parent in [tax_group1, tax_group2]:
			rawsqls.append((aqab.sql_filtered_aln_query(align_id, parent), Taxgroups.objects.get(pk=parent).groupname))
		nogap_tupaln = dict()
		max_alnposition = 0
		for rawsql, parent in rawsqls:
			nogap_tupaln, max_alnposition= aqab.query_to_dict_structure(rawsql, parent, nogap_tupaln, max_alnposition)

		print(nogap_tupaln)
		#### __________________Build alignment__________________ ####
		fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_alnposition)
	
	concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)
	#print(concat_fasta)
	
	#### ______________Trim down the alignment______________ ####
	if anchor_structure != '':
		alignment = trim_alignment(concat_fasta, filter_strain)
	else:
		alignment = concat_fasta
	
	#### ________________Calculate TwinCons_________________ ####
	list_for_topology_viewer = calculate_twincons(alignment)
	
	return JsonResponse(list_for_topology_viewer, safe = False)

def minmaxIndex_handler(minIndex, maxIndex):
	if (minIndex == ''):
		minIndex = str(0)
	else:
		minIndex = str(minIndex)
	if (maxIndex == ''):
		maxIndex = str(100000)
	else:
		maxIndex = str(maxIndex)
	return minIndex, maxIndex

def twincons_handler(request, anchor_structure, chain, align_name='', tax_group1='', tax_group2='', minIndex = '', maxIndex = ''):
	from django.urls import resolve
	current_url = resolve(request.path_info).url_name
	minIndex, maxIndex = minmaxIndex_handler(minIndex, maxIndex)
	context = dict()
	context = {
		'pdbid': anchor_structure, 
		'chainid': chain,
		'minIndex' : minIndex,
		'maxIndex' : maxIndex
	}
	if current_url == 'twc_with_upload':
		context['entropy_address'] = "upload/twc-api/"+str(anchor_structure)
	elif current_url == 'custom_csv_data_viewer':
		upload_custom_data_for_mapping(request)
		context['entropy_address'] = "custom-csv-data"
	elif current_url == 'twincons':
		context['entropy_address'] = "twc-api/"+align_name+"/"+str(tax_group1)+"/"+str(tax_group2)+"/"+str(anchor_structure)
	return render(request, 'alignments/twc_detail.html', context)

def entropy(request, align_name, tax_group, anchor_structure):
	from alignments import Shannon
	taxid = pdbid_to_strainid(anchor_structure)
	filter_strain = Species.objects.filter(strain_id = taxid)[0].strain
	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	polymerid = PolymerData.objects.values("pdata_id").filter(polymeralignments__aln = align_id, strain = taxid)[0]["pdata_id"]
	chainid = Chainlist.objects.values("chainname").filter(polymer = polymerid)[0]["chainname"]
	
	rawsql_result = aqab.sql_filtered_aln_query(align_id, tax_group)
	nogap_tupaln, max_aln_length = aqab.query_to_dict_structure(rawsql_result, Taxgroups.objects.get(pk=tax_group).groupname)
	fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_aln_length)
	#fastastring,max_aln_length = aqab.sql_filtered_aln_query(align_id,tax_group)
	aln_shannon_list = Shannon.main(['-a',fastastring,'-f','fastastring','--return_within','-s',filter_strain])
	#print(aln_shannon_list)
	context = {
		'pdbid': anchor_structure, 
		'chainid': chainid, 
		'shannon_dictionary': aln_shannon_list, 
		'entropy_address':"entropy-api/"+align_name+"/"+str(tax_group)+"/"+str(anchor_structure)
	}
	return render(request, 'alignments/entropy_detail.html', context)

def api_entropy(request, align_name, tax_group, anchor_structure):
	from alignments import Shannon
	import os
	from django.conf import settings
	taxid = pdbid_to_strainid(anchor_structure)
	filter_strain = Species.objects.filter(strain_id = taxid)[0].strain
	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	
	rawsql_result = aqab.sql_filtered_aln_query(align_id, tax_group)
	nogap_tupaln, max_aln_length = aqab.query_to_dict_structure(rawsql_result, Taxgroups.objects.get(pk=tax_group).groupname)
	fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_aln_length)
	#fastastring,max_aln_length = aqab.sql_filtered_aln_query(align_id,tax_group)
	aln_shannon_list = Shannon.main(['-a',fastastring,'-f','fastastring','--return_within','-s',filter_strain])
	return JsonResponse(aln_shannon_list, safe = False)

def index(request):
	some_Alignments = Alignment.objects.all()
	superKingdoms = Taxgroups.objects.raw('SELECT * FROM SEREB.TaxGroups WHERE\
		 SEREB.TaxGroups.groupLevel = "superkingdom";')
	
	context = {
		'props': list(Taxgroups.objects.values('taxgroup_id', 'groupname')),
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
		'threeDstructures': three_d_structures,
		'file_name' : None
	}
	return render(request, 'alignments/index_orthologs.html', context)

# def visualizerHelper(request, urlSuffix):
# 	data = json.load(urllib.request.urlopen("http://127.0.0.1:8000/orthologs/twc-api/" + urlSuffix))
# 	xyPairs = []
# 	context = {
# 		"xyPairs" : xyPairs
# 	}
# 	for xyPair in data:
# 		xyPairs.append(xyPair)
# 	return render(request, 'alignments/simpleVisualization.html', context)

def visualizer(request, align_name, tax_group1, tax_group2, anchor_structure = ''):
	twc_api_url = "http://127.0.0.1:8000/orthologs/twc-api/" + align_name + "/" + str(tax_group1) + "/" + str(tax_group2) + "/" + anchor_structure
	context = {
		"twc_api_url" : twc_api_url
	}
	return render(request, 'alignments/simpleVisualization.html', context)
	#return visualizerHelper(request, align_name + "/" + str(tax_group1) + "/" + str(tax_group2) + "/" + anchor_structure)

def upload_custom_data(request):
	return render(request, 'alignments/upload_custom_data.html')

def paralog_entry_form(request):
	return render(request, 'alignments/index_paralogs.html')

def paralog_display_entropy(request, align_name, fold1, fold2):
	pass
	#return render(request, 'alignments/twc_detail.html', context)



def simple_fasta(request, aln_id, tax_group, internal=False):
	rawsqls = []
	if type(tax_group) == int:
		tax_group = str(tax_group)
	for parent in tax_group.split(','):
		rawsqls.append((aqab.sql_filtered_aln_query(aln_id, parent), Taxgroups.objects.get(pk=parent).groupname))

	nogap_tupaln = dict()
	max_alnposition = 0

	for rawsql, parent in rawsqls:
		nogap_tupaln, max_alnposition= aqab.query_to_dict_structure(rawsql, parent, nogap_tupaln, max_alnposition)
	
	fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_alnposition)
	unf_species_list = [x.split('\\')[0] for x in fastastring.split('>')[1:]]
	unf_seq_list = [x.split('\\n')[1] for x in fastastring.split('>')[1:]]
	list_for_intersect = list()
	for sequence in unf_seq_list:
		iterator = re.finditer('-', sequence)
		gap_positions = [m.start(0) for m in iterator]
		list_for_intersect.append(gap_positions)
	gap_only_cols = list(set(list_for_intersect[0]).intersection(*list_for_intersect))
	filtered_spec_list = [re.sub('_',' ', re.sub(r'^.*?_', '', x)) for x in unf_species_list]
	if internal:
		return fastastring
	concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)
	return JsonResponse([concat_fasta,filtered_spec_list,gap_only_cols,frequency_list], safe = False)

# def simple_fasta(request, aln_id, tax_group, internal=False):
# 	rawsql_result = aqab.sql_filtered_aln_query(aln_id, tax_group)
# 	nogap_tupaln = dict()
# 	nogap_tupaln, max_aln_length = aqab.query_to_dict_structure(rawsql_result, Taxgroups.objects.get(pk=tax_group).groupname, nogap_tupaln)
# 	fastastring = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_aln_length)
# 	if internal:
# 		return fastastring
# 	concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)
# 	return JsonResponse(concat_fasta, safe = False)

def rProtein(request, align_name, tax_group):
	#if tax_group == 0 - no filter
	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	fastastring = simple_fasta(request, align_id, tax_group, internal=True)
	print(fastastring)
	#fastastring,max_aln_length = aqab.sql_filtered_aln_query(align_id,tax_group)
	context = {'fastastring': fastastring, 'aln_name':str(Alignment.objects.filter(aln_id = align_id)[0].name)}
	return render(request, 'alignments/detail.html', context)

def rRNA(request, align_name, tax_group):
	from Bio.SeqUtils import IUPACData
	align_id = Alignment.objects.filter(name = align_name)[0].aln_id
	rawsql_result = aqab.sql_filtered_aln_query(align_id, tax_group)
	nogap_tupaln = dict()
	nogap_tupaln, max_aln_length = aqab.query_to_dict_structure(rawsql_result, Taxgroups.objects.get(pk=tax_group).groupname, nogap_tupaln)
	fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_aln_length, IUPACData.unambiguous_rna_letters)
	#fastastring,max_aln_length = aqab.sql_filtered_aln_query(align_id,tax_group)
	context = {'fastastring': fastastring, 'aln_name':str(Alignment.objects.filter(aln_id = align_id)[0].name)}
	return render(request, 'alignments/rRNA.html', context)