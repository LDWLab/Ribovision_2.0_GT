import re

import urllib.request
import json

import datetime

import subprocess
from subprocess import Popen, PIPE

import os

from django.shortcuts import render
from django.http import HttpResponse, Http404, JsonResponse, HttpResponseServerError
from django.urls import reverse_lazy
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.views.generic import ListView, CreateView, UpdateView
from django.views.decorators.csrf import csrf_exempt

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
	from TwinCons.bin import TwinCons
	list_for_phymeas = ['-as',alignment.format("fasta"), '-r', '-mx', 'blosum62']
	alnindex_score, sliced_alns, number_of_aligned_positions, gp_mapping = TwinCons.main(list_for_phymeas)
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

def constructEbiAlignmentString(fasta, ebi_sequence, startIndex):
	now = datetime.datetime.now()
	fileNameSuffix = "_" + str(now.year) + "_" + str(now.month) + "_" + str(now.day) + "_" + str(now.hour) + "_" + str(now.minute) + "_" + str(now.second) + "_" + str(now.microsecond)
	alignmentFileName = "./static/alignment" + fileNameSuffix + ".txt"
	ebiFileName = "./static/ebi_sequence" + fileNameSuffix + ".txt"
	mappingFileName = ebiFileName + ".map"

	fh = open(alignmentFileName, "w")
	fh.write(fasta)
	fh.close()

	fh = open(ebiFileName, "w")
	fh.write(">ebi_sequence\n")
	fh.write(ebi_sequence)
	fh.close()

	shiftIndexBy = 0
	if startIndex > 1:
		shiftIndexBy = startIndex - 1

	pipe = Popen("mafft --addfull " + ebiFileName + " --mapout " + alignmentFileName + "; cat " + mappingFileName, stdout=PIPE, shell=True)
	output = pipe.communicate()[0]
	text = output.decode("ascii").split('\n#')[1]

	mapping = dict()
	firstLine = True
	for line in text.split('\n'):
		if firstLine:
			firstLine = False
			continue
		row = line.split(', ')
		if len(row) < 3:
			continue
		if row[2] == '-':
			continue
		if row[1] == '-':
			raise Http404("Mapping did not work properly.")
		mapping[int(row[2])] = int(row[1]) + shiftIndexBy

	for removeFile in [alignmentFileName, ebiFileName, mappingFileName]:
		os.remove(removeFile)
	return mapping

def request_post_data(post_data):
	fasta = post_data["fasta"]
	ebi_sequence = post_data["ebi_sequence"]
	startIndex = int(post_data["startIndex"])
	return fasta, ebi_sequence, startIndex

def make_map_from_alnix_to_sequenceix(request):
	fasta, ebi_sequence, startIndex = request_post_data(request.POST)
	mapping = constructEbiAlignmentString(fasta, ebi_sequence, startIndex)
	return JsonResponse(mapping, safe = False)

def api_twc_parameterless(request):
	fasta, ebi_sequence, startIndex = request_post_data(request.POST)

	mapping = constructEbiAlignmentString(fasta, ebi_sequence, startIndex)
	concat_fasta = re.sub(r'\\n','\n', fasta,flags=re.M)
	list_for_topology_viewer = calculate_twincons(concat_fasta)
	return JsonResponse([list_for_topology_viewer, mapping], safe = False)

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

def propensity_sample_fasta(request):
    pass

def propensity_data(request):
    test = [['eS04c_Ferroplasma_acidiphi', 0.0, 0.04081632653061224, 0.12244897959183673, 0.02040816326530612, 0.08163265306122448, 0.1836734693877551, 0.22448979591836735, 0.32653061224489793], ['eS04c_Picrophilus_torridus', 0.0, 0.04081632653061224, 0.12244897959183673, 0.02040816326530612, 0.09183673469387756, 0.1836734693877551, 0.23469387755102042, 0.30612244897959184], ['eS04c_Thermoplasma_volcani', 0.0, 0.034013605442176874, 0.10884353741496598, 0.034013605442176874, 0.09523809523809523, 0.1836734693877551, 0.2108843537414966, 0.3333333333333333], ['eS4e_Ashbya_gossypii', 0.0, 0.04081632653061224, 0.10204081632653061, 0.04081632653061224, 0.10204081632653061, 0.17346938775510204, 0.21428571428571427, 0.32653061224489793], ['eS4e_Yarrowia_lipolytica_', 0.0, 0.044897959183673466, 0.09795918367346938, 0.044897959183673466, 0.10612244897959183, 0.15918367346938775, 0.2163265306122449, 0.3306122448979592], ['eS4e_Saccharomyces_cerevi', 0.0, 0.047619047619047616, 0.09523809523809523, 0.047619047619047616, 0.10884353741496598, 0.1564625850340136, 0.21768707482993196, 0.32653061224489793], ['eS4e_Schizosaccharomyces_', 0.0029154518950437317, 0.04956268221574344, 0.09037900874635568, 0.04956268221574344, 0.11078717201166181, 0.15160349854227406, 0.21865889212827988, 0.32653061224489793], ['eS4e_Homo_sapiens', 0.00510204081632653, 0.04846938775510204, 0.08673469387755102, 0.05102040816326531, 0.11224489795918367, 0.14795918367346939, 0.2193877551020408, 0.32908163265306123], ['eS4e_Pan_troglodytes', 0.006802721088435374, 0.047619047619047616, 0.08390022675736962, 0.05215419501133787, 0.11337868480725624, 0.1473922902494331, 0.2199546485260771, 0.3287981859410431], ['eS4e_Gallus_gallus', 0.00816326530612245, 0.04897959183673469, 0.07959183673469387, 0.053061224489795916, 0.11632653061224489, 0.14489795918367346, 0.22040816326530613, 0.32857142857142857], ['eS4e_Latimeria_chalumnae', 0.00927643784786642, 0.05009276437847866, 0.07606679035250463, 0.05380333951762523, 0.11688311688311688, 0.14471243042671614, 0.22077922077922077, 0.32838589981447125], ['eS4e_Monodelphis_domestic', 0.01020408163265306, 0.05102040816326531, 0.07142857142857142, 0.05442176870748299, 0.11904761904761904, 0.1445578231292517, 0.22108843537414966, 0.3282312925170068], ['eS4e_Mus_musculus', 0.01098901098901099, 0.05180533751962323, 0.06907378335949764, 0.054945054945054944, 0.12087912087912088, 0.14285714285714285, 0.22135007849293564, 0.3281004709576138], ['eS4e_Rattus_norvegicus', 0.011661807580174927, 0.052478134110787174, 0.06705539358600583, 0.05539358600583091, 0.12244897959183673, 0.141399416909621, 0.22157434402332363, 0.32798833819241985], ['eS4e_Xenopus_laevis', 0.012244897959183673, 0.053061224489795916, 0.0653061224489796, 0.055782312925170066, 0.12380952380952381, 0.1414965986394558, 0.2217687074829932, 0.32653061224489793], ['eS4e_Danio_rerio', 0.012755102040816327, 0.05229591836734694, 0.06377551020408163, 0.05612244897959184, 0.125, 0.13903061224489796, 0.22321428571428573, 0.3278061224489796], ['eS4e_Drosophila_melanogas', 0.013205282112845138, 0.05282112845138055, 0.06362545018007203, 0.056422569027611044, 0.12484993997599039, 0.13925570228091236, 0.22208883553421369, 0.3277310924369748], ['eS4e_Aedes_albopictus', 0.012471655328798186, 0.05442176870748299, 0.06235827664399093, 0.05782312925170068, 0.12471655328798185, 0.13945578231292516, 0.2222222222222222, 0.32653061224489793], ['eS4e_Cyanidioschyzon_mero', 0.01288936627282492, 0.05477980665950591, 0.06337271750805586, 0.05800214822771214, 0.12567132116004295, 0.13748657357679914, 0.22234156820622986, 0.32545649838882923], ['eS4e_Caenorhabditis_brigg', 0.013265306122448979, 0.05510204081632653, 0.06224489795918367, 0.058163265306122446, 0.12448979591836734, 0.14081632653061224, 0.22040816326530613, 0.32551020408163267], ['eS4e_Caenorhabditis_elega', 0.013605442176870748, 0.05539358600583091, 0.061224489795918366, 0.05830903790087463, 0.12342079689018465, 0.14480077745383868, 0.21865889212827988, 0.32458697764820216], ['eS4e_Dictyostelium_discoi', 0.012987012987012988, 0.055658627087198514, 0.061224489795918366, 0.05844155844155844, 0.1261595547309833, 0.1456400742115028, 0.21706864564007422, 0.3228200371057514], ['eS4e_Oryza_sativa_Japonic', 0.013297872340425532, 0.057624113475177305, 0.061170212765957445, 0.05851063829787234, 0.12411347517730496, 0.15070921985815602, 0.2154255319148936, 0.3191489361702128], ['eS4e_Arabidopsis_thaliana', 0.013593882752761258, 0.059473237043330504, 0.06117247238742566, 0.059473237043330504, 0.12489379779099405, 0.15123194562446898, 0.2141036533559898, 0.3160577740016992], ['eS4e_Thalassiosira_pseudo', 0.013866231647634585, 0.05954323001631321, 0.06035889070146819, 0.05954323001631321, 0.12561174551386622, 0.15089722675367048, 0.2128874388254486, 0.3172920065252855]]
    return JsonResponse(test, safe = False)

def propensities(request):
    # how to not hardcode this?
    propensity_data = "http://127.0.0.1:8000/propensity-data"

    # where does the context variable come from? what does it do?
    context = {"propensity_data" : propensity_data}
    
    return render(request, 'alignments/propensities.html', context)