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
    with open('sample_fasta.fas') as sample:
        test = sample.read()
    return JsonResponse(test, safe = False)

def propensity_data(request):
    test = [['eS04c_Ferroplasma_acidiphi', 0.0, 0.04081632653061224, 0.12244897959183673, 0.02040816326530612, 0.08163265306122448, 0.1836734693877551, 0.22448979591836735, 0.32653061224489793], ['eS04c_Picrophilus_torridus', 0.0, 0.04081632653061224, 0.12244897959183673, 0.02040816326530612, 0.09183673469387756, 0.1836734693877551, 0.23469387755102042, 0.30612244897959184], ['eS04c_Thermoplasma_volcani', 0.0, 0.034013605442176874, 0.10884353741496598, 0.034013605442176874, 0.09523809523809523, 0.1836734693877551, 0.2108843537414966, 0.3333333333333333], ['eS4e_Ashbya_gossypii', 0.0, 0.04081632653061224, 0.10204081632653061, 0.04081632653061224, 0.10204081632653061, 0.17346938775510204, 0.21428571428571427, 0.32653061224489793], ['eS4e_Yarrowia_lipolytica_', 0.0, 0.044897959183673466, 0.09795918367346938, 0.044897959183673466, 0.10612244897959183, 0.15918367346938775, 0.2163265306122449, 0.3306122448979592], ['eS4e_Saccharomyces_cerevi', 0.0, 0.047619047619047616, 0.09523809523809523, 0.047619047619047616, 0.10884353741496598, 0.1564625850340136, 0.21768707482993196, 0.32653061224489793], ['eS4e_Schizosaccharomyces_', 0.0029154518950437317, 0.04956268221574344, 0.09037900874635568, 0.04956268221574344, 0.11078717201166181, 0.15160349854227406, 0.21865889212827988, 0.32653061224489793], ['eS4e_Homo_sapiens', 0.00510204081632653, 0.04846938775510204, 0.08673469387755102, 0.05102040816326531, 0.11224489795918367, 0.14795918367346939, 0.2193877551020408, 0.32908163265306123], ['eS4e_Pan_troglodytes', 0.006802721088435374, 0.047619047619047616, 0.08390022675736962, 0.05215419501133787, 0.11337868480725624, 0.1473922902494331, 0.2199546485260771, 0.3287981859410431], ['eS4e_Gallus_gallus', 0.00816326530612245, 0.04897959183673469, 0.07959183673469387, 0.053061224489795916, 0.11632653061224489, 0.14489795918367346, 0.22040816326530613, 0.32857142857142857], ['eS4e_Latimeria_chalumnae', 0.00927643784786642, 0.05009276437847866, 0.07606679035250463, 0.05380333951762523, 0.11688311688311688, 0.14471243042671614, 0.22077922077922077, 0.32838589981447125], ['eS4e_Monodelphis_domestic', 0.01020408163265306, 0.05102040816326531, 0.07142857142857142, 0.05442176870748299, 0.11904761904761904, 0.1445578231292517, 0.22108843537414966, 0.3282312925170068], ['eS4e_Mus_musculus', 0.01098901098901099, 0.05180533751962323, 0.06907378335949764, 0.054945054945054944, 0.12087912087912088, 0.14285714285714285, 0.22135007849293564, 0.3281004709576138], ['eS4e_Rattus_norvegicus', 0.011661807580174927, 0.052478134110787174, 0.06705539358600583, 0.05539358600583091, 0.12244897959183673, 0.141399416909621, 0.22157434402332363, 0.32798833819241985], ['eS4e_Xenopus_laevis', 0.012244897959183673, 0.053061224489795916, 0.0653061224489796, 0.055782312925170066, 0.12380952380952381, 0.1414965986394558, 0.2217687074829932, 0.32653061224489793], ['eS4e_Danio_rerio', 0.012755102040816327, 0.05229591836734694, 0.06377551020408163, 0.05612244897959184, 0.125, 0.13903061224489796, 0.22321428571428573, 0.3278061224489796], ['eS4e_Drosophila_melanogas', 0.013205282112845138, 0.05282112845138055, 0.06362545018007203, 0.056422569027611044, 0.12484993997599039, 0.13925570228091236, 0.22208883553421369, 0.3277310924369748], ['eS4e_Aedes_albopictus', 0.012471655328798186, 0.05442176870748299, 0.06235827664399093, 0.05782312925170068, 0.12471655328798185, 0.13945578231292516, 0.2222222222222222, 0.32653061224489793], ['eS4e_Cyanidioschyzon_mero', 0.01288936627282492, 0.05477980665950591, 0.06337271750805586, 0.05800214822771214, 0.12567132116004295, 0.13748657357679914, 0.22234156820622986, 0.32545649838882923], ['eS4e_Caenorhabditis_brigg', 0.013265306122448979, 0.05510204081632653, 0.06224489795918367, 0.058163265306122446, 0.12448979591836734, 0.14081632653061224, 0.22040816326530613, 0.32551020408163267], ['eS4e_Caenorhabditis_elega', 0.013605442176870748, 0.05539358600583091, 0.061224489795918366, 0.05830903790087463, 0.12342079689018465, 0.14480077745383868, 0.21865889212827988, 0.32458697764820216], ['eS4e_Dictyostelium_discoi', 0.012987012987012988, 0.055658627087198514, 0.061224489795918366, 0.05844155844155844, 0.1261595547309833, 0.1456400742115028, 0.21706864564007422, 0.3228200371057514], ['eS4e_Oryza_sativa_Japonic', 0.013297872340425532, 0.057624113475177305, 0.061170212765957445, 0.05851063829787234, 0.12411347517730496, 0.15070921985815602, 0.2154255319148936, 0.3191489361702128], ['eS4e_Arabidopsis_thaliana', 0.013593882752761258, 0.059473237043330504, 0.06117247238742566, 0.059473237043330504, 0.12489379779099405, 0.15123194562446898, 0.2141036533559898, 0.3160577740016992], ['eS4e_Thalassiosira_pseudo', 0.013866231647634585, 0.05954323001631321, 0.06035889070146819, 0.05954323001631321, 0.12561174551386622, 0.15089722675367048, 0.2128874388254486, 0.3172920065252855], ['eS4e_Plasmodium_falciparu', 0.01411764705882353, 0.058823529411764705, 0.0596078431372549, 0.058823529411764705, 0.12627450980392158, 0.14980392156862746, 0.21490196078431373, 0.3176470588235294], ['eS4e_Tetrahymena_thermoph', 0.013595166163141994, 0.05966767371601209, 0.05891238670694864, 0.05966767371601209, 0.12613293051359517, 0.15030211480362538, 0.21374622356495468, 0.31797583081570996], ['eS4e_Guillardia_theta', 0.013838310269482883, 0.05972323379461034, 0.05899490167516387, 0.05899490167516387, 0.12454479242534595, 0.15149308084486526, 0.21267297887836853, 0.3197378004369993], ['eS4e_Trypanosoma_brucei_b', 0.013361462728551337, 0.05977496483825598, 0.05977496483825598, 0.05907172995780591, 0.12376933895921238, 0.15260196905766527, 0.21237693389592124, 0.31926863572433195], ['eS4e_Leishmania_brazilien', 0.013596193065941536, 0.05982324949014276, 0.05982324949014276, 0.05914343983684568, 0.12236573759347383, 0.15363698164513936, 0.2107409925220938, 0.32087015635622024], ['eS04c_Theionarchaea_archae', 0.012919896640826873, 0.06136950904392765, 0.060723514211886306, 0.056847545219638244, 0.12467700258397933, 0.1511627906976744, 0.21382428940568476, 0.3184754521963824], ['eS04c_Nitrosopumilus_marit', 0.013157894736842105, 0.06015037593984962, 0.06140350877192982, 0.05576441102756892, 0.12343358395989974, 0.15538847117794485, 0.2136591478696742, 0.31704260651629074], ['eS04c_Cenarchaeum_symbiosu', 0.013990267639902677, 0.05900243309002433, 0.06265206812652069, 0.05474452554744526, 0.12347931873479319, 0.15815085158150852, 0.21289537712895376, 0.3150851581508516], ['eS04c_Haloredivivus_sp.', 0.013609467455621301, 0.06035502958579882, 0.06331360946745562, 0.05325443786982249, 0.12603550295857988, 0.16153846153846155, 0.21124260355029587, 0.3106508875739645], ['eS04c_Hadesarchaea_archaeo', 0.013801035077630822, 0.060379528464634846, 0.06382978723404255, 0.05290396779758482, 0.12535940195514664, 0.16101207590569291, 0.21161587119033928, 0.3110983323749281], ['eS04c_Methanothermobacter_', 0.013989927252378288, 0.0604364857302742, 0.06491326245103525, 0.05204252937884723, 0.12534974818130945, 0.16228315612758815, 0.21096810296586457, 0.3100167879127029], ['eS04c_Methanosphaera_stadt', 0.013623978201634877, 0.05994550408719346, 0.0664850136239782, 0.05122615803814714, 0.12534059945504086, 0.16348773841961853, 0.2098092643051771, 0.31008174386920984], ['eS04c_Lokiarchaeum_sp._GC1', 0.013742071881606765, 0.05919661733615222, 0.0660676532769556, 0.051268498942917545, 0.12526427061310783, 0.16543340380549684, 0.20983086680761098, 0.3091966173361522], ['eS04c_Korarchaeum_cryptofi', 0.013395157135497167, 0.05821741370427615, 0.0674909840288511, 0.05100463678516229, 0.125193199381762, 0.16537867078825347, 0.21020092735703247, 0.3091190108191654], ['eS04c_Methanopyrus_kandler', 0.013513513513513514, 0.05805805805805806, 0.06756756756756757, 0.05105105105105105, 0.12712712712712712, 0.16516516516516516, 0.2087087087087087, 0.3088088088088088], ['eS04c_Bathyarchaeota_archa', 0.013625304136253041, 0.057907542579075426, 0.06715328467153285, 0.05060827250608273, 0.12846715328467154, 0.16642335766423358, 0.2077858880778589, 0.30802919708029197], ['eS04c_Thermococcus', 0.013289036544850499, 0.05790223065970574, 0.0664451827242525, 0.05030849549121975, 0.12766967252017086, 0.16753678215472234, 0.20787850023730423, 0.3089700996677741], ['eS04c_Pyrococcus_furiosus', 0.012968967114404817, 0.058823529411764705, 0.06577119036591014, 0.05002315886984715, 0.12691060676239, 0.1681333950903196, 0.2079666512274201, 0.3094025011579435], ['eS4_PF', 0.012663952962460425, 0.05970149253731343, 0.06512890094979647, 0.04975124378109453, 0.12618724559023067, 0.16870194482134782, 0.20805065581184984, 0.3098145635459068], ['eS04c_Staphylothermus_mari', 0.012345679012345678, 0.059082892416225746, 0.06525573192239859, 0.04982363315696649, 0.12610229276895943, 0.17019400352733685, 0.20767195767195767, 0.30952380952380953], ['eS04c_Hyperthermus_butylic', 0.012037833190025795, 0.058469475494411005, 0.06577815993121239, 0.05030094582975064, 0.12596732588134135, 0.1706792777300086, 0.20765262252794497, 0.30911435941530524], ['eS04c_Aeropyrum_pernix_K1', 0.01174989509022241, 0.057490558120016785, 0.0658833403273185, 0.05077633235417541, 0.12673101133025597, 0.17079311791859, 0.2073017205203525, 0.3092740243390684], ['eS04c_Ignicoccus_hospitali', 0.01188037689471528, 0.0565342072920934, 0.06636624334289226, 0.050798852929127405, 0.1269971323228185, 0.17042195821384679, 0.2072920934043425, 0.3097091356001639], ['eS04c_Methanococcus_aeolic', 0.012409927942353884, 0.055644515612489995, 0.06645316253002402, 0.0500400320256205, 0.12730184147317855, 0.1709367493995196, 0.20736589271417133, 0.3098478783026421], ['eS04c_Methanocaldococcus_j', 0.012524461839530333, 0.0547945205479452, 0.06653620352250489, 0.05009784735812133, 0.12759295499021525, 0.17142857142857143, 0.2070450097847358, 0.30998043052837576], ['eS04c_Aenigmarchaeota_arch', 0.01225114854517611, 0.053981623277182235, 0.06623277182235834, 0.05053598774885146, 0.1274885145482389, 0.1726646248085758, 0.20673813169984687, 0.3101071975497703], ['eS04c_Odinarchaeota_archae', 0.012364181341326339, 0.0532034469838891, 0.06594230048707381, 0.0502060696890221, 0.12701386286998875, 0.1738478831022855, 0.20644436118396403, 0.3109778943424504], ['eS04c_Thermofilum_pendens_', 0.01210564930300807, 0.05392516507703595, 0.06603081438004402, 0.04952311078503301, 0.12655906089508437, 0.17461482024944974, 0.2057960381511372, 0.31144534115920763], ['eS04c_Pyrobaculum_calidifo', 0.011879049676025918, 0.05399568034557235, 0.0665946724262059, 0.04967602591792657, 0.12670986321094313, 0.1742260619150468, 0.20518358531317496, 0.3117350611951044], ['eS04c_Caldivirga_maquiling', 0.012014134275618375, 0.05335689045936396, 0.06713780918727916, 0.04876325088339223, 0.12756183745583038, 0.17420494699646644, 0.20459363957597174, 0.31236749116607776], ['eS04c_Sulfolobus_acidocald', 0.011789181692094313, 0.05339805825242718, 0.06692094313453537, 0.04854368932038835, 0.1276005547850208, 0.1754507628294036, 0.2035367545076283, 0.3127600554785021], ['eS04c_Sulfolobus_tokodaii', 0.011572498298162015, 0.05343771272974813, 0.06637168141592921, 0.04901293396868618, 0.12763784887678692, 0.17665078284547311, 0.20285908781484002, 0.3124574540503744], ['eS04c_Metallosphaera_sedul', 0.011363636363636364, 0.053475935828877004, 0.0661764705882353, 0.04879679144385027, 0.12733957219251338, 0.17780748663101603, 0.20220588235294118, 0.31283422459893045], ['eS04c_Woesearchaeota_archa', 0.01149047931713723, 0.05285620485883125, 0.0659881812212738, 0.048588312541037425, 0.12705187130663165, 0.17793827971109652, 0.20157583716349312, 0.314510833880499], ['eS04c_NANEQ', 0.011312217194570135, 0.05397543632837751, 0.06690368455074337, 0.048157724628312866, 0.12669683257918551, 0.17808661926308986, 0.20135746606334842, 0.31351001939237233], ['eS04c_Halobacterium_salina', 0.011135857461024499, 0.053770283168946865, 0.06776964683423481, 0.04740693604836144, 0.12758510976773782, 0.17944638880050906, 0.19980909958638243, 0.31307667833280306], ['eS04c_Natronomonas_pharaon', 0.01096147823363608, 0.05324146570623239, 0.0685875352333229, 0.046664578766050735, 0.12840588787973692, 0.1807077983088005, 0.19793297839022864, 0.3134982774819919], ['eS04c_Haloferax_volcanii_D', 0.010792476102374344, 0.05303731113166821, 0.06938020351526364, 0.04594511255010793, 0.12920135676842429, 0.18193031144002467, 0.19642306506321308, 0.31329016342892385], ['eS04c_Halorubrum_lacusprof', 0.010628606134224112, 0.0525356817491649, 0.06984512602490131, 0.045551169146674765, 0.13027634375948982, 0.18281202550865472, 0.194959003947768, 0.3133920437291224], ['eS04c_Haloquadratum_walsby', 0.01046337817638266, 0.05201793721973094, 0.06995515695067264, 0.04514200298953662, 0.13064275037369208, 0.18445440956651718, 0.1937219730941704, 0.3136023916292975], ['eS04c_Haloarcula_marismort', 0.010309278350515464, 0.05154639175257732, 0.07098674521354933, 0.04447717231222386, 0.13166421207658321, 0.1849779086892489, 0.1923416789396171, 0.3136966126656848], ['eS04c_Methanospirillum_hun', 0.010165553296543712, 0.051118210862619806, 0.07115887307580598, 0.04385710136508859, 0.13099041533546327, 0.1864652918966018, 0.19227417949462677, 0.31397037467325006], ['eS04c_Methanoculleus_maris', 0.01031223145230593, 0.05127470638785448, 0.071326267545116, 0.04325408192494987, 0.13004869664852478, 0.18733887138355773, 0.1922085362360355, 0.3142366084216557], ['eS04c_Methanoregula_boonei', 0.010454930771404351, 0.05114439107092399, 0.07120655552415937, 0.042949985871715175, 0.1294150890081944, 0.18818875388527834, 0.19186210794009606, 0.3147781859282283], ['eS04c_Methanocorpusculum_l', 0.010593810984109284, 0.05073877892389183, 0.0719264008921104, 0.042375243936437136, 0.12879843880680233, 0.18901589071647618, 0.19124616671313074, 0.3153052690270421], ['eS04c_Methanosarcina_aceti', 0.010729023383768913, 0.05061898211829436, 0.07235213204951857, 0.042090784044016505, 0.1281980742778542, 0.19119669876203577, 0.1903713892709766, 0.3144429160935351], ['eS04c_Methanosarcina_mazei', 0.010860711376595167, 0.05023079011675265, 0.07276676622318762, 0.04181373879989139, 0.1276133586749932, 0.19332066250339397, 0.18951941352158566, 0.3138745587836003], ['eS04c_Methanosarcina_barke', 0.01098901098901099, 0.04985258643795229, 0.07290270704904851, 0.04154382203162691, 0.12677566336102922, 0.19565800053604931, 0.18868935942106674, 0.313588850174216], ['eS04c_Methanococcoides_bur', 0.011114051336332363, 0.04948399047367028, 0.07356443503572374, 0.04128076210637735, 0.1264884890182588, 0.19661286054511776, 0.18814501190791214, 0.31331039957660756], ['eS04c_Methanosaeta_thermop', 0.010974653775803502, 0.04886333943036321, 0.0739482623464855, 0.04154690357982754, 0.12646981970211654, 0.19754376796446302, 0.18735301802978835, 0.31330023517115235], ['eS04c_Archaeoglobus_fulgid', 0.010838709677419355, 0.04903225806451613, 0.07406451612903225, 0.04129032258064516, 0.12645161290322582, 0.19793548387096774, 0.1870967741935484, 0.31329032258064515], ['eS04c_Pacearchaeota_archae', 0.010719754977029096, 0.048749361919346604, 0.07401735579377233, 0.0408371618172537, 0.1260847371107708, 0.1978050025523226, 0.18734047983665136, 0.3144461459928535], ['eS04c_Diapherotrites_archa', 0.010842158345940494, 0.04815935451336359, 0.07362581946545638, 0.040342914775592535, 0.12556732223903178, 0.19919314170448815, 0.1875945537065053, 0.3146747352496218], ['eS04c_Micrarchaeum_acidiph', 0.010712506228201295, 0.04808171400099651, 0.0744892874937718, 0.0401096163428002, 0.12580966616841055, 0.19855505729945191, 0.1878425510712506, 0.3143996013951171], ['eS04c_Nanopusillus_acidilo', 0.010583312823037165, 0.04799409303470342, 0.07482156042333252, 0.03962589219788334, 0.1260152596603495, 0.1983755845434408, 0.18828451882845187, 0.31429977848880136], ['eS04c_Heimdallarchaeota_ar', 0.010666666666666666, 0.04896969696969697, 0.07466666666666667, 0.039757575757575755, 0.12581818181818183, 0.19975757575757574, 0.1886060606060606, 0.31175757575757573], ['eS28c_Odinarchaeota_archa', 0.010523798134417603, 0.04831380052618991, 0.07534082755321693, 0.03970342023439369, 0.12556804592202822, 0.19899545563262377, 0.18847165749820616, 0.3130829944989237], ['eS28c_Ferroplasma_acidiph', 0.01038470616001888, 0.04767524191645032, 0.0759971678074109, 0.03965069624734482, 0.1260325702147746, 0.19801746518763275, 0.1878687750767052, 0.3143733773896625], ['eS28c_Thermoplasma_volcan', 0.010249242953645469, 0.047053342650826925, 0.07663638481248544, 0.03959934777544841, 0.12601910086186816, 0.19753086419753085, 0.18704868390402982, 0.3158630328441649], ['eS28c_Picrophilus_torridu', 0.010117268337548863, 0.04644745918601977, 0.07725914003219131, 0.03954932168314555, 0.1262359163025983, 0.19682685674867786, 0.18647965049436652, 0.31708438721545185], ['eS28c_Archaeoglobus_fulgi', 0.009988649262202044, 0.04608399545970488, 0.07786606129398412, 0.0395005675368899, 0.1264472190692395, 0.19591373439273552, 0.18615209988649262, 0.3180476730987514], ['eS28e_Yarrowia_lipolytica', 0.009863259358888142, 0.04572965702757229, 0.07868190988567586, 0.03945303743555257, 0.1268773817529702, 0.1947993723380408, 0.1862811028917283, 0.31831427930957185], ['eS28c_Pyrococcus_furiosus', 0.00974097852557007, 0.04538410449413328, 0.07947752933362852, 0.03940668585344255, 0.1272968784591543, 0.19371264113349568, 0.18640690723931813, 0.31857427496125745], ['eS28c_Thermococcus', 0.009621692543188278, 0.045047015088563305, 0.08047233763393834, 0.03936146949486114, 0.1274874261972447, 0.19243385086376558, 0.1865296304395364, 0.31904657773890227], ['eS28c_Thorarchaeota_archa', 0.009505292719809894, 0.04450205227910996, 0.08101101749837979, 0.03953337653920933, 0.1272413048174552, 0.1918340894361633, 0.18643335493627133, 0.3199395117736012], ['eS28c_Thermofilum_pendens', 0.009391675560298827, 0.0439701173959445, 0.08153681963713981, 0.03948772678762007, 0.1272145144076841, 0.1912486659551761, 0.18633938100320172, 0.3208110992529349], ['eS28c_Ignicoccus_hospital', 0.009280742459396751, 0.04366167475216199, 0.08226112634465303, 0.03944315545243619, 0.12718835688673275, 0.1904661463826197, 0.1862476270828939, 0.32145117063910567], ['eS28c_Sulfolobus_tokodaii', 0.009172399416301855, 0.04356889722743381, 0.08276005836981447, 0.039399624765478425, 0.12737127371273713, 0.1897018970189702, 0.1861580154263081, 0.321867834062956], ['eS28c_Sulfolobus_acidocal', 0.00906655676900886, 0.043478260869565216, 0.08324747578817226, 0.03935709870183392, 0.12754996909128374, 0.18895528539048012, 0.18607047187306822, 0.3222748815165877], ['eS28c_Metallosphaera_sedu', 0.00896312894683235, 0.043389692401711144, 0.08372377266245672, 0.03931554288042371, 0.12772458749236096, 0.18822570788347934, 0.18598492564677124, 0.3226726420859646], ['eS28c_Aeropyrum_pernix_K1', 0.008862034239677744, 0.0433031218529708, 0.08439073514602216, 0.03927492447129909, 0.12769385699899294, 0.18751258811681773, 0.18569989929506545, 0.32326283987915405], ['eS28c_Caldivirga_maquilin', 0.00876319458275244, 0.04282015534754033, 0.08484365664210317, 0.03903604859589723, 0.1276638119896435, 0.1866162119099781, 0.1856203943437562, 0.32463652658832903], ['eS28c_Nitrosopumilus_mari', 0.008666535355524916, 0.04234784321449675, 0.0852865865668702, 0.03860547567461099, 0.12763442978136694, 0.18613354343116015, 0.185739610005909, 0.32558597597006106], ['eS28c_Staphylothermus_mar', 0.008571985193843756, 0.04208065458796026, 0.08591466978375219, 0.0385739333722969, 0.12760568868108318, 0.1852717708942139, 0.18566140658484317, 0.32631989090200664], ['eS28c_Cenarchaeum_symbios', 0.008479475814222393, 0.04162651763345539, 0.0863364810175371, 0.03815764116400077, 0.12757756793216418, 0.1850067450375795, 0.18558489111582194, 0.32723068028521873]]
    return JsonResponse(test, safe = False)

def propensities(request):
    # how to not hardcode this?
    propensity_data = "http://127.0.0.1:8000/propensity-data"

    # where does the context variable come from? what does it do?
    context = {"propensity_data" : propensity_data}
    
    return render(request, 'alignments/propensities.html', context)