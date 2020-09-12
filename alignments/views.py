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

def propensity_data(request):
    test = [[1, 5.980671641791045], [2, 2.1377815468113983], [3, -0.9052035278154678], [4, 1.0471438263229307], [5, -0.047062415196743354], [6, 0.19905698778833122], [7, 0.1305359565807329], [8, -0.5349999999999999], [9, 0.2888805970149256], [10, 0.43440298507462705], [11, -0.3238059701492536], [12, -0.4573880597014924], [13, 0.13216417910447778], [14, 0.015000000000000124], [15, -0.0036567164179102524], [16, 0.18141791044776143], [17, 0.15305970149253753], [18, 0.03888059701492559], [19, 0.06574626865671665], [20, -0.05514925373134308], [21, 0.002313432835821085], [22, 0.1597761194029853], [23, 0.09335820895522409], [24, -0.21186567164179076], [25, 0.24932835820895538], [26, 0.1978358208955226], [27, -0.14171641791044753], [28, 0.19932835820895542], [29, -0.492462686567164], [30, 0.15380597014925396], [31, 0.22843283582089574], [32, 0.23291044776119427], [33, -0.17902985074626843], [34, 0.07012890094979672], [35, 0.1331411126187248], [36, -0.12216417910447745], [37, 0.042310719131614896], [38, 0.2275508819538672], [39, 1.2058412483039354], [40, 0.6148643147896881], [41, -0.7953120759837174], [42, 0.8824355495251022], [43, 1.9668181818181818], [44, 0.8431275440976935], [45, 1.5326526458616012], [46, 0.9046200814111265], [47, 4.566818181818181], [48, 1.7483107191316152], [49, -1.552910447761194], [50, 2.221309362279512], [51, 0.23976255088195397], [52, 1.2432903663500678], [53, 2.0717842605156034], [54, 1.519274084124831], [55, 1.157537313432836], [56, -0.47584124830393465], [57, 1.8650000000000002], [58, 3.4351492537313435], [59, 2.153873812754411], [60, -0.4151221166892807], [61, 0.7143080054274086], [62, 1.717605156037992], [63, 1.8277272727272729], [64, 0.7339009497964725], [65, 0.13977611940298526], [66, 1.5322727272727275], [67, 0.2240909090909094], [68, 1.6196811397557669], [69, 0.3731818181818184], [70, 0.5489077340569879], [71, 1.2130746268656716], [72, 0.784267299864315], [73, -0.19307327001356825], [74, 1.7800067842605163], [75, -0.03816146540027117], [76, 2.1493962008141114], [77, 1.3844029850746271], [78, 4.722259158751696], [79, 0.227279511533243], [80, 0.07829715061058357], [81, 3.949124830393487], [82, 0.5637788331071913], [83, 2.5922727272727273], [84, 5.065], [85, 6.428636363636364], [86, 0.2234803256445048], [87, 4.810454545454546], [88, 2.0101831750339216], [89, 3.7271438263229313], [90, -1.426451831750339], [91, 3.823208955223881], [92, 3.894036635006785], [93, 0.5784328358208958], [94, 3.507062415196744], [95, 0.15699457259158756], [96, 2.98331750339213], [97, 0.6850814111261874], [98, 3.8092876526458617], [99, 5.421309362279512], [100, 2.4281207598371783], [101, 1.9607937584803257], [102, 0.8297218453188604], [103, 3.7922727272727275], [104, 3.1719199457259157], [105, 2.1857598371777485], [106, 3.4655427408412485], [107, 2.844647218453188], [108, 1.5561804613297154], [109, 4.348582089552239], [110, 2.764864314789688], [111, 0.2886092265943014], [112, 2.465542740841249], [113, 6.319545454545455], [114, 2.244104477611941], [115, 1.0147964721845322], [116, 3.9010922659430123], [117, 0.38955902306648593], [118, -0.1822184531886022], [119, 5.766221166892809], [120, 0.5939009497964723], [121, 1.456424694708277], [122, 2.8457327001356862], [123, 1.2995725915875171], [124, 2.8769402985074626], [125, 2.4462754409769336], [126, 1.0131682496607872], [127, 6.9754477611940295], [128, 5.85523066485753], [129, 0.8617435549525105], [130, 3.5024491180461332], [131, 4.877211668928087], [132, 4.099735413839892], [133, 0.8620149253731345], [134, 1.6332496607869749], [135, 4.203670284938942], [136, 7.6104545454545445], [137, -0.7868317503392128], [138, 6.8104545454545455], [139, 2.752109905020353], [140, -0.04626187245590212], [141, 4.508962008141113], [142, 0.2191383989145185], [143, 4.234063772048847], [144, 4.502449118046133], [145, 6.065], [146, -0.023466757123473327], [147, 1.3187313432835823], [148, 4.794172320217096], [149, 5.299192672998644], [150, -0.030251017639077175], [151, 7.065], [152, 4.465271370420624], [153, 2.3816892808683856], [154, 1.9483107191316145], [155, 2.5390841248303935], [156, 5.792272727272727], [157, 3.4017706919945727], [158, 2.374090909090909], [159, 7.065], [160, 1.8082835820895524], [161, 1.7892876526458614], [162, 3.17544776119403], [163, 0.13067164179104498], [164, 3.3898303934871103], [165, 2.5423405698778834], [166, 1.610997286295794], [167, 1.9906445047489827], [168, 2.2441044776119403], [169, 1.2264654002713704], [170, 1.1246472184531888], [171, -0.42428086838534573], [172, 2.1458683853459974], [173, 0.828636363636364], [174, 3.4443758480325646], [175, 2.7124898236092263], [176, 4.293222523744912], [177, 2.8685278154681133], [178, 5.065], [179, 5.2579443690637735], [180, 2.6210379918588873], [181, 1.6150678426051561], [182, 2.5377272727272735], [183, 0.17789009497964736], [184, 2.3844029850746273], [185, 0.3496675712347355], [186, 0.051160108548168506], [187, 1.1176458616010858], [188, 7.106519674355495], [189, 2.8099118046132974], [190, 4.81588195386703], [191, 4.480739484396201], [192, 4.347496607869743], [193, 7.065], [194, -0.24843283582089537], [195, 4.267713704206242], [196, 1.503263229308006], [197, 1.0650000000000002], [198, 0.19580054274084147], [199, -0.5170895522388058], [200, 0.7534667571234738], [201, 0.30896200814111285], [202, -0.26335820895522355], [203, 0.5651356852103124], [204, -0.42753731343283563], [205, 1.5122184531886023], [206, -2.2685142469470825], [207, 5.883181818181818], [208, 5.010454545454546], [209, 7.065], [210, 1.2598439620081412], [211, 1.4432903663500678], [212, -0.8343215739484394], [213, 1.2807394843962008], [214, 0.8191383989145185], [215, -0.444090909090909], [216, 0.6536024423337858], [217, -1.4788263229308005], [218, 6.846818181818181], [219, 0.9203595658073271], [220, 3.5029918588873823], [221, 4.665000000000001], [222, 6.265], [223, 2.3702917232021714], [224, 2.705976933514248], [225, 5.50489145183175], [226, 6.275583446404341], [227, 1.104348710990502], [228, 3.0527883310719135], [229, 4.270698778833108], [230, 9.065], [231, 7.919545454545454], [232, 2.9488534599728635], [233, 7.065], [234, 6.7357734056987795], [235, 5.627252374491181], [236, 0.6555834464043421], [237, 1.847388059701493], [238, 0.3599796472184534], [239, 0.8593012211668929], [240, 1.5787042062415204], [241, -0.37367028493894155], [242, -0.4419063772048846], [243, 0.18440298507462705], [244, -0.09863636363636354], [245, 0.8424762550881958], [246, 2.3165603799185894], [247, 3.7678493894165532], [248, -2.2172795115332424], [249, 4.119545454545454], [250, 0.5013636363636365], [251, 2.2924084124830393], [252, 0.38174355495251033], [253, -0.019341926729986015], [254, -0.7873744911804611], [255, 0.06500000000000017], [256, -0.19741519674355476], [257, 0.5195454545454548], [258, 2.283181818181818], [259, 1.2663839891451834], [260, 0.5825033921302578], [261, 1.5184599728629582], [262, 0.9741180461329717], [263, 0.4235888738127547], [264, -0.23730664857530512], [265, 1.2642808683853461], [266, -0.4820284938941653], [267, -0.662272727272727], [268, -0.17415875169606482], [269, 0.2306987788331073], [270, 0.36248982360922677], [271, 0.44641112618724577], [272, 0.3522727272727275], [273, 0.3126879240162825]] 
    return JsonResponse(test, safe = False)

def propensities(request):
    # how to not hardcode this?
    propensity_data = "http://127.0.0.1:8000/propensity-data"

    # where does the context variable come from? what does it do?
    context = {"propensity_data" : propensity_data}
    
    return render(request, 'alignments/propensities.html', context)