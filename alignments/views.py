import re, os
import datetime
from subprocess import Popen, PIPE
from Bio import AlignIO
from io import StringIO
from Bio.SeqUtils import IUPACData

from django.shortcuts import render
from django.http import HttpResponse, Http404, JsonResponse, HttpResponseServerError
from django.urls import reverse

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

    pipe = Popen("mafft --quiet --addfull " + ebiFileName + " --mapout " + alignmentFileName + "; cat " + mappingFileName, stdout=PIPE, shell=True)
    output = pipe.communicate()[0]
    text = output.decode("ascii").split('\n#')[1]

    mapping, firstLine, badMapping = dict(), True, 0
    for line in text.split('\n'):
        if firstLine:
            firstLine = False
            continue
        row = line.split(', ')
        if len(row) < 3:
            continue
        if row[2] == '-':
            badMapping += 1
            continue
        if row[1] == '-':
            raise Http404("Mapping did not work properly.")
        mapping[int(row[2])] = int(row[1]) + shiftIndexBy

    if badMapping > 0:
        mapping['BadMappingPositions'] = badMapping

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
    fasta = request.POST["fasta"]
    concat_fasta = re.sub(r'\\n','\n', fasta,flags=re.M)
    list_for_topology_viewer = calculate_twincons(concat_fasta)
    return JsonResponse(list_for_topology_viewer, safe = False)

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

def index_test(request):
    some_Alignments = Alignment.objects.all()
    superKingdoms = Taxgroups.objects.raw('SELECT * FROM SEREB.TaxGroups WHERE\
         SEREB.TaxGroups.groupLevel = "superkingdom";')
    
    context = {
        'props': list(Taxgroups.objects.values('taxgroup_id', 'groupname')),
        'some_Alignments': some_Alignments,
        'superKingdoms': superKingdoms
    }
    return render(request, 'alignments/index_test.html', context)

def index_orthologs(request):
    print ("request.method == 'POST': " + str(request.method == 'POST'))
    print ("'custom_propensity_data' in request.FILES: " + str('custom_propensity_data' in request.FILES))
    print ("request.method: " + str(request.method))
    for x in request.FILES:
        print ("x: " + str(x))
    if request.method == 'GET' and 'custom_propensity_data' in request.FILES:
        propensity_indices_file = request.FILES['custom_propensity_data']
        propensity_indices_string = ''
        for propensity_part in propensity_indices_file.chunks():
            propensity_indices_string += propensity_part.decode()
        print ("propensity_indices_string: " + propensity_indices_string)
    # if request.method == 'POST' and 'custom_propensity_data' in request.FILES:
    #     propensity_indices_file = request.FILES['custom_propensity_data']
    #     propensity_indices_string = ''
    #     for propensity_part in propensity_indices_file.chunks():
    #         propensity_indices_string += propensity_part.decode()
    #     print ("propensity_indices_string: " + propensity_indices_string)
    return render(request, 'alignments/index_orthologs.html')
def index(request):
    return render(request, 'alignments/index.html')

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

def extract_species_list(fastastring):
    '''Filters out species from a fastastring'''
    unf_species_list = [x.split('\\')[0] for x in fastastring.split('>')[1:]]
    filtered_spec_list = [re.sub('_',' ', re.sub(r'^.*?_', '', x)) for x in unf_species_list]
    return filtered_spec_list

def extract_gap_only_cols(fastastring):
    '''Extracts positions in the fastastring that are only gaps'''
    unf_seq_list = [x.split('\\n')[1] for x in fastastring.split('>')[1:]]
    list_for_intersect = list()
    for sequence in unf_seq_list:
        iterator = re.finditer('-', sequence)
        gap_positions = [m.start(0) for m in iterator]
        list_for_intersect.append(gap_positions)
    gap_only_cols = list(set(list_for_intersect[0]).intersection(*list_for_intersect))
    return gap_only_cols

def construct_dict_for_json_response(response_data):
    '''Takes list of datas for Json response.
    Check their types and assigns names to each.
    Returns them as a dictionary.'''
    response_dict = dict()
    for entry in response_data:
        if type(entry) == str:
            response_dict['Alignment'] = entry
            continue
        if type(entry) == bool:
            response_dict['TwinCons'] = entry
            continue
        if type(entry) == list:
            if all(isinstance(item, int) for item in entry):
                response_dict['Gap-only columns'] = entry
                continue
            if all(isinstance(item, list) for item in entry):
                response_dict['AA frequencies'] = entry
                continue
            if all(isinstance(item, str) for item in entry):
                response_dict['Sequence names'] = entry
                continue
    return response_dict

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
    if internal:
        return fastastring
    
    gap_only_cols = extract_gap_only_cols(fastastring)
    filtered_spec_list = extract_species_list(fastastring)

    concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)
    response_dict = construct_dict_for_json_response([concat_fasta,filtered_spec_list,gap_only_cols,frequency_list])

    return JsonResponse(response_dict, safe = False)

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

def validate_fasta_string(fastaString):
    malicious_strings = [
        "eval\(unescape",
        "base64_decode\(",
        "substr\(md5\(strrev\(",
        "cwd = @getcwd\(\);",
        "chr\(\(ord\(",
        "gzinflate\(base64_decode\(",
        "php_uname\(\)\" \"] = chr\(ord\(",
        "cwd\[strlen\(\$cwd\)",
        "ini_get\('safe_mode'\);",
        "=\"\x62\"",
        "\"+ r + \"&r=\" + document.referrer;\"",
        "if\(strtoupper\(substr\(PHP_OS, 0, 3\) \) == \"WIN\"\)",
        "window.top.location.href=\"http://",
        "@ini_get\(\"disable_functions\"\)",
        "\$g3='';\$g3.=\$r;\$g3.=\$h;\$g3.=\$y",
    ]
    for mstring in malicious_strings:
        regex = re.compile(mstring)
        if re.search(regex, fastaString):
            return True
    return False

def handle_custom_upload_alignment(request):
	if request.method == 'POST' and 'custom_aln_file' in request.FILES:
		aln_file = request.FILES['custom_aln_file']
		alignment_string = ''
		for aln_part in aln_file.chunks():
			alignment_string += aln_part.decode()
		alignments = list(AlignIO.parse(StringIO(alignment_string), 'fasta'))
		if len(alignments) == 0:
			return HttpResponseServerError("Wasn't able to parse the alignment file! Is your file in fasta format?")
		if len(alignments) > 1:
			return HttpResponseServerError("Alignment file had more than one alignments!\nPlease upload a single alignment.")
		fastastring = alignments[0].format("fasta")
		request.session['custom_alignment_file'] = fastastring
		return HttpResponse('Success!')
	if request.method == 'GET':
		from alignments.Shannon import gap_adjusted_frequency
		from TwinCons.bin.TwinCons import slice_by_name
		fastastring = request.session.get('custom_alignment_file')
		alignment_obj = AlignIO.read(StringIO(fastastring), 'fasta')
		sliced_alns = slice_by_name(alignment_obj)
		twc = False
		if len(sliced_alns.keys()) == 2:
			twc = True			
		fastastring = fastastring.replace('\n','\\n')
		gap_only_cols = extract_gap_only_cols(fastastring)
		filtered_spec_list = extract_species_list(fastastring)
		concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)
		frequency_list = list()
		for i in range(0, alignment_obj.get_alignment_length()):
			frequency_list.append(gap_adjusted_frequency(alignment_obj[:,i], IUPACData.protein_letters))
		
		response_dict = construct_dict_for_json_response([concat_fasta,filtered_spec_list,gap_only_cols,frequency_list,twc])
		return JsonResponse(response_dict, safe = False)

# trims fasta by a list of indices
def trim_fasta_by_index(input_file, indices):
    from Bio import AlignIO
    align = AlignIO.read(input_file, "fasta")
    trimmed_align = align[:,indices[0]:indices[0]+1] # initialize align object
    for i in indices[1:]:
        trimmed_align += align[:,i:i+1]
    return trimmed_align

# TODO: change this back
def propensity_data(request, aln_id, tax_group):
    from io import StringIO
    import alignments.propensities as propensities

    fastastring = simple_fasta(request, aln_id, tax_group, internal=True).replace('\\n', '\n')
    fasta = StringIO(fastastring)

    if request.method == 'POST' and 'indices' in request.POST:
        indices = request.POST['indices']
        trimmed_fasta = trim_fasta_by_index(fasta, indices)
        fasta = StringIO(trimmed_fasta.format('fasta'))

    aa = propensities.aa_composition(fasta, reduced = False)
    fasta.seek(0) # reload the fasta object
    red_aa = propensities.aa_composition(fasta, reduced = False)

    data = {
        'aln_id' : aln_id,
        'tax_group' : tax_group,
        'reduced alphabet' : red_aa,
        'amino acid' : aa
    }
    return JsonResponse(data)

def propensities(request, align_name, tax_group):
    aln_id = Alignment.objects.filter(name = align_name)[0].aln_id
    propensity_data = reverse('alignments:propensity_data', kwargs={'aln_id': aln_id, 'tax_group' : tax_group})

    names = []
    if type(tax_group) == int:
        tax_group = str(tax_group)
    for group in tax_group.split(','):
        names.append(Taxgroups.objects.get(pk=group).groupname)

    context = {
        "propensity_data" : propensity_data, 
        "align_name" : align_name,
        "tax_name" : ', '.join(names)
    }
    
    return render(request, 'alignments/propensities.html', context)