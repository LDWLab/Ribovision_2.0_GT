import contextlib
import re, os, warnings, io, base64, json
import datetime
import urllib.request
import shutil
from subprocess import Popen, PIPE
from pdbecif.mmcif_io import CifFileReader
from Bio import AlignIO, BiopythonDeprecationWarning, PDB
from io import StringIO

from django.shortcuts import render
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from django.urls import reverse
from django.contrib.sites.shortcuts import get_current_site

from alignments.models import *
from alignments.taxonomy_views import *
from alignments.residue_api import *
from alignments.structure_api import *
from alignments.fold_api import *
from alignments.runal2co import executeAl2co
import alignments.alignment_query_and_build as aqab
from twincons.TwinCons import slice_by_name
from django.db import connection
import time
from xml.dom import minidom

class c:
    structureObj = None

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
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonDeprecationWarning)
    from twincons import TwinCons
    list_for_phymeas = ['-as',alignment, '-nc', '-r', '-mx', 'blastn']
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
    ### BE CAREFUL WHEN MERGING THE FOLLOWING LINES TO PUBLIC; PATHS ARE HARDCODED FOR THE APACHE SERVER ###
    alignmentFileName = "/home/sumon/repos/Ribovision_2.0_GT/static/alignment" + fileNameSuffix + ".txt"
    ebiFileName = "/home/sumon/repos/Ribovision_2.0_GT/static/ebi_sequence" + fileNameSuffix + ".txt"
    #alignmentFileName = "./static/alignment" + fileNameSuffix + ".txt"
    #ebiFileName = "./static/ebi_sequence" + fileNameSuffix + ".txt"
    mappingFileName = ebiFileName + ".map"
    fasta = re.sub('>Structure sequence[\s\S]*?>','>',fasta)
    fh = open(alignmentFileName, "w")
    fh.write(fasta)
    fh.close()

    fh = open(ebiFileName, "w")
    fh.write(">Structure sequence\n")
    fh.write(ebi_sequence)
    fh.close()

    shiftIndexBy = 0
    if startIndex > 1:
        shiftIndexBy = startIndex - 1

    print("Mafft")
    pipe = Popen("mafft --preservecase --quiet --addfull " + ebiFileName + " --mapout " + alignmentFileName + "; cat " + mappingFileName, stdout=PIPE, shell=True)
    output = pipe.communicate()[0]
    decoded_text = output.decode("ascii")
    
    if len(decoded_text) <= 0:
        for removeFile in [alignmentFileName, ebiFileName]:
            os.remove(removeFile)
        return HttpResponseServerError("Failed mapping the polymer sequence to the alignment!\nTry a different structure.")

    text = decoded_text.split('\n#')[1]
    amendedAln = re.sub('>Structure sequence$','',decoded_text.split('\n#')[0])
    outputDict, mapping, firstLine, badMapping = dict(), dict(), True, 0
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
            return HttpResponseServerError("Failed mapping the polymer sequence to the alignment!\nTry a different structure.")
        mapping[int(row[2])] = int(row[1]) + shiftIndexBy

    outputDict["structureMapping"] = mapping
    if badMapping > 0:
        outputDict['BadMappingPositions'] = badMapping

    for removeFile in [alignmentFileName, ebiFileName, mappingFileName]:
        os.remove(removeFile)
    outputDict["amendedAln"] = f'>Structure sequence{amendedAln.split(">Structure sequence")[1]}{amendedAln.split(">Structure sequence")[0]}'
    return outputDict

def request_post_data(post_data):
    fasta = post_data["fasta"]
    ebi_sequence = post_data["ebi_sequence"]
    startIndex = int(post_data["startIndex"])
    return fasta, ebi_sequence, startIndex

def make_map_from_alnix_to_sequenceix(request):
    fasta, ebi_sequence, startIndex = request_post_data(request.POST)
    mapping = constructEbiAlignmentString(fasta, ebi_sequence, startIndex)
    if type(mapping) != dict:
        return mapping
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
    superKingdoms = Taxgroups.objects.raw('SELECT * FROM TaxGroups WHERE\
         TaxGroups.groupLevel = "superkingdom";')
    
    context = {
        'props': list(Taxgroups.objects.values('taxgroup_id', 'groupname')),
        'some_Alignments': some_Alignments,
        'superKingdoms': superKingdoms
    }
    return render(request, 'alignments/index_test.html', context)

def proteinTypes(request):
    if request.method == 'POST' and 'taxIDs' in request.POST:
        print("request.POST['taxIDs']", request.POST['taxIDs'])
        taxIDs = request.POST['taxIDs']
        return proteinTypesDirect(request, taxIDs)
    else:
        return []

def allProteinTypes(request):
    allProteinTypes = []
    with connection.cursor() as cursor:
        sql = "select distinct(MoleculeGroup) from Nomenclature order by MoleculeGroup ASC;"
        cursor.execute(sql)
        for row in cursor.fetchall():
            allProteinTypes.append(row[0])
    context = {
        "allProteinTypes" : allProteinTypes
    }
    return JsonResponse(context)

def proteinTypesDirect(request, concatenatedTaxIds):
    results = []
    with connection.cursor() as cursor:
        if concatenatedTaxIds.endswith(','):
            concatenatedTaxIds = concatenatedTaxIds[:-1]
        for taxID in concatenatedTaxIds.split(','):
            if not taxID:
                continue
            taxID = int(taxID)
            proteinTypesList = []
            sql = "use DESIRE;"
            cursor.execute(sql)
            # sql = 'SET sql_mode=(SELECT REPLACE(@@sql_mode,\'ONLY_FULL_GROUP_BY\',\'\'));'
            sql = 'select Nomenclature.MoleculeGroup from TaxGroups join Species_TaxGroup on Species_TaxGroup.taxgroup_id = TaxGroups.taxgroup_id join Species on Species.strain_id = Species_TaxGroup.strain_id join Species_Polymer on Species.strain_id = Species_Polymer.strain_id join Polymer_Data on Polymer_Data.strain_id = Species_Polymer.strain_id and Polymer_Data.GI = Species_Polymer.GI and Polymer_Data.nomgd_id = Species_Polymer.nomgd_id join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id where TaxGroups.taxgroup_id = ' + str(taxID) + ' group by Nomenclature.MoleculeGroup;'

            cursor.execute(sql)
            # results = Taxgroups.objects.raw(sql)
            for row in cursor.fetchall():
                proteinTypesList.append(row[0])
            filteredList = [s for s in proteinTypesList if s.find("RNA") !=-1]
            results.append(filteredList)
    context = {'results' : results}
    return JsonResponse(context)

def getAlignmentsFilterByProteinTypeAndTaxIds(request):

    if request.method == 'POST' and 'selectedProteinType' in request.POST and 'taxIDs' in request.POST:
        return getAlignmentsFilterByProteinTypeAndTaxIdsDirect(request, request.POST['selectedProteinType'], request.POST['taxIDs'])
    else:
        return []

def getAlignmentsFilterByProteinTypeAndTaxIdsDirect(request, concatenatedProteinTypes, concatenatedTaxIds):
    results = []
    with connection.cursor() as cursor:
        if concatenatedProteinTypes.endswith(','):
            concatenatedProteinTypes = concatenatedProteinTypes[:-1]
        if concatenatedTaxIds.endswith(','):
            concatenatedTaxIds = concatenatedTaxIds[:-1]
        proteinTypes = concatenatedProteinTypes.split(',')
        concatenatedProteinTypes = '\'' + proteinTypes[0] + '\''
        for i in range(1, len(proteinTypes)):
            concatenatedProteinTypes += ', \'' + proteinTypes[i] + '\''
        for taxID in concatenatedTaxIds.split(','):
            alignmentNamesAndPrimaryKeys = []
            sql = "select Alignment.Name, Alignment.Aln_id from Nomenclature join Polymer_Data on Polymer_Data.nomgd_id = Nomenclature.nom_id join Polymer_Alignments on Polymer_Alignments.PData_id = Polymer_Data.PData_id join Alignment on Alignment.Aln_id = Polymer_Alignments.Aln_id join Species_Polymer on Species_Polymer.strain_id = Polymer_Data.strain_id and Species_Polymer.GI = Polymer_Data.GI and Species_Polymer.nomgd_id = Polymer_Data.nomgd_id join Species on Species.strain_id = Species_Polymer.strain_id join Species_TaxGroup on Species.strain_id = Species_TaxGroup.strain_id join TaxGroups on TaxGroups.taxgroup_id = Species_TaxGroup.taxgroup_id where Nomenclature.MoleculeGroup in (" + concatenatedProteinTypes + ") and Alignment.Method = 'GSD_LSD_rRNA' and TaxGroups.taxgroup_id = " + taxID + " group by Alignment.Aln_id;"

            cursor.execute(sql)
            for row in cursor.fetchall():
                alignmentNamesAndPrimaryKeys.append([row[0], row[1]])
            results.append(alignmentNamesAndPrimaryKeys)

    context = {'results' : results}
    return JsonResponse(context)


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
    print("fastastring views", fastastring)
    unf_seq_list = [x.split('\\n')[1] for x in fastastring.split('>')[1:]]
    list_for_intersect = list()
    for sequence in unf_seq_list:
        iterator = re.finditer('-', sequence)
        gap_positions = [m.start(0) for m in iterator]
        list_for_intersect.append(gap_positions)
    print('list_for_intersect', list_for_intersect)
    gap_only_cols = sorted(list(set(list_for_intersect[0]).intersection(*list_for_intersect)))
    print("views gap_only_cols", gap_only_cols)
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
    print("Simple Fasta tax_group", tax_group)
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
    # print('fastastring', fastastring)
    
    if internal:
        return fastastring
    
    # concat_fasta, twc, gap_only_cols, filtered_spec_list, alignment_obj = calculateFastaProps(fastastring)
    # print('gap_only_cols', gap_only_cols)
    # # print('alignment_obj', alignment_obj)
    # response_dict = construct_dict_for_json_response([concat_fasta,filtered_spec_list,gap_only_cols,frequency_list,twc])
    
    concat_fasta, twc, gap_only_cols, filtered_spec_list, alignment_obj, fr_list = calculateFastaProps(fastastring, frequency_list)
    print('GOC', gap_only_cols)
    response_dict = construct_dict_for_json_response([concat_fasta,filtered_spec_list,gap_only_cols,fr_list,twc])

    return JsonResponse(response_dict, safe = False)

def calculateFastaProps(fastastring, frequency_list):
    concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)
    alignment_obj = AlignIO.read(StringIO(concat_fasta), 'fasta')
    twc = False
    if (len(alignment_obj) < 1000):
        sliced_alns = slice_by_name(alignment_obj)
        if len(sliced_alns.keys()) == 2:
            twc = True
    gap_only_cols = extract_gap_only_cols(fastastring)
    
    removed_gaps = gap_only_cols[:] 
    print(gap_only_cols)
    mapped_dict = {}
    pos = 0
    
    for i in range(len(alignment_obj[0])):
        if i in gap_only_cols:
            continue
        mapped_dict[i] = pos
        pos += 1
    # print(mapped_dict)
    
    for gap in removed_gaps[::-1]:
        alignment_obj = alignment_obj[:, :gap] + alignment_obj[:, gap+1:]
        gap_only_cols.remove(gap)
        frequency_list.pop(gap)
        
    # for i in range(max(gap_only_cols)+1)[::-1]:
    #     if i in gap_only_cols:
    #         alignment_obj = alignment_obj[:, :i] + alignment_obj[:, i+1:]
    #         gap_only_cols.remove(i)
    #         frequency_list.pop(i)
    #         deleted_gaps += 1
        # mapped_dict[i] = (i-deleted_gaps)
    
    
    file_path = "../data_mapping.json"
    # Removing the file if it exists
    if os.path.exists(file_path):
        os.remove(file_path)
    
    # Convert dictionary to JSON string
    json_string = json.dumps(mapped_dict)
    
    with open(file_path, "w") as json_file:
        json_file.write(json_string)
        
    # for gap in removed_gaps[::-1]:
    #     # print(f'Removing: {gap}')
    #     alignment_obj = alignment_obj[:, :gap] + alignment_obj[:, gap+1:]
    #     gap_only_cols.remove(gap)
    #     frequency_list.pop(gap)
    # print("Post removal")
    # print(gap_only_cols)
    
    # print('gap_only_cols', gap_only_cols)
    # print(10 * '-')
    string_writer = StringIO()
    # AlignmentWriter(string_writer).write_alignment(alignment_obj)

    AlignIO.write(alignment_obj, string_writer, "fasta")
    alignment_string = string_writer.getvalue()
    
    # print('alignment_obj', alignment_string)
    
    filtered_spec_list = extract_species_list(fastastring)
    # return concat_fasta, twc, gap_only_cols, filtered_spec_list, alignment_obj
    return alignment_string, twc, gap_only_cols, filtered_spec_list, alignment_obj, frequency_list

def rProtein(request, align_name, tax_group):
    #if tax_group == 0 - no filter
    align_id = Alignment.objects.filter(name = align_name)[0].aln_id
    fastastring = simple_fasta(request, align_id, tax_group, internal=True)
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
            return False
    return True
def getGenusFromStrainIdsDirect(request, concatenatedStrainIds):
    strainIds = concatenatedStrainIds.split(',')
    concatenatedStrainIds = '\'' + strainIds[0] + '\''
    for i in range(1, len(strainIds)):
        concatenatedStrainIds += ', \'' + strainIds[i] + '\''
    sql = "use DESIRE;"
    cursor.execute(sql)
    sql = "select TaxGroups.taxgroup_id from Species join Species_TaxGroup on Species.strain_id = Species_TaxGroup.strain_id join TaxGroups on TaxGroups.taxgroup_id = Species_TaxGroup.taxgroup_id where Species.strain_id in (" + concatenatedStrainIds + ") and TaxGroups.groupLevel = 'genus';"
    genusList = []
    with connection.cursor() as cursor:
        cursor.execute(sql)
        for row in cursor.fetchall():
            genusList.append(row[0])
    context = {
        'genusList' : genusList
    }
    return JsonResponse(context)

def string_fasta_two_strains(request, protein_type, aln_id, strain_id_0, strain_id_1):
    if type(strain_id_0) == int:
        strain_id_0 = str(strain_id_0)
    if type(strain_id_1) == int:
        strain_id_1 = str(strain_id_1)
    protein_type = protein_type.split(',')
    for i in range(len(protein_type)):
        protein_type[i] = '\'' + protein_type[i] + '\''
    protein_type=','.join(protein_type)
    sql = ""
    with connection.cursor() as cursor:
        cursor.execute(sql)

# trims fasta by a list of indices
def trim_fasta_by_index(input_file, indices):
    from Bio import AlignIO
    align = AlignIO.read(input_file, "fasta")
    trimmed_align = align[:,int(indices[0]):int(indices[0])] # initialize align object
    for i in indices.split(','):
        if (int(i)-1 < 0):
            continue
        trimmed_align += align[:,int(i)-1:int(i)]
    return trimmed_align

def flushSession (request):
    try:
        request.session.flush()
    except:
        return HttpResponseServerError ("Failed to flush the session!")
    return HttpResponse ("Success!")

def parse_string_structure(request, stringData, strucID):
    #if(c.structureObj):
    #    return c.structureObj
    from Bio.PDB import MMCIFParser
    parser = MMCIFParser()
    strucFile = io.StringIO(stringData)
 
    return parser.get_structure(strucID,strucFile)

def getAlignmentsFilterByProteinTypeDirect(request, concatenatedProteinTypes):
    proteinTypes = concatenatedProteinTypes.split(',')
    concatenatedProteinTypes = '\'' + proteinTypes[0] + '\''
    
    for i in range(1, len(proteinTypes)):
        concatenatedProteinTypes += ', \'' + proteinTypes[i] + '\''
    sql = "select distinct(Name) from Alignment join Polymer_Alignments on Polymer_Alignments.Aln_id = Alignment.Aln_id join Polymer_Data on Polymer_Data.PData_id = Polymer_Alignments.PData_id join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id where Nomenclature.MoleculeGroup in (" + concatenatedProteinTypes + ");"
    with connection.cursor() as cursor:
        cursor.execute(sql)
        results = []
        for row in cursor.fetchall():
            results.append(row[0])
    context = {
        'results' : results
    }
    return JsonResponse(context)

def getPairwiseAlignmentDirect(request, moleculeType, alignmentName, strainId0, strainId1):
    strainIds = [strainId0, strainId1]
    concatenatedStrainIds = '\'' + strainIds[0] + '\''
    for i in range(1, len(strainIds)):
        concatenatedStrainIds += ', \'' + strainIds[i] + '\''
    sql = "select TaxGroups.taxgroup_id, Species.strain from Species join Species_TaxGroup on Species.strain_id = Species_TaxGroup.strain_id join TaxGroups on TaxGroups.taxgroup_id = Species_TaxGroup.taxgroup_id where Species.strain_id in (" + concatenatedStrainIds + ") and TaxGroups.groupLevel = 'genus';"
    genusList = []
    speciesNames = []
    with connection.cursor() as cursor:
        cursor.execute(sql)
        for row in cursor.fetchall():
            genusList.append(row[0])
            speciesNames.append(row[1])
    if len(speciesNames) == 1:
        speciesName0 = speciesName1 = speciesNames[0]
    else:
        speciesName0, speciesName1 = speciesNames

    alignment = string_fasta(request, moleculeType, alignmentName, strainId0 + ',' + strainId1, internal=True)
    alignment = alignment.replace('\\n', '\n')
    PairwiseAlignment = ''
    if (alignment.endswith('\n')) :
        alignment = alignment[0:-1]
        alignmentLines = alignment.splitlines()
        modifiedStrainName0 = speciesName0.replace(' ', '_')
        modifiedStrainName1 = speciesName1.replace(' ', '_')
        for i in range(1, len(alignmentLines), 2):
            titleLine = alignmentLines[i - 1]
            alignmentLine = alignmentLines[i]
            if modifiedStrainName0 in titleLine or modifiedStrainName1 in titleLine:
                PairwiseAlignment += alignmentLine + '\n'
    context = {
        'PairwiseAlignment' : PairwiseAlignment
    }
    return JsonResponse(context)

def getStrainsFilterByMoleculeGroupAndAlignmentDirect(request, concatenatedMoleculeGroups, concatenatedAlignmentNames):
    moleculeGroups = concatenatedMoleculeGroups.split(',')
    concatenatedMoleculeGroups = '\'' + moleculeGroups[0] + '\''
    for i in range(1, len(moleculeGroups)):
        concatenatedMoleculeGroups += ', \'' + moleculeGroups[i] + '\''
    alignmentNames = concatenatedAlignmentNames.split(',')
    concatenatedAlignmentNames = '\'' + alignmentNames[0] + '\''
    for i in range(1, len(alignmentNames)):
        concatenatedAlignmentNames += ', \'' + alignmentNames[i] + '\''
    sql = 'select Species.strain, Species.strain_id from Species join Species_Polymer on Species_Polymer.strain_id = Species.strain_id join Polymer_Data on Polymer_Data.GI = Species_Polymer.GI and Polymer_Data.nomgd_id = Species_Polymer.nomgd_id join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id join Polymer_Alignments on Polymer_Alignments.PData_id = Polymer_Data.PData_id join Alignment on Alignment.Aln_id = Polymer_Alignments.Aln_id where Nomenclature.MoleculeGroup in (' + concatenatedMoleculeGroups + ') and Alignment.name in (' + concatenatedAlignmentNames + ');'
    with connection.cursor() as cursor:
        cursor.execute(sql)
        results = []
        for row in cursor.fetchall():
            results.append([row[0], row[1]])
    context = {
        'results' : results
    }
    return JsonResponse(context)

def getAlignmentFilterByNameAndMoleculeGroupTrimByStrainIdDirect(request, concatenatedMoleculeGroups, concatenatedAlignmentNames, concatenatedStrainIds):
    moleculeGroups = concatenatedMoleculeGroups.split(',')
    concatenatedMoleculeGroups = '\'' + moleculeGroups[0] + '\''
    for i in range(1, len(moleculeGroups)):
        concatenatedMoleculeGroups += ', \'' + moleculeGroups[i] + '\''
    alignmentNames = concatenatedAlignmentNames.split(',')
    concatenatedAlignmentNames = '\'' + alignmentNames[0] + '\''
    for i in range(1, len(alignmentNames)):
        concatenatedAlignmentNames += ', \'' + alignmentNames[i] + '\''
    sql = 'select Polymer_metadata.Fullseq, Species.strain from Polymer_metadata join Polymer_Data on Polymer_metadata.polymer_id = Polymer_Data.PData_id join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id join Species_Polymer on Polymer_Data.GI = Species_Polymer.GI and Polymer_Data.nomgd_id = Species_Polymer.nomgd_id join Species on Species.strain_id = Species_Polymer.strain_id join Polymer_Alignments on Polymer_Alignments.PData_id = Polymer_Data.PData_id join Alignment on Polymer_Alignments.Aln_id = Alignment.Aln_id where Nomenclature.MoleculeGroup in (' + concatenatedMoleculeGroups + ') and Alignment.Name in (' + concatenatedAlignmentNames + ') and Species.strain_id in (' + concatenatedStrainIds + ');'
    with connection.cursor() as cursor:
        cursor.execute(sql)
        results = []
        for row in cursor.fetchall():
            results.append(row[0])
    context = {
        'results' : results
    }
    return JsonResponse(context)

def string_fasta(request, protein_type, aln_name, tax_group, internal=False):
    if type(tax_group) == int:
        tax_group = str(tax_group)
    elif type(tax_group) == list:
        tax_group = ','.join(tax_group)
    protein_type = protein_type.split(',')
    for i in range(len(protein_type)):
        protein_type[i] = '\'' + protein_type[i] + '\''
    protein_type=','.join(protein_type)

    with connection.cursor() as cursor:
        sql = 'SET sql_mode=(SELECT REPLACE(@@sql_mode,\'ONLY_FULL_GROUP_BY\',\'\'));'
        cursor.execute(sql)
        sql = "select Alignment.*, Polymer_Data.*, Nomenclature.* from Alignment join Polymer_Alignments on Polymer_Alignments.Aln_id = Alignment.Aln_id join Polymer_Data on Polymer_Data.PData_id = Polymer_Alignments.PData_id join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id join Species on Polymer_Data.strain_id = Species.strain_id join Species_TaxGroup on Species_TaxGroup.strain_id = Species.strain_id join TaxGroups on Species_TaxGroup.taxgroup_id = Species_TaxGroup.taxgroup_id where Alignment.Name = '" + aln_name + "' and MoleculeGroup in (" + protein_type + ") and TaxGroups.taxgroup_id in (" + tax_group + ") group by Alignment.Name;"
        cursor.execute(sql)
        raw_result = aqab.dictfetchall(cursor)
    return simple_fasta(request, raw_result[0]['Aln_id'], tax_group, internal)

def custom_modified_residues(request, entity_id, cif_file_path):
    import Bio.PDB.MMCIF2Dict
    mmcdata = Bio.PDB.MMCIF2Dict.MMCIF2Dict("/tmp/" + cif_file_path)
    index = 0
    try:
        index = mmcdata['_entity_poly.pdbx_strand_id'].index(entity_id)
    except:
        i = 0
        for item in mmcdata['_entity_poly.pdbx_strand_id']:
            if entity_id in item.split(','):
                index = i
            i += 1
    pattern = '\([a-zA-Z0-9]*\)'
    sequence = mmcdata['_entity_poly.pdbx_seq_one_letter_code'][index]
    sequence_without_spaces = ''.join(sequence.split())
    #modiifed_residues = re.findall('\([a-zA-Z0-9]*\)', mmcdata['_entity_poly.pdbx_seq_one_letter_code'][0])
    #iter = re.finditer('\([a-zA-Z0-9]*\)', mmcdata['_entity_poly.pdbx_seq_one_letter_code'][0])
    modified_residues = []
    for match in re.finditer(pattern, sequence_without_spaces):
        s = match.start() + 1
        e = match.end() + 1
        modified_residues.append([sequence_without_spaces[s:e - 2], s, e])
    #indices = [m.start(0) for m in iter]
    context = {
        'Modified' : modified_residues
    }
    return JsonResponse(context)

def modified_residues(request, pdbid, chain_id):
    import Bio.PDB.MMCIF2Dict
    mmcdata = Bio.PDB.MMCIF2Dict.MMCIF2Dict("/tmp/PDB/" + str(pdbid) + ".cif")
    index = 0
    try:
        index = mmcdata['_entity_poly.pdbx_strand_id'].index(chain_id)
    except:
        i = 0
        for item in mmcdata['_entity_poly.pdbx_strand_id']:
            if chain_id in item.split(','):
                index = i
            i += 1
    pattern = '\([a-zA-Z0-9]*\)'
    sequence = mmcdata['_entity_poly.pdbx_seq_one_letter_code'][index]
    sequence_without_spaces = ''.join(sequence.split())
    #modiifed_residues = re.findall('\([a-zA-Z0-9]*\)', mmcdata['_entity_poly.pdbx_seq_one_letter_code'][0])
    #iter = re.finditer('\([a-zA-Z0-9]*\)', mmcdata['_entity_poly.pdbx_seq_one_letter_code'][0])
    modified_residues = []
    for match in re.finditer(pattern, sequence_without_spaces):
        s = match.start() + 1
        e = match.end() + 1
        modified_residues.append([sequence_without_spaces[s:e - 2], s, e])
    #indices = [m.start(0) for m in iter]
    context = {
        'Modified' : modified_residues
    }
    return JsonResponse(context)

def protein_contacts(request, pdbid, chain_id):
    pdbl = PDB.PDBList(pdb='/tmp/PDB')
    pdbl.retrieve_pdb_file(pdbid, pdir='/tmp/PDB')
    parser = PDB.MMCIFParser()
    structure = parser.get_structure(pdbid, "/tmp/PDB/" + str(pdbid) + ".cif")
    atom_list_23S = PDB.Selection.unfold_entities(structure[0][chain_id], 'A')
    neighbor_23S = PDB.NeighborSearch(atom_list_23S)
    neighbors = {}
    for chain in structure[0]:
        isProtein = False
        for child in chain.child_list:
            if child.resname != 'A' and child.resname != 'U' and child.resname != 'C' and child.resname != 'G':
                isProtein = True
                break
            elif child.resname == 'U':
                break
        if isProtein:
            target_atoms_L2=list(at for at in structure[0][chain.id].get_atoms() if at.parent.id[0] == ' ')
            neighbors_L2_all=set()
            for atom in target_atoms_L2:
                point = atom.get_coord()
                neighbors_L2 = neighbor_23S.search(point, 3.5, level='R')
                for neighbor in neighbors_L2:
                    if not neighbor.id[0] == 'W':
                        neighbors_L2_all.add(neighbor.id[1])

            if len(neighbors_L2_all) > 0:
                neighbors[chain.id] = list(neighbors_L2_all)
    
    return JsonResponse(neighbors)
def full_RNA_seq(request, pdbid, chain_id):
    import alignments.config
    
    RNA_full_sequence={}
    print('cid',chain_id)
    cif_fileNameSuffix=alignments.config.cif_fileNameSuffix_share
    print(alignments.config.cif_fileNameSuffix_share)
    import Bio.PDB.MMCIF2Dict
    #mmcdata = Bio.PDB.MMCIF2Dict.MMCIF2Dict('/tmp/cust2{cif_fileNameSuffix}.cif')
    mmcdata = Bio.PDB.MMCIF2Dict.MMCIF2Dict("/tmp/cust2"+str(cif_fileNameSuffix)+".cif")
    #print(mmcdata)
    index = 0
    try:
        index = mmcdata['_entity_poly.pdbx_strand_id'].index(chain_id)
        
    except:
        i = 0
        for item in mmcdata['_entity_poly.pdbx_strand_id']:
            if chain_id in item.split(','):
                index = i
            i += 1
    pattern = '\([a-zA-Z0-9]*\)'
    sequence = mmcdata['_entity_poly.pdbx_seq_one_letter_code_can'][index]
    #FULL RNA SEQ IS PARSED AND STORED HERE
    RNA_full_sequence = ''.join(sequence.split())
    context = {
        'RNAseq' : RNA_full_sequence
    }
    
    return JsonResponse(context)
file_r2dt_counter_dict={}    
def r2dt(request, entity_id):
    import os
    import datetime
    import alignments.config
    cwd = os.getcwd()
    now = datetime.datetime.now()
    sequence_file = request.FILES["custom_seq_file"]
    as_bytes = sequence_file.read()
    sequence_file_text = as_bytes.decode().replace("\r", "")
    sequence_file_lines = sequence_file_text.split("\n")
    keys = json.loads(sequence_file_lines[0])
    sequence = "\n".join(sequence_file_lines[1:])
    
    RIBODIR=os.environ['RIBODIR']
    os.chdir('/home/sumon/repos/Ribovision_2.0_GT/rna/R2DT-master')
    fileNameSuffix = "_" + str(now.year) + "_" + str(now.month) + "_" + str(now.day) + "_" + str(now.hour) + "_" + str(now.minute) + "_" + str(now.second) + "_" + str(now.microsecond)
    if request.method == "POST":
        
        cif_mode_flag = keys["cif_mode_flag"]#request.POST["cif_mode_flag"]
        parsed_cif_mode_flag = cif_mode_flag
        if cif_mode_flag == "true":
            parsed_cif_mode_flag = True
        elif cif_mode_flag == "false":
            parsed_cif_mode_flag = False
        elif cif_mode_flag == "":
            parsed_cif_mode_flag = None
        cif_mode_flag = parsed_cif_mode_flag
    else:
        cif_mode_flag = None

    newcwd = os.getcwd()
    with open('sequence10'+str(fileNameSuffix)+'.fasta', 'w') as f:
        f.write('>Sequence\n')
        f.write(sequence)
        f.close()       
    output = f"/home/sumon/repos/Ribovision_2.0_GT/rna/R2DT-master/R2DT-test20{fileNameSuffix}"    
     
    
    cmd = f'python3 r2dt.py draw  {newcwd}/sequence10{fileNameSuffix}.fasta {output}' # --skip_ribovore_filters 
    os.system(cmd)
    
    filename = '' 
    print('output_dir', output)
    
    files = os.listdir(os.path.join(output, "results/json"))
    
    print('json_file',  files)
    
    if len(files) == 0:
        print('reached here')
        # clean up
        shutil.rmtree(output)
        
        cmd = f'python3 r2dt.py draw --skip_ribovore_filters {newcwd}/sequence10{fileNameSuffix}.fasta {output}' 
        os.system(cmd)
        
    files = os.listdir(os.path.join(output, "results/json"))
    filename = os.path.join(output, "results/json", files[0])
    
    files_to_remove = []

    if (cif_mode_flag is None) or (not cif_mode_flag):
        #FOR NONE OR PDB modes
        #cmd = f'/usr/bin/python3 {newcwd}/json2json_split2.py -i {filename} -o1 {output}/results/json/RNA_2D_json.json -o2 {output}/results/json/BP_json.json'
        cmd = f'/usr/bin/python3 /home/sumon/repos/Ribovision_2.0_GT/rna/R2DT-master/json2json_split2.py -i {filename} -o1 {output}/results/json/RNA_2D_json.json -o2 {output}/results/json/BP_json.json'
    else:
        #FOR CIF MODE
        cif_file_path = keys["cif_file_path"]#request.POST["cif_file_path"]
        #files_to_remove.append(cif_file_path)
        print('r2dt results parsing cif')
        print('cif_file_path')
        cmd = f'python3 /home/sumon/repos/Ribovision_2.0_GT/rna/R2DT-master/parse_cif4.py -ij {filename} -ic {cif_file_path} -ie {entity_id} -o1 {output}/results/json/RNA_2D_json.json -o2 {output}/results/json/BP_json.json'
        print('r2dt cif parsed')
    os.system(cmd)

    with open(f'{output}/results/json/RNA_2D_json.json', 'r') as f:
        data = json.loads(f.read())
        f.close()
        
    with open(f'{output}/results/json/BP_json.json', 'r') as f:
        BP = json.loads(f.read())
        f.close()   

    r2dt_json = {
        'RNA_2D_json' : data, 
        'RNA_BP_json' : BP
        
    }

    for file_path in files_to_remove:
        if os.path.isfile(file_path):
            os.remove(file_path)

    if os.path.isfile('./sequence10'+str(fileNameSuffix)+'.fasta'):
        os.remove('./sequence10'+str(fileNameSuffix)+'.fasta')
    else:
        # If it fails, inform the user.
        print("Error: %s file not found" % './sequence10'+str(fileNameSuffix)+'.fasta')
        
    if os.path.isfile('./sequence10'+str(fileNameSuffix)+'.fasta.ssi'):
        os.remove('./sequence10'+str(fileNameSuffix)+'.fasta.ssi')
    else:
        # If it fails, inform the user.
        print("Error: %s file not found" % './sequence10'+str(fileNameSuffix)+'.fasta.ssi') 
    
    if cif_mode_flag is True:
        cif_fileNameSuffix=alignments.config.cif_fileNameSuffix_share
        cif_file_name="/tmp/cust2"+str(cif_fileNameSuffix)+".cif"
        if not cif_file_name in file_r2dt_counter_dict:
            file_r2dt_counter_dict[cif_file_name] = 0
        file_r2dt_counter_dict[cif_file_name] += 1

        if file_r2dt_counter_dict[cif_file_name] == 2:  

            # Delete file
            if os.path.isfile(cif_file_name):
                os.remove(cif_file_name)
                file_r2dt_counter_dict[cif_file_name] = 0
            else:
            # If it fails, inform the user.
                print("Error: %s file not found" % cif_file_name)       
    
    dir_path = str(output)
    if os.path.isdir(dir_path):
        # remove directory and all its content
        shutil.rmtree(dir_path)
    else:
        raise ValueError("Path {} is not a file or dir.".format(dir_path))
     
    return JsonResponse(r2dt_json)
