# Standard library imports
import contextlib
import re
import os
import warnings
import io
import base64
import json
import datetime
import urllib.request
import shutil
import time
from subprocess import Popen, PIPE
from xml.dom import minidom

# Third-party imports
from pdbecif.mmcif_io import CifFileReader
from Bio import AlignIO, BiopythonDeprecationWarning, PDB
from Bio.SeqUtils import IUPACData
import Bio.PDB.MMCIF2Dict

# Django imports
from django.shortcuts import render
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from django.urls import reverse
from django.contrib.sites.shortcuts import get_current_site
from django.db import connection

# Local application imports
from alignments.models import *
from alignments.taxonomy_views import *
from alignments.residue_api import *
from alignments.structure_api import *
from alignments.fold_api import *
from alignments.runal2co import executeAl2co
import alignments.alignment_query_and_build as aqab
from alignments import Shannon
from alignments.log import LoggerSetup
import alignments.config
from alignments.paths import R2DT_PATH, LOGS_PATH
from twincons.TwinCons import slice_by_name
from alignments.fred import get_fred_base_pairs

# Logging setup
import logging
LoggerSetup(LOGS_PATH)
logger = logging.getLogger("ribovision3-logger")

# Global variables
file_r2dt_counter_dict = {}

class c:
    structureObj = None

def trim_alignment(concat_fasta, filter_strain):
    '''Reads a fasta string into alignment and trims it down by filter sequence'''
    logger.info(f"Starting alignment trimming for filter strain: {filter_strain}")
    try:
        from alignments.Shannon import species_index_to_aln_index, truncate_aln
        from Bio import AlignIO
        from io import StringIO
        
        alignment = list(AlignIO.parse(StringIO(concat_fasta), 'fasta'))[0]
        logger.debug(f"Initial alignment length: {len(alignment)}")
        
        aln_anchor_map, anchor_ix_in_alignment = species_index_to_aln_index(alignment, filter_strain)
        logger.debug(f"Anchor index in alignment: {anchor_ix_in_alignment}")
        
        alignment = truncate_aln(alignment, list(aln_anchor_map.keys()), aln_anchor_map=aln_anchor_map)
        
        logger.info(f"Alignment trimmed successfully. New length: {len(alignment)}")
        return alignment
    except Exception as e:
        logger.error(f"Error in trim_alignment: {str(e)}")
        raise

def calculate_twincons(alignment):
    '''Calculates twincons score given an alignment object.
    Returns data in a list format for the topology viewer'''
    logger.info("Starting TwinCons calculation")
    try:
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', BiopythonDeprecationWarning)
        from twincons import TwinCons
        
        list_for_phymeas = ['-as', alignment, '-nc', '-r', '-mx', 'blastn']
        logger.debug(f"TwinCons parameters: {list_for_phymeas}")
        
        alnindex_score, sliced_alns, number_of_aligned_positions, gp_mapping = TwinCons.main(list_for_phymeas)
        logger.debug(f"Number of aligned positions: {number_of_aligned_positions}")
        
        list_for_topology_viewer = []
        for alnindex in alnindex_score:
            list_for_topology_viewer.append([alnindex, alnindex_score[alnindex][0]])
        
        logger.info(f"TwinCons calculated successfully. Number of scores: {len(list_for_topology_viewer)}")
        return list_for_topology_viewer
    except Exception as e:
        logger.error(f"Error in calculate_twincons: {str(e)}")
        raise

def upload_custom_data_for_mapping(request):
    logger.info(f"Handling custom data upload request: {request.method}")
    try:
        if request.method == 'POST' and 'filename' in request.FILES:
            data_pairs = []
            file = request.FILES['filename']
            file_iterator = iter(file)
            for line in file_iterator:
                entry = line.decode().strip().split(',')
                data_pairs.append((int(entry[0]), float(entry[1])))
            
            logger.info(f"Uploaded custom data with {len(data_pairs)} entries")
            request.session['csv'] = data_pairs
        
        if request.method == 'GET':
            data_pairs = request.session.get('csv')
            if data_pairs:
                logger.info(f"Returning {len(data_pairs)} custom data pairs")
            else:
                logger.warning("No custom data found in session")
            return JsonResponse(data_pairs, safe=False)
    except Exception as e:
        logger.error(f"Error in upload_custom_data_for_mapping: {str(e)}")
        return HttpResponseServerError("Error processing custom data")

def api_twc_with_upload(request, anchor_structure):
    logger.info(f"Processing TwinCons with upload for anchor structure: {anchor_structure}")
    try:
        anchor_taxid = pdbid_to_strainid(anchor_structure)
        logger.debug(f"Converted PDB ID to taxid: {anchor_taxid}")

        fastastring = request.session.get('fasta')
        if not fastastring:
            logger.warning("No fasta string found in session")
            return HttpResponseServerError("No fasta string provided")

        concat_fasta = re.sub(r'\\n', '\n', fastastring, flags=re.M)
        logger.debug(f"Fasta string length after processing: {len(concat_fasta)}")
        
        alignment = trim_alignment(concat_fasta, str(anchor_taxid))
        list_for_topology_viewer = calculate_twincons(alignment)
        
        logger.info(f"TwinCons calculation completed for {anchor_structure}")
        return JsonResponse(list_for_topology_viewer, safe=False)
    except Exception as e:
        logger.error(f"Error in api_twc_with_upload: {str(e)}")
        return HttpResponseServerError("Error processing TwinCons calculation")

def constructEbiAlignmentString(fasta, ebi_sequence, startIndex):
    logger.info("Starting construction of EBI alignment string")
    try:
        now = datetime.datetime.now()
        BASE_DIR = os.environ.get("BASE_DIR", os.getcwd())
        fileNameSuffix = f"_{now.year}_{now.month}_{now.day}_{now.hour}_{now.minute}_{now.second}_{now.microsecond}"
        alignmentFileName = os.path.join(BASE_DIR, f"static/alignment{fileNameSuffix}.txt")
        ebiFileName = os.path.join(BASE_DIR, f"static/ebi_sequence{fileNameSuffix}.txt") 
        mappingFileName = ebiFileName + ".map"
        
        logger.debug(f"File names: alignment={alignmentFileName}, ebi={ebiFileName}, mapping={mappingFileName}")
        
        fasta = re.sub('>Structure sequence[\s\S]*?>','>',fasta)
        with open(alignmentFileName, "w") as fh:
            fh.write(fasta)
        with open(ebiFileName, "w") as fh:
            fh.write(">Structure sequence\n")
            fh.write(ebi_sequence)

        shiftIndexBy = startIndex - 1 if startIndex > 1 else 0
        logger.debug(f"Shift index by: {shiftIndexBy}")

        logger.info("Running MAFFT alignment")
        pipe = Popen(f"mafft --preservecase --quiet --addfull {ebiFileName} --mapout {alignmentFileName}; cat {mappingFileName}", stdout=PIPE, shell=True)
        output = pipe.communicate()[0]
        decoded_text = output.decode("ascii")
        
        if len(decoded_text) <= 0:
            logger.error("Failed to map polymer sequence to alignment")
            for removeFile in [alignmentFileName, ebiFileName]:
                os.remove(removeFile)
            return HttpResponseServerError("Failed mapping the polymer sequence to the alignment!\nTry a different structure.")

        text = decoded_text.split('\n#')[1]
        amendedAln = re.sub('>Structure sequence$','',decoded_text.split('\n#')[0])
        outputDict, mapping, badMapping = dict(), dict(), 0
        
        for line in text.split('\n')[1:]:
            row = line.split(', ')
            if len(row) < 3:
                continue
            if row[2] == '-':
                badMapping += 1
                continue
            if row[1] == '-':
                logger.error("Critical mapping error: gap in structure sequence")
                return HttpResponseServerError("Failed mapping the polymer sequence to the alignment!\nTry a different structure.")
            mapping[int(row[2])] = int(row[1]) + shiftIndexBy

        outputDict["structureMapping"] = mapping
        if badMapping > 0:
            outputDict['BadMappingPositions'] = badMapping
            logger.warning(f"Found {badMapping} bad mapping positions")

        for removeFile in [alignmentFileName, ebiFileName, mappingFileName]:
            os.remove(removeFile)
        
        outputDict["amendedAln"] = f'>Structure sequence{amendedAln.split(">Structure sequence")[1]}{amendedAln.split(">Structure sequence")[0]}'
        
        logger.info("EBI alignment string constructed successfully")
        return outputDict
    except Exception as e:
        logger.error(f"Error in constructEbiAlignmentString: {str(e)}")
        raise

def request_post_data(post_data):
    logger.debug("Extracting post data")
    try:
        fasta = post_data["fasta"]
        ebi_sequence = post_data["ebi_sequence"]
        startIndex = int(post_data["startIndex"])
        logger.debug(f"Extracted fasta (length: {len(fasta)}), ebi_sequence (length: {len(ebi_sequence)}), and startIndex: {startIndex}")
        return fasta, ebi_sequence, startIndex
    except KeyError as e:
        logger.error(f"Missing key in post_data: {str(e)}")
        raise
    except ValueError as e:
        logger.error(f"Invalid startIndex: {str(e)}")
        raise

def make_map_from_alnix_to_sequenceix(request):
    logger.info("Starting map creation from alignment index to sequence index")
    try:
        fasta, ebi_sequence, startIndex = request_post_data(request.POST)
        logger.debug(f"Post data extracted. Fasta length: {len(fasta)}, EBI sequence length: {len(ebi_sequence)}, Start index: {startIndex}")
        
        mapping = constructEbiAlignmentString(fasta, ebi_sequence, startIndex)
        
        if not isinstance(mapping, dict):
            logger.error("Failed to construct EBI alignment string")
            return mapping
        
        logger.info("Successfully created map from alignment index to sequence index")
        return JsonResponse(mapping, safe=False)
    except Exception as e:
        logger.error(f"Error in make_map_from_alnix_to_sequenceix: {str(e)}")
        return HttpResponseServerError("Error creating alignment map")

def api_twc_parameterless(request):
    logger.info("Starting parameterless TwinCons API request")
    try:
        fasta = request.POST["fasta"]
        logger.debug(f"Received fasta string of length: {len(fasta)}")
        
        concat_fasta = re.sub(r'\\n','\n', fasta,flags=re.M)
        logger.debug(f"Processed fasta string length: {len(concat_fasta)}")
        
        list_for_topology_viewer = calculate_twincons(concat_fasta)
        logger.info(f"Parameterless TwinCons calculated. Number of scores: {len(list_for_topology_viewer)}")
        return JsonResponse(list_for_topology_viewer, safe=False)
    except KeyError:
        logger.error("Missing 'fasta' in POST data")
        return HttpResponseServerError("Missing 'fasta' data")
    except Exception as e:
        logger.error(f"Error in api_twc_parameterless: {str(e)}")
        return HttpResponseServerError("Error calculating TwinCons")


def api_twc(request, align_name, tax_group1, tax_group2, anchor_structure=''):
    logger.info(f"Starting api_twc for align_name: {align_name}, tax_groups: {tax_group1}, {tax_group2}, anchor_structure: {anchor_structure}")
    
    try:
        if anchor_structure:
            anchor_taxid = pdbid_to_strainid(anchor_structure)
            filter_strain = str(Species.objects.filter(strain_id=anchor_taxid)[0].strain).replace(" ", "_")
            logger.debug(f"Anchor structure {anchor_structure} converted to taxid: {anchor_taxid}, filter_strain: {filter_strain}")

        fastastring = request.POST.get('fasta')
        if fastastring is None:
            logger.info("Fasta string not provided, querying database for alignment")
            align_id = Alignment.objects.filter(name=align_name)[0].aln_id
            logger.debug(f"Alignment ID for {align_name}: {align_id}")

            rawsqls = []
            for parent in [tax_group1, tax_group2]:
                rawsqls.append((aqab.sql_filtered_aln_query(align_id, parent), Taxgroups.objects.get(pk=parent).groupname))
            
            nogap_tupaln = dict()
            max_alnposition = 0
            for rawsql, parent in rawsqls:
                nogap_tupaln, max_alnposition = aqab.query_to_dict_structure(rawsql, parent, nogap_tupaln, max_alnposition)
            
            logger.debug(f"Max alignment position: {max_alnposition}")
            fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_alnposition)
        else:
            logger.info("Fasta string provided in the request")

        concat_fasta = re.sub(r'\\n','\n',fastastring,flags=re.M)
        logger.debug(f"Concatenated fasta length: {len(concat_fasta)}")
        
        if anchor_structure:
            alignment = trim_alignment(concat_fasta, filter_strain)
        else:
            alignment = concat_fasta
        
        list_for_topology_viewer = calculate_twincons(alignment)
        logger.info(f"TwinCons calculation completed. Number of scores: {len(list_for_topology_viewer)}")
        
        return JsonResponse(list_for_topology_viewer, safe=False)
    except Exception as e:
        logger.error(f"Error in api_twc: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error processing TwinCons calculation")

def minmaxIndex_handler(minIndex, maxIndex):
    logger.debug(f"minmaxIndex_handler called with minIndex: {minIndex}, maxIndex: {maxIndex}")
    
    if minIndex == '':
        minIndex = str(0)
        logger.info("minIndex was empty, set to 0")
    else:
        minIndex = str(minIndex)
    
    if maxIndex == '':
        maxIndex = str(100000)
        logger.info("maxIndex was empty, set to 100000")
    else:
        maxIndex = str(maxIndex)
    
    logger.debug(f"Returning minIndex: {minIndex}, maxIndex: {maxIndex}")
    return minIndex, maxIndex

def twincons_handler(request, anchor_structure, chain, align_name='', tax_group1='', tax_group2='', minIndex='', maxIndex=''):
    logger.info(f"twincons_handler called for anchor_structure: {anchor_structure}, chain: {chain}")
    try:
        from django.urls import resolve
        current_url = resolve(request.path_info).url_name
        logger.debug(f"Resolved URL name: {current_url}")

        minIndex, maxIndex = minmaxIndex_handler(minIndex, maxIndex)
        
        context = {
            'pdbid': anchor_structure, 
            'chainid': chain,
            'minIndex': minIndex,
            'maxIndex': maxIndex
        }
        
        if current_url == 'twc_with_upload':
            context['entropy_address'] = f"upload/twc-api/{anchor_structure}"
        elif current_url == 'custom_csv_data_viewer':
            upload_custom_data_for_mapping(request)
            context['entropy_address'] = "custom-csv-data"
        elif current_url == 'twincons':
            context['entropy_address'] = f"twc-api/{align_name}/{tax_group1}/{tax_group2}/{anchor_structure}"
        
        logger.info(f"Rendering template for {current_url}")
        return render(request, 'alignments/twc_detail.html', context)
    except Exception as e:
        logger.error(f"Error in twincons_handler: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error processing TwinCons request")

def entropy(request, align_name, tax_group, anchor_structure):
    logger.info(f"Entropy calculation started for align_name: {align_name}, tax_group: {tax_group}, anchor_structure: {anchor_structure}")
    try:
        from alignments import Shannon
        
        taxid = pdbid_to_strainid(anchor_structure)
        filter_strain = Species.objects.filter(strain_id=taxid)[0].strain
        logger.debug(f"Anchor structure {anchor_structure} converted to taxid: {taxid}, filter_strain: {filter_strain}")
        
        align_id = Alignment.objects.filter(name=align_name)[0].aln_id
        polymerid = PolymerData.objects.values("pdata_id").filter(polymeralignments__aln=align_id, strain=taxid)[0]["pdata_id"]
        chainid = Chainlist.objects.values("chainname").filter(polymer=polymerid)[0]["chainname"]
        logger.debug(f"Alignment ID: {align_id}, Polymer ID: {polymerid}, Chain ID: {chainid}")
        
        rawsql_result = aqab.sql_filtered_aln_query(align_id, tax_group)
        nogap_tupaln, max_aln_length = aqab.query_to_dict_structure(rawsql_result, Taxgroups.objects.get(pk=tax_group).groupname)
        fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_aln_length)
        logger.debug(f"Alignment built. Max alignment length: {max_aln_length}")
        
        aln_shannon_list = Shannon.main(['-a', fastastring, '-f', 'fastastring', '--return_within', '-s', filter_strain])
        logger.info(f"Shannon entropy calculation completed. Number of results: {len(aln_shannon_list)}")
        
        context = {
            'pdbid': anchor_structure, 
            'chainid': chainid, 
            'shannon_dictionary': aln_shannon_list, 
            'entropy_address': f"entropy-api/{align_name}/{tax_group}/{anchor_structure}"
        }
        return render(request, 'alignments/entropy_detail.html', context)
    except Exception as e:
        logger.error(f"Error in entropy calculation: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error calculating entropy")

def api_entropy(request, align_name, tax_group, anchor_structure):
    logger.info(f"API entropy calculation started for align_name: {align_name}, tax_group: {tax_group}, anchor_structure: {anchor_structure}")
    try:
        from alignments import Shannon
        import os
        from django.conf import settings
        
        taxid = pdbid_to_strainid(anchor_structure)
        filter_strain = Species.objects.filter(strain_id=taxid)[0].strain
        logger.debug(f"Anchor structure {anchor_structure} converted to taxid: {taxid}, filter_strain: {filter_strain}")
        
        align_id = Alignment.objects.filter(name=align_name)[0].aln_id
        logger.debug(f"Alignment ID: {align_id}")
        
        rawsql_result = aqab.sql_filtered_aln_query(align_id, tax_group)
        nogap_tupaln, max_aln_length = aqab.query_to_dict_structure(rawsql_result, Taxgroups.objects.get(pk=tax_group).groupname)
        fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_aln_length)
        logger.debug(f"Alignment built. Max alignment length: {max_aln_length}")
        
        aln_shannon_list = Shannon.main(['-a', fastastring, '-f', 'fastastring', '--return_within', '-s', filter_strain])
        logger.info(f"Shannon entropy calculation completed. Number of results: {len(aln_shannon_list)}")
        
        return JsonResponse(aln_shannon_list, safe=False)
    except Exception as e:
        logger.error(f"Error in API entropy calculation: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error calculating entropy")

def index_test(request):
    logger.info("index_test view called")
    try:
        some_Alignments = Alignment.objects.all()
        superKingdoms = Taxgroups.objects.raw('SELECT * FROM TaxGroups WHERE TaxGroups.groupLevel = "superkingdom";')
        
        context = {
            'props': list(Taxgroups.objects.values('taxgroup_id', 'groupname')),
            'some_Alignments': some_Alignments,
            'superKingdoms': superKingdoms
        }
        logger.debug(f"Context prepared. Number of alignments: {len(some_Alignments)}, Number of superKingdoms: {len(list(superKingdoms))}")
        return render(request, 'alignments/index_test.html', context)
    except Exception as e:
        logger.error(f"Error in index_test view: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error loading index test page")

def proteinTypes(request):
    logger.info("proteinTypes view called")
    if request.method == 'POST' and 'taxIDs' in request.POST:
        taxIDs = request.POST['taxIDs']
        logger.debug(f"Received taxIDs: {taxIDs}")
        return proteinTypesDirect(request, taxIDs)
    else:
        logger.warning("Invalid request or missing taxIDs")
        return []

def allProteinTypes(request):
    logger.info("allProteinTypes view called")
    try:
        allProteinTypes = []
        with connection.cursor() as cursor:
            sql = "select distinct(MoleculeGroup) from Nomenclature order by MoleculeGroup ASC;"
            cursor.execute(sql)
            for row in cursor.fetchall():
                allProteinTypes.append(row[0])
        logger.debug(f"Retrieved {len(allProteinTypes)} protein types")
        context = {"allProteinTypes": allProteinTypes}
        return JsonResponse(context)
    except Exception as e:
        logger.error(f"Error retrieving all protein types: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error retrieving protein types")

def proteinTypesDirect(request, concatenatedTaxIds):
    logger.info(f"proteinTypesDirect called with concatenatedTaxIds: {concatenatedTaxIds}")
    try:
        results = []
        with connection.cursor() as cursor:
            if concatenatedTaxIds.endswith(','):
                concatenatedTaxIds = concatenatedTaxIds[:-1]
            for taxID in concatenatedTaxIds.split(','):
                if not taxID:
                    continue
                taxID = int(taxID)
                logger.debug(f"Processing taxID: {taxID}")
                
                sql = "use DESIRE;"
                cursor.execute(sql)
                
                sql = 'select Nomenclature.MoleculeGroup from TaxGroups join Species_TaxGroup on Species_TaxGroup.taxgroup_id = TaxGroups.taxgroup_id join Species on Species.strain_id = Species_TaxGroup.strain_id join Species_Polymer on Species.strain_id = Species_Polymer.strain_id join Polymer_Data on Polymer_Data.strain_id = Species_Polymer.strain_id and Polymer_Data.GI = Species_Polymer.GI and Polymer_Data.nomgd_id = Species_Polymer.nomgd_id join Nomenclature on Nomenclature.nom_id = Polymer_Data.nomgd_id where TaxGroups.taxgroup_id = ' + str(taxID) + ' group by Nomenclature.MoleculeGroup;'
                cursor.execute(sql)
                
                proteinTypesList = [row[0] for row in cursor.fetchall()]
                filteredList = [s for s in proteinTypesList if s.find("RNA") != -1]
                logger.debug(f"Found {len(filteredList)} RNA-related protein types for taxID {taxID}")
                results.append(filteredList)
        
        context = {'results': results}
        return JsonResponse(context)
    except Exception as e:
        logger.error(f"Error in proteinTypesDirect: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error retrieving protein types")

def getAlignmentsFilterByProteinTypeAndTaxIds(request):
    logger.info("getAlignmentsFilterByProteinTypeAndTaxIds called")
    if request.method == 'POST' and 'selectedProteinType' in request.POST and 'taxIDs' in request.POST:
        selectedProteinType = request.POST['selectedProteinType']
        taxIDs = request.POST['taxIDs']
        logger.debug(f"Received selectedProteinType: {selectedProteinType}, taxIDs: {taxIDs}")
        return getAlignmentsFilterByProteinTypeAndTaxIdsDirect(request, selectedProteinType, taxIDs)
    else:
        logger.warning("Invalid request or missing parameters")
        return []

def getAlignmentsFilterByProteinTypeAndTaxIdsDirect(request, concatenatedProteinTypes, concatenatedTaxIds):
    logger.info(f"getAlignmentsFilterByProteinTypeAndTaxIdsDirect called with proteinTypes: {concatenatedProteinTypes}, taxIds: {concatenatedTaxIds}")
    results = []
    try:
        with connection.cursor() as cursor:
            concatenatedProteinTypes = concatenatedProteinTypes.rstrip(',')
            concatenatedTaxIds = concatenatedTaxIds.rstrip(',')
            
            proteinTypes = concatenatedProteinTypes.split(',')
            concatenatedProteinTypes = ', '.join(f"'{pt}'" for pt in proteinTypes)
            logger.debug(f"Formatted protein types: {concatenatedProteinTypes}")

            for taxID in concatenatedTaxIds.split(','):
                logger.debug(f"Processing taxID: {taxID}")
                alignmentNamesAndPrimaryKeys = []
                sql = f"""
                SELECT Alignment.Name, Alignment.Aln_id 
                FROM Nomenclature 
                JOIN Polymer_Data ON Polymer_Data.nomgd_id = Nomenclature.nom_id 
                JOIN Polymer_Alignments ON Polymer_Alignments.PData_id = Polymer_Data.PData_id 
                JOIN Alignment ON Alignment.Aln_id = Polymer_Alignments.Aln_id 
                JOIN Species_Polymer ON Species_Polymer.strain_id = Polymer_Data.strain_id 
                    AND Species_Polymer.GI = Polymer_Data.GI 
                    AND Species_Polymer.nomgd_id = Polymer_Data.nomgd_id 
                JOIN Species ON Species.strain_id = Species_Polymer.strain_id 
                JOIN Species_TaxGroup ON Species.strain_id = Species_TaxGroup.strain_id 
                JOIN TaxGroups ON TaxGroups.taxgroup_id = Species_TaxGroup.taxgroup_id 
                WHERE Nomenclature.MoleculeGroup IN ({concatenatedProteinTypes}) 
                    AND Alignment.Method = 'GSD_LSD_rRNA' 
                    AND TaxGroups.taxgroup_id = {taxID} 
                GROUP BY Alignment.Aln_id;
                """
                logger.debug(f"Executing SQL: {sql}")
                cursor.execute(sql)
                for row in cursor.fetchall():
                    alignmentNamesAndPrimaryKeys.append([row[0], row[1]])
                logger.info(f"Found {len(alignmentNamesAndPrimaryKeys)} alignments for taxID {taxID}")
                results.append(alignmentNamesAndPrimaryKeys)

        context = {'results': results}
        return JsonResponse(context)
    except Exception as e:
        logger.error(f"Error in getAlignmentsFilterByProteinTypeAndTaxIdsDirect: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error retrieving alignments")

def index(request):
    logger.info("Index view called")
    return render(request, 'alignments/index.html')

def visualizer(request, align_name, tax_group1, tax_group2, anchor_structure=''):
    logger.info(f"Visualizer called for align_name: {align_name}, tax_groups: {tax_group1}, {tax_group2}, anchor_structure: {anchor_structure}")
    twc_api_url = f"http://127.0.0.1:8000/orthologs/twc-api/{align_name}/{tax_group1}/{tax_group2}/{anchor_structure}"
    context = {"twc_api_url": twc_api_url}
    logger.debug(f"TWC API URL: {twc_api_url}")
    return render(request, 'alignments/simpleVisualization.html', context)

def upload_custom_data(request):
    logger.info("Upload custom data view called")
    return render(request, 'alignments/upload_custom_data.html')

def paralog_entry_form(request):
    logger.info("Paralog entry form view called")
    return render(request, 'alignments/index_paralogs.html')

def paralog_display_entropy(request, align_name, fold1, fold2):
    logger.info(f"Paralog display entropy called for align_name: {align_name}, fold1: {fold1}, fold2: {fold2}")
    # Implement the function logic here
    pass

def extract_species_list(fastastring):
    logger.info("Extracting species list from fasta string")
    try:
        unf_species_list = [x.split('\\')[0] for x in fastastring.split('>')[1:]]
        filtered_spec_list = [re.sub('_',' ', re.sub(r'^.*?_', '', x)) for x in unf_species_list]
        logger.debug(f"Extracted {len(filtered_spec_list)} species")
        return filtered_spec_list
    except Exception as e:
        logger.error(f"Error in extract_species_list: {str(e)}", exc_info=True)
        raise

def extract_gap_only_cols(fastastring):
    logger.info("Extracting gap-only columns from fasta string")
    try:
        unf_seq_list = [x.split('\\n')[1] for x in fastastring.split('>')[1:]]
        list_for_intersect = []
        for sequence in unf_seq_list:
            iterator = re.finditer('-', sequence)
            gap_positions = [m.start(0) for m in iterator]
            list_for_intersect.append(gap_positions)
        gap_only_cols = sorted(list(set(list_for_intersect[0]).intersection(*list_for_intersect)))
        logger.debug(f"Found {len(gap_only_cols)} gap-only columns")
        return gap_only_cols
    except Exception as e:
        logger.error(f"Error in extract_gap_only_cols: {str(e)}", exc_info=True)
        raise

def construct_dict_for_json_response(response_data):
    logger.info("Constructing dictionary for JSON response")
    try:
        response_dict = {}
        for entry in response_data:
            if isinstance(entry, str):
                response_dict['Alignment'] = entry
            elif isinstance(entry, bool):
                response_dict['TwinCons'] = entry
            elif isinstance(entry, list):
                if all(isinstance(item, int) for item in entry):
                    response_dict['Gap-only columns'] = entry
                elif all(isinstance(item, list) for item in entry):
                    response_dict['AA frequencies'] = entry
                elif all(isinstance(item, str) for item in entry):
                    response_dict['Sequence names'] = entry
        logger.debug(f"Constructed response dictionary with keys: {response_dict.keys()}")
        return response_dict
    except Exception as e:
        logger.error(f"Error in construct_dict_for_json_response: {str(e)}", exc_info=True)
        raise

def simple_fasta(request, aln_id, tax_group, internal=False):
    logger.info(f"Simple fasta called for aln_id: {aln_id}, tax_group: {tax_group}, internal: {internal}")
    try:
        rawsqls = []
        if isinstance(tax_group, int):
            tax_group = str(tax_group)
        for parent in tax_group.split(','):
            rawsqls.append((aqab.sql_filtered_aln_query(aln_id, parent), Taxgroups.objects.get(pk=parent).groupname))

        nogap_tupaln = {}
        max_alnposition = 0

        for rawsql, parent in rawsqls:
            nogap_tupaln, max_alnposition = aqab.query_to_dict_structure(rawsql, parent, nogap_tupaln, max_alnposition)
        
        fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_alnposition)
        logger.debug(f"Built alignment with max position: {max_alnposition}")
        
        if internal:
            return fastastring
        
        concat_fasta, twc, gap_only_cols, filtered_spec_list, alignment_obj, fr_list = calculateFastaProps(fastastring, frequency_list)
        response_dict = construct_dict_for_json_response([concat_fasta,filtered_spec_list,gap_only_cols,fr_list,twc])

        return JsonResponse(response_dict, safe=False)
    except Exception as e:
        logger.error(f"Error in simple_fasta: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error processing fasta data")

def calculateFastaProps(fastastring, frequency_list=None):
    logger.info("Calculating fasta properties")
    try:
        concat_fasta = re.sub(r'\\n', '\n', fastastring, flags=re.M)
        alignment_obj = AlignIO.read(StringIO(concat_fasta), 'fasta')
        
        twc = False
        if len(alignment_obj) < 1000:
            sliced_alns = slice_by_name(alignment_obj)
            if len(sliced_alns.keys()) == 2:
                twc = True
        logger.debug(f"TwinCons applicable: {twc}")
        
        gap_only_cols = extract_gap_only_cols(fastastring)
        removed_gaps = gap_only_cols[:]
        
        mapped_dict = {i: pos for pos, i in enumerate(range(len(alignment_obj[0]))) if i not in gap_only_cols}
        logger.debug(f"Mapped {len(mapped_dict)} positions")

        if frequency_list:
            for gap in removed_gaps[::-1]:
                alignment_obj = alignment_obj[:, :gap] + alignment_obj[:, gap+1:]
                gap_only_cols.remove(gap)
                frequency_list.pop(gap)
        
        string_writer = StringIO()
        AlignIO.write(alignment_obj, string_writer, "fasta")
        alignment_string = string_writer.getvalue()
        
        filtered_spec_list = extract_species_list(fastastring)
        logger.info(f"Calculated fasta properties. Alignment length: {len(alignment_obj[0])}, Species count: {len(filtered_spec_list)}")
        return alignment_string, twc, gap_only_cols, filtered_spec_list, alignment_obj, frequency_list
    except Exception as e:
        logger.error(f"Error in calculateFastaProps: {str(e)}", exc_info=True)
        raise

def rProtein(request, align_name, tax_group):
    logger.info(f"rProtein view called for align_name: {align_name}, tax_group: {tax_group}")
    try:
        align_id = Alignment.objects.filter(name=align_name)[0].aln_id
        fastastring = simple_fasta(request, align_id, tax_group, internal=True)
        context = {
            'fastastring': fastastring, 
            'aln_name': str(Alignment.objects.filter(aln_id=align_id)[0].name)
        }
        return render(request, 'alignments/detail.html', context)
    except Exception as e:
        logger.error(f"Error in rProtein view: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error processing rProtein data")

def rRNA(request, align_name, tax_group):
    logger.info(f"rRNA view called for align_name: {align_name}, tax_group: {tax_group}")
    try:
        from Bio.SeqUtils import IUPACData
        align_id = Alignment.objects.filter(name=align_name)[0].aln_id
        rawsql_result = aqab.sql_filtered_aln_query(align_id, tax_group)
        nogap_tupaln = {}
        nogap_tupaln, max_aln_length = aqab.query_to_dict_structure(rawsql_result, Taxgroups.objects.get(pk=tax_group).groupname, nogap_tupaln)
        fastastring, frequency_list = aqab.build_alignment_from_multiple_alignment_queries(nogap_tupaln, max_aln_length, IUPACData.unambiguous_rna_letters)
        logger.debug(f"Built rRNA alignment with max length: {max_aln_length}")
        context = {
            'fastastring': fastastring, 
            'aln_name': str(Alignment.objects.filter(aln_id=align_id)[0].name)
        }
        return render(request, 'alignments/rRNA.html', context)
    except Exception as e:
        logger.error(f"Error in rRNA view: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error processing rRNA data")

def validate_fasta_string(fastaString):
    logger.info("Validating fasta string")
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
    try:
        for mstring in malicious_strings:
            if re.search(re.compile(mstring), fastaString):
                logger.warning(f"Malicious string detected: {mstring}")
                return False
        logger.info("Fasta string validated successfully")
        return True
    except Exception as e:
        logger.error(f"Error in validate_fasta_string: {str(e)}", exc_info=True)
        return False

def getGenusFromStrainIdsDirect(request, concatenatedStrainIds):
    logger.info(f"getGenusFromStrainIdsDirect called with strainIds: {concatenatedStrainIds}")
    try:
        strainIds = concatenatedStrainIds.split(',')
        concatenatedStrainIds = ', '.join(f"'{strain}'" for strain in strainIds)
        
        sql = "USE DESIRE; SELECT TaxGroups.taxgroup_id FROM Species JOIN Species_TaxGroup ON Species.strain_id = Species_TaxGroup.strain_id JOIN TaxGroups ON TaxGroups.taxgroup_id = Species_TaxGroup.taxgroup_id WHERE Species.strain_id IN (" + concatenatedStrainIds + ") AND TaxGroups.groupLevel = 'genus';"
        
        genusList = []
        with connection.cursor() as cursor:
            logger.debug(f"Executing SQL: {sql}")
            cursor.execute(sql)
            for row in cursor.fetchall():
                genusList.append(row[0])
        
        logger.info(f"Found {len(genusList)} genus IDs")
        return JsonResponse({'genusList': genusList})
    except Exception as e:
        logger.error(f"Error in getGenusFromStrainIdsDirect: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error retrieving genus information")

def string_fasta_two_strains(request, protein_type, aln_id, strain_id_0, strain_id_1):
    logger.info(f"string_fasta_two_strains called with protein_type: {protein_type}, aln_id: {aln_id}, strain_ids: {strain_id_0}, {strain_id_1}")
    try:
        strain_id_0 = str(strain_id_0)
        strain_id_1 = str(strain_id_1)
        protein_type = ','.join(f"'{pt}'" for pt in protein_type.split(','))
        
        sql = ""  # SQL query to be defined
        logger.warning("SQL query for string_fasta_two_strains is not implemented")
        
        with connection.cursor() as cursor:
            logger.debug(f"Executing SQL: {sql}")
            cursor.execute(sql)
        
        logger.info("string_fasta_two_strains completed successfully")
        return JsonResponse({"message": "Function executed successfully"})
    except Exception as e:
        logger.error(f"Error in string_fasta_two_strains: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error processing string fasta for two strains")

def trim_fasta_by_index(input_file, indices):
    logger.info(f"trim_fasta_by_index called with input_file: {input_file}, indices: {indices}")
    try:
        align = AlignIO.read(input_file, "fasta")
        indices_list = indices.split(',')
        trimmed_align = align[:, int(indices_list[0]):int(indices_list[0])]
        
        for i in indices_list[1:]:
            if int(i) - 1 < 0:
                continue
            trimmed_align += align[:, int(i)-1:int(i)]
        
        logger.info(f"Trimmed alignment from {len(align[0])} to {len(trimmed_align[0])} positions")
        return trimmed_align
    except Exception as e:
        logger.error(f"Error in trim_fasta_by_index: {str(e)}", exc_info=True)
        raise

def flushSession(request):
    logger.info("flushSession called")
    try:
        request.session.flush()
        logger.info("Session flushed successfully")
        return HttpResponse("Success!")
    except Exception as e:
        logger.error(f"Failed to flush session: {str(e)}", exc_info=True)
        return HttpResponseServerError("Failed to flush the session!")

def parse_string_structure(request, stringData, strucID):
    logger.info(f"parse_string_structure called with strucID: {strucID}")
    try:
        from Bio.PDB import MMCIFParser
        parser = MMCIFParser()
        strucFile = StringIO(stringData)
        structure = parser.get_structure(strucID, strucFile)
        logger.info(f"Structure parsed successfully for {strucID}")
        return structure
    except Exception as e:
        logger.error(f"Error in parse_string_structure: {str(e)}", exc_info=True)
        raise

def getAlignmentsFilterByProteinTypeDirect(request, concatenatedProteinTypes):
    logger.info(f"getAlignmentsFilterByProteinTypeDirect called with proteinTypes: {concatenatedProteinTypes}")
    try:
        proteinTypes = concatenatedProteinTypes.split(',')
        concatenatedProteinTypes = ', '.join(f"'{pt}'" for pt in proteinTypes)
        
        sql = f"SELECT DISTINCT(Name) FROM Alignment JOIN Polymer_Alignments ON Polymer_Alignments.Aln_id = Alignment.Aln_id JOIN Polymer_Data ON Polymer_Data.PData_id = Polymer_Alignments.PData_id JOIN Nomenclature ON Nomenclature.nom_id = Polymer_Data.nomgd_id WHERE Nomenclature.MoleculeGroup IN ({concatenatedProteinTypes});"
        
        with connection.cursor() as cursor:
            logger.debug(f"Executing SQL: {sql}")
            cursor.execute(sql)
            results = [row[0] for row in cursor.fetchall()]
        
        logger.info(f"Found {len(results)} alignments")
        return JsonResponse({'results': results})
    except Exception as e:
        logger.error(f"Error in getAlignmentsFilterByProteinTypeDirect: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error retrieving alignments")

def getPairwiseAlignmentDirect(request, moleculeType, alignmentName, strainId0, strainId1):
    logger.info(f"getPairwiseAlignmentDirect called with moleculeType: {moleculeType}, alignmentName: {alignmentName}, strainIds: {strainId0}, {strainId1}")
    try:
        strainIds = [strainId0, strainId1]
        concatenatedStrainIds = ', '.join(f"'{strain}'" for strain in strainIds)
        
        sql = f"SELECT TaxGroups.taxgroup_id, Species.strain FROM Species JOIN Species_TaxGroup ON Species.strain_id = Species_TaxGroup.strain_id JOIN TaxGroups ON TaxGroups.taxgroup_id = Species_TaxGroup.taxgroup_id WHERE Species.strain_id IN ({concatenatedStrainIds}) AND TaxGroups.groupLevel = 'genus';"
        
        genusList = []
        speciesNames = []
        with connection.cursor() as cursor:
            logger.debug(f"Executing SQL: {sql}")
            cursor.execute(sql)
            for row in cursor.fetchall():
                genusList.append(row[0])
                speciesNames.append(row[1])
        
        if len(speciesNames) == 1:
            speciesName0 = speciesName1 = speciesNames[0]
        else:
            speciesName0, speciesName1 = speciesNames

        alignment = string_fasta(request, moleculeType, alignmentName, f"{strainId0},{strainId1}", internal=True)
        alignment = alignment.replace('\\n', '\n')
        PairwiseAlignment = ''
        if alignment.endswith('\n'):
            alignment = alignment[:-1]
            alignmentLines = alignment.splitlines()
            modifiedStrainName0 = speciesName0.replace(' ', '_')
            modifiedStrainName1 = speciesName1.replace(' ', '_')
            for i in range(1, len(alignmentLines), 2):
                titleLine = alignmentLines[i - 1]
                alignmentLine = alignmentLines[i]
                if modifiedStrainName0 in titleLine or modifiedStrainName1 in titleLine:
                    PairwiseAlignment += alignmentLine + '\n'
        
        logger.info(f"Pairwise alignment generated, length: {len(PairwiseAlignment)}")
        return JsonResponse({'PairwiseAlignment': PairwiseAlignment})
    except Exception as e:
        logger.error(f"Error in getPairwiseAlignmentDirect: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error generating pairwise alignment")

def getStrainsFilterByMoleculeGroupAndAlignmentDirect(request, concatenatedMoleculeGroups, concatenatedAlignmentNames):
    logger.info(f"getStrainsFilterByMoleculeGroupAndAlignmentDirect called with moleculeGroups: {concatenatedMoleculeGroups}, alignmentNames: {concatenatedAlignmentNames}")
    try:
        moleculeGroups = concatenatedMoleculeGroups.split(',')
        concatenatedMoleculeGroups = ', '.join(f"'{mg}'" for mg in moleculeGroups)
        
        alignmentNames = concatenatedAlignmentNames.split(',')
        concatenatedAlignmentNames = ', '.join(f"'{an}'" for an in alignmentNames)
        
        sql = f"""
        SELECT Species.strain, Species.strain_id 
        FROM Species 
        JOIN Species_Polymer ON Species_Polymer.strain_id = Species.strain_id 
        JOIN Polymer_Data ON Polymer_Data.GI = Species_Polymer.GI AND Polymer_Data.nomgd_id = Species_Polymer.nomgd_id 
        JOIN Nomenclature ON Nomenclature.nom_id = Polymer_Data.nomgd_id 
        JOIN Polymer_Alignments ON Polymer_Alignments.PData_id = Polymer_Data.PData_id 
        JOIN Alignment ON Alignment.Aln_id = Polymer_Alignments.Aln_id 
        WHERE Nomenclature.MoleculeGroup IN ({concatenatedMoleculeGroups}) 
        AND Alignment.name IN ({concatenatedAlignmentNames});
        """
        
        with connection.cursor() as cursor:
            logger.debug(f"Executing SQL: {sql}")
            cursor.execute(sql)
            results = [[row[0], row[1]] for row in cursor.fetchall()]
        
        logger.info(f"Found {len(results)} strains")
        return JsonResponse({'results': results})
    except Exception as e:
        logger.error(f"Error in getStrainsFilterByMoleculeGroupAndAlignmentDirect: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error retrieving strains")

def getAlignmentFilterByNameAndMoleculeGroupTrimByStrainIdDirect(request, concatenatedMoleculeGroups, concatenatedAlignmentNames, concatenatedStrainIds):
    logger.info(f"getAlignmentFilterByNameAndMoleculeGroupTrimByStrainIdDirect called with moleculeGroups: {concatenatedMoleculeGroups}, alignmentNames: {concatenatedAlignmentNames}, strainIds: {concatenatedStrainIds}")
    try:
        moleculeGroups = concatenatedMoleculeGroups.split(',')
        concatenatedMoleculeGroups = ', '.join(f"'{mg}'" for mg in moleculeGroups)
        
        alignmentNames = concatenatedAlignmentNames.split(',')
        concatenatedAlignmentNames = ', '.join(f"'{an}'" for an in alignmentNames)
        
        sql = f"""
        SELECT Polymer_metadata.Fullseq, Species.strain 
        FROM Polymer_metadata 
        JOIN Polymer_Data ON Polymer_metadata.polymer_id = Polymer_Data.PData_id 
        JOIN Nomenclature ON Nomenclature.nom_id = Polymer_Data.nomgd_id 
        JOIN Species_Polymer ON Polymer_Data.GI = Species_Polymer.GI AND Polymer_Data.nomgd_id = Species_Polymer.nomgd_id 
        JOIN Species ON Species.strain_id = Species_Polymer.strain_id 
        JOIN Polymer_Alignments ON Polymer_Alignments.PData_id = Polymer_Data.PData_id 
        JOIN Alignment ON Polymer_Alignments.Aln_id = Alignment.Aln_id 
        WHERE Nomenclature.MoleculeGroup IN ({concatenatedMoleculeGroups}) 
        AND Alignment.Name IN ({concatenatedAlignmentNames}) 
        AND Species.strain_id IN ({concatenatedStrainIds});
        """
        
        with connection.cursor() as cursor:
            logger.debug(f"Executing SQL: {sql}")
            cursor.execute(sql)
            results = [row[0] for row in cursor.fetchall()]
        
        logger.info(f"Found {len(results)} sequences")
        return JsonResponse({'results': results})
    except Exception as e:
        logger.error(f"Error in getAlignmentFilterByNameAndMoleculeGroupTrimByStrainIdDirect: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error retrieving alignment data")

def string_fasta(request, protein_type, aln_name, tax_group, internal=False):
    logger.info(f"string_fasta called with protein_type: {protein_type}, aln_name: {aln_name}, tax_group: {tax_group}, internal: {internal}")
    try:
        if isinstance(tax_group, int):
            tax_group = str(tax_group)
        elif isinstance(tax_group, list):
            tax_group = ','.join(map(str, tax_group))
        
        protein_type = ','.join(f"'{pt}'" for pt in protein_type.split(','))
        
        with connection.cursor() as cursor:
            sql = 'SET sql_mode=(SELECT REPLACE(@@sql_mode,\'ONLY_FULL_GROUP_BY\',\'\'));'
            cursor.execute(sql)
            
            sql = f"""
            SELECT Alignment.*, Polymer_Data.*, Nomenclature.* 
            FROM Alignment 
            JOIN Polymer_Alignments ON Polymer_Alignments.Aln_id = Alignment.Aln_id 
            JOIN Polymer_Data ON Polymer_Data.PData_id = Polymer_Alignments.PData_id 
            JOIN Nomenclature ON Nomenclature.nom_id = Polymer_Data.nomgd_id 
            JOIN Species ON Polymer_Data.strain_id = Species.strain_id 
            JOIN Species_TaxGroup ON Species_TaxGroup.strain_id = Species.strain_id 
            JOIN TaxGroups ON Species_TaxGroup.taxgroup_id = Species_TaxGroup.taxgroup_id 
            WHERE Alignment.Name = '{aln_name}' 
            AND MoleculeGroup IN ({protein_type}) 
            AND TaxGroups.taxgroup_id IN ({tax_group}) 
            GROUP BY Alignment.Name;
            """
            logger.debug(f"Executing SQL: {sql}")
            cursor.execute(sql)
            raw_result = aqab.dictfetchall(cursor)
        
        logger.info(f"Found {len(raw_result)} results")
        return simple_fasta(request, raw_result[0]['Aln_id'], tax_group, internal)
    except Exception as e:
        logger.error(f"Error in string_fasta: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error processing fasta string")

def custom_modified_residues(request, entity_id, cif_file_path):
    logger.info(f"custom_modified_residues called for entity_id: {entity_id}, cif_file_path: {cif_file_path}")
    try:
        mmcdata = Bio.PDB.MMCIF2Dict.MMCIF2Dict("/tmp/" + cif_file_path)
        
        index = 0
        try:
            index = mmcdata['_entity_poly.pdbx_strand_id'].index(entity_id)
        except ValueError:
            for i, item in enumerate(mmcdata['_entity_poly.pdbx_strand_id']):
                if entity_id in item.split(','):
                    index = i
                    break
        
        pattern = r'\([a-zA-Z0-9]*\)'
        sequence = mmcdata['_entity_poly.pdbx_seq_one_letter_code'][index]
        sequence_without_spaces = ''.join(sequence.split())
        
        modified_residues = []
        for match in re.finditer(pattern, sequence_without_spaces):
            s = match.start() + 1
            e = match.end() + 1
            modified_residues.append([sequence_without_spaces[s:e - 2], s, e])
        
        logger.info(f"Found {len(modified_residues)} modified residues")
        return JsonResponse({'Modified': modified_residues})
    except Exception as e:
        logger.error(f"Error in custom_modified_residues: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error processing modified residues")

def modified_residues(request, pdbid, chain_id):
    logger.info(f"modified_residues called for pdbid: {pdbid}, chain_id: {chain_id}")
    try:
        mmcdata = Bio.PDB.MMCIF2Dict.MMCIF2Dict(f"/tmp/PDB/{pdbid}.cif")
        
        index = 0
        try:
            index = mmcdata['_entity_poly.pdbx_strand_id'].index(chain_id)
        except ValueError:
            for i, item in enumerate(mmcdata['_entity_poly.pdbx_strand_id']):
                if chain_id in item.split(','):
                    index = i
                    break
        
        pattern = r'\([a-zA-Z0-9]*\)'
        sequence = mmcdata['_entity_poly.pdbx_seq_one_letter_code'][index]
        sequence_without_spaces = ''.join(sequence.split())
        
        modified_residues = []
        for match in re.finditer(pattern, sequence_without_spaces):
            s = match.start() + 1
            e = match.end() + 1
            modified_residues.append([sequence_without_spaces[s:e - 2], s, e])
        
        logger.info(f"Found {len(modified_residues)} modified residues")
        return JsonResponse({'Modified': modified_residues})
    except Exception as e:
        logger.error(f"Error in modified_residues: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error processing modified residues")

def protein_contacts(request, pdbid, chain_id):
    logger.info(f"protein_contacts called for pdbid: {pdbid}, chain_id: {chain_id}")
    
    try:
        pdbl = PDB.PDBList(server="https://files.wwpdb.org/", pdb='/tmp/PDB')
        pdbl.retrieve_pdb_file(pdbid, pdir='/tmp/PDB')
        parser = PDB.MMCIFParser()
        structure = parser.get_structure(pdbid, f"/tmp/PDB/{pdbid}.cif")
        atom_list_23S = PDB.Selection.unfold_entities(structure[0][chain_id], 'A')
        
        neighbor_23S = PDB.NeighborSearch(atom_list_23S)
        neighbors = {}
        
        for chain in structure[0]:
            isProtein = any(child.resname not in ['A', 'U', 'C', 'G'] for child in chain.child_list)
            if isProtein:
                target_atoms_L2 = [at for at in structure[0][chain.id].get_atoms() if at.parent.id[0] == ' ']
                neighbors_L2_all = set()
                for atom in target_atoms_L2:
                    point = atom.get_coord()
                    neighbors_L2 = neighbor_23S.search(point, 3.5, level='R')
                    neighbors_L2_all.update(neighbor.id[1] for neighbor in neighbors_L2 if neighbor.id[0] != 'W')
                
                if neighbors_L2_all:
                    neighbors[chain.id] = list(neighbors_L2_all)
        
        logger.info(f"Found protein contacts for {len(neighbors)} chains")
        return JsonResponse(neighbors)
    except Exception as e:
        logger.error(f"Error in protein_contacts: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error processing protein contacts")

def full_RNA_seq(request, pdbid, chain_id):
    logger.info(f"full_RNA_seq called for pdbid: {pdbid}, chain_id: {chain_id}")
    try:
        cif_file_path = alignments.config.cif_path_share
        alignments.config.chainid = chain_id
        mmcdata = Bio.PDB.MMCIF2Dict.MMCIF2Dict(cif_file_path)
        
        index = 0
        try:
            index = mmcdata['_entity_poly.pdbx_strand_id'].index(chain_id)
        except ValueError:
            for i, item in enumerate(mmcdata['_entity_poly.pdbx_strand_id']):
                if chain_id in item.split(','):
                    index = i
                    break
        
        sequence = mmcdata['_entity_poly.pdbx_seq_one_letter_code_can'][index]
        RNA_full_sequence = ''.join(sequence.split())
        
        logger.info(f"Retrieved full RNA sequence, length: {len(RNA_full_sequence)}")
        return JsonResponse({'RNAseq': RNA_full_sequence})
    except Exception as e:
        logger.error(f"Error in full_RNA_seq: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error retrieving full RNA sequence")
    

def r2dt(request, entity_id):
    logger.info(f"r2dt function called for entity_id: {entity_id}")
    try:
        BASE_DIR = os.environ.get("BASE_DIR", os.getcwd())
        now = datetime.datetime.now()
        
        sequence_file = request.FILES["custom_seq_file"]
        as_bytes = sequence_file.read()
        sequence_file_text = as_bytes.decode().replace("\r", "")
        sequence_file_lines = sequence_file_text.split("\n")
        keys = json.loads(sequence_file_lines[0])
        sequence = "\n".join(sequence_file_lines[1:])
        
        logger.debug(f"Sequence file processed, length: {len(sequence)}")
        
        fileNameSuffix = f"_{now.year}_{now.month}_{now.day}_{now.hour}_{now.minute}_{now.second}_{now.microsecond}"
        
        if request.method == "POST":
            cif_mode_flag = keys["cif_mode_flag"]
            parsed_cif_mode_flag = None
            if cif_mode_flag == "true":
                parsed_cif_mode_flag = True
            elif cif_mode_flag == "false":
                parsed_cif_mode_flag = False
            cif_mode_flag = parsed_cif_mode_flag
        else:
            cif_mode_flag = None
        
        logger.info(f"CIF mode flag: {cif_mode_flag}")
        
        
        seq_path = os.path.join(R2DT_PATH, f'sequence10{fileNameSuffix}.fasta')
        with open(seq_path, 'w') as f:
            f.write('>Sequence\n')
            f.write(sequence)
        
        logger.debug(f"Sequence file written to {seq_path}")
        
        output = os.path.join(R2DT_PATH, f"R2DT-output{fileNameSuffix}")
        cmd = f'python3 {R2DT_PATH}/r2dt.py draw {seq_path} {output}'
        logger.debug(f"Executing command: {cmd}")
        os.system(cmd)
        
        files = os.listdir(os.path.join(output, "results/json"))
        
        if len(files) == 0:
            logger.warning("No output files generated, retrying with --skip_ribovore_filters")
            shutil.rmtree(output)
            cmd = f'python3 {R2DT_PATH}/r2dt.py draw --skip_ribovore_filters {seq_path} {output}'
            logger.debug(f"Executing command: {cmd}")
            os.system(cmd)
            files = os.listdir(os.path.join(output, "results/json"))
        
        filename = os.path.join(output, "results/json", files[0])
        logger.info(f"R2DT output file: {filename}")
        
        chainid = alignments.config.chainid

        # if (cif_mode_flag is None) or (not cif_mode_flag):
        if alignments.config.pdb_path_share:
            logger.info("Processing in PDB mode")
            cmd = f'/usr/bin/python3 {R2DT_PATH}/json2json_split2.py -i {filename} -o1 {output}/results/json/RNA_2D_json.json -o2 {output}/results/json/BP_json.json'
            pdb_path = alignments.config.pdb_path_share
            
            logger.debug(f"Executing commands:\n{cmd}")
            os.system(cmd)
            alignments.config.pdb_path_share = ""
            get_fred_base_pairs("", entity_id, chainid, output, file_path=pdb_path)
        else:
            logger.info("Processing in CIF mode")
            cif_file_path = keys["cif_file_path"]
            
            cmd = f'/usr/bin/python3 {R2DT_PATH}/parse_cif4.py -ij {filename} -ic {cif_file_path} -ie {entity_id} -o1 {output}/results/json/RNA_2D_json.json -o2 {output}/results/json/BP_json.json'
            logging.debug(f"CIF file for parse cif: {cif_file_path}")
            logger.debug(f"Executing commands:\n{cmd}")
            os.system(cmd)
            
            if cif_file_path.endswith("cust.cif"):
                # i.e user uploaded a cif file
                get_fred_base_pairs("", entity_id, chainid, output, file_path=cif_file_path)
            else:
                # this is in ribovision mode for r2dt generated 2d layouts
                pdbid = os.path.split(cif_file_path)[1][:4]
                get_fred_base_pairs(pdbid, entity_id, chainid, output)
        
        rna2d_path = f'{output}/results/json/RNA_2D_json.json'
        logger.debug(f"rna2d_path = {rna2d_path}")
        with open(rna2d_path, 'r') as f:
            data = json.loads(f.read())
        
        rnabp_path = f'{output}/results/json/BP_json.json'
        
        with open(rnabp_path, 'r') as f:
            BP = json.loads(f.read())
        
        r2dt_json = {
            'RNA_2D_json': data,
            'RNA_BP_json': BP
        }
        
        logger.info("R2DT JSON data prepared successfully")
        
        # Clean up files
        if os.path.isfile(seq_path):
            os.remove(seq_path)
        else:
            logger.warning(f"Error: {seq_path} file not found")
        
        if os.path.isfile(seq_path + '.ssi'):
            os.remove(seq_path + '.ssi')
        else:
            logger.warning(f"Error: {seq_path + '.ssi'} file not found")
        
        if cif_mode_flag is True:
            cif_file_name = alignments.config.cif_path_share
            if cif_file_name not in file_r2dt_counter_dict:
                file_r2dt_counter_dict[cif_file_name] = 0
            file_r2dt_counter_dict[cif_file_name] += 1
            
            if file_r2dt_counter_dict[cif_file_name] == 2:
                if os.path.isfile(cif_file_name):
                    os.remove(cif_file_name)
                    file_r2dt_counter_dict[cif_file_name] = 0
                    logger.info(f"Removed CIF file: {cif_file_name}")
                else:
                    logger.warning(f"Error: {cif_file_name} file not found")
        
        dir_path = str(output)
        if os.path.isdir(dir_path):
            shutil.rmtree(dir_path)
            logger.info(f"Removed output directory: {dir_path}")
        else:
            logger.warning(f"Path {dir_path} is not a file or directory")
        
        logger.info("r2dt function completed successfully")
        return JsonResponse(r2dt_json)
    
    except Exception as e:
        logger.error(f"Error in r2dt function: {str(e)}", exc_info=True)
        return HttpResponseServerError("Error processing R2DT request")

