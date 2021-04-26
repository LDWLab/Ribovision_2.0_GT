from django.http import JsonResponse, HttpResponseServerError
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
from Bio.SeqRecord import SeqRecord
import re

from alignments.views import parse_string_structure

def request_post_data(post_data):
    fasta = post_data["fasta"]
    struc_id = post_data["struc_id"]
    return fasta, struc_id

def make_map_from_alnix_to_sequenceix_new(request):
    fasta, struc_id = request_post_data(request.POST)
    serializeData = request.session[struc_id]
    strucObj = parse_string_structure(serializeData, struc_id)
    seq_ix_mapping, struc_seq, gapsInStruc = constructStrucSeqMap(strucObj)
    mapping = create_aln_struc_mapping_with_mafft(fasta, struc_seq, seq_ix_mapping)
    mapping["gapsInStruc"] = gapsInStruc
    if type(mapping) != dict:
        return mapping
    return JsonResponse(mapping, safe = False)

def constructStrucSeqMap(structure):
    chains = list()
    for chain in structure.get_chains():
        chains.append(chain)
    sequence, gapsInStruc, seq_ix_mapping = str(), list(), dict()
    old_resi, untrue_seq_ix = 0, 1

    residues = list(chains[0].get_residues())
    for resi in residues:
        resi_id = resi.get_id()
        if (old_resi == 0):
            old_resi = resi_id[1]
        if (resi_id[1] - old_resi > 1):
            gapsInStruc.append((old_resi,resi_id[1]))
        if not re.match(r' ', resi_id[2]):
            continue
        if re.match(r'^H_', resi_id[0]):
            continue
        sequence += resi.get_resname().replace(' ','')
        seq_ix_mapping[untrue_seq_ix] = int(resi_id[1])
        untrue_seq_ix += 1
        old_resi = resi_id[1]
    if len(seq1(residues[0].get_resname().replace(' ',''))) != 0:
        sequence = seq1(sequence)

    return seq_ix_mapping, SeqRecord(Seq(sequence)), gapsInStruc

def create_aln_struc_mapping_with_mafft(fasta, struc_seq, seq_ix_mapping):
    from subprocess import Popen, PIPE
    from os import remove, path
    from warnings import warn
    import datetime
    
    fasta = re.sub('>Structure sequence[\s\S]*?>','>',fasta)
    now = datetime.datetime.now()
    fileNameSuffix = "_" + str(now.year) + "_" + str(now.month) + "_" + str(now.day) + "_" + str(now.hour) + "_" + str(now.minute) + "_" + str(now.second) + "_" + str(now.microsecond)
    ### BE CAREFUL WHEN MERGING THE FOLLOWING LINES TO PUBLIC; PATHS ARE HARDCODED FOR THE APACHE SERVER ###
    aln_group_path = "/home/Desire-DEV/PVDev/static/alignment" + fileNameSuffix + ".txt"
    pdb_seq_path = "/home/Desire-DEV/PVDev/static/ebi_sequence" + fileNameSuffix + ".txt"
    mappingFileName = pdb_seq_path + ".map"
    tempfiles = [aln_group_path, pdb_seq_path, mappingFileName]
    for tempf in tempfiles:
        if path.isfile(tempf):
            warn(f"When using mafft to make structural mapping the working directory must be free of file {tempf}. Trying to delete the file.")
            remove(tempf)
            if path.isfile(tempf):
                raise IOError(f"Couldn't delete the file {tempf} please remove it manually!")
    

    fh = open(aln_group_path, "w")
    fh.write(fasta)
    fh.close()

    fh = open(pdb_seq_path, "w")
    fh.write(">Structure sequence\n")
    fh.write(str(struc_seq.seq))
    fh.close()
    
    pipe = Popen(f"mafft --quiet --addfull {pdb_seq_path} --mapout {aln_group_path}; cat {mappingFileName}", stdout=PIPE, shell=True)
    output = pipe.communicate()[0]

    if len(output.decode("ascii")) <= 0:
        for removeFile in tempfiles:
            remove(removeFile)
        return HttpResponseServerError("Failed mapping the polymer sequence to the alignment!\nTry a different structure.")

    mapping_file = output.decode("ascii").split('\n#')[1]
    amendedAln = re.sub('>Structure sequence$','',output.decode("ascii").split('\n#')[0])
    groupName = output.decode('ascii').split('>')[1].split('_')[0]
    firstLine = True
    outputDict, mapping, bad_map_positions, fail_map = dict(), dict(), 0, False
    for line in mapping_file.split('\n'):
        if firstLine:
            firstLine = False
            continue
        row = line.split(', ')
        if len(row) < 3:
            continue
        if row[2] == '-':
            bad_map_positions += 1
            continue
        if row[1] == '-':
            fail_map = True
        mapping[int(row[2])] = seq_ix_mapping[int(row[1])]
    for tempf in tempfiles:
        remove(tempf)
    if fail_map:
        return HttpResponseServerError("Failed mapping the polymer sequence to the alignment!\nTry a different structure.")
    if bad_map_positions > 0:
        outputDict['BadMappingPositions'] = bad_map_positions
    outputDict["amendedAln"] = f'>Structure sequence{amendedAln.split(">Structure sequence")[1]}{amendedAln.split(">Structure sequence")[0]}'
    outputDict["structureMapping"] = mapping
    return outputDict
