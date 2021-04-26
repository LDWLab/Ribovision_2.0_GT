from django.http import JsonResponse, HttpResponseServerError, HttpResponse
from io import StringIO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.SeqUtils import IUPACData

from alignments.views import validate_fasta_string, calculateFastaProps, construct_dict_for_json_response

def handle_custom_upload_alignment(request):
    if request.method == 'POST' and 'custom_aln_file' in request.FILES:
        aln_file = request.FILES['custom_aln_file']
        alignment_string = ''
        for aln_part in aln_file.chunks():
            alignment_string += aln_part.decode()
        try:
            alignments = list(AlignIO.parse(StringIO(alignment_string), 'fasta'))
        except ValueError as e:
            return HttpResponseServerError(f"Wasn't able to parse the alignment file with error: {e.args[0]}")
        except:
            return HttpResponseServerError("Wasn't able to parse the alignment file! Is the file in FASTA format?")
        if len(alignments) == 0:
            return HttpResponseServerError("Wasn't able to parse the alignment file! Is the file in FASTA format?")
        if len(alignments) > 1:
            return HttpResponseServerError("Alignment file had more than one alignments!\nPlease upload a single alignment.")
        fastastring = format(alignments[0], "fasta")
        if validate_fasta_string(fastastring):
            cdHitTruncatedAln, cdHitReport = handleCDhit(alignments[0])
            request.session['custom_alignment_file'] = cdHitTruncatedAln
            request.session['cdHitUnTruncatedAln'] = fastastring
            request.session['cdHitReport'] = cdHitReport
            if cdHitReport:
                return HttpResponse('Success!')
            else:
                return HttpResponse('No CDHITS')
        else:
            return HttpResponseServerError("Alignment file had forbidden characters!\nWhat are you trying to do?")
    if request.method == 'GET':
        fastastring = request.session.get('custom_alignment_file')
        response_dict = handleCustomAlnGETRequest(fastastring)
        response_dict['cdHitReport'] = request.session['cdHitReport']
        return JsonResponse(response_dict, safe = False)

def getUntruncAln(request):
    fastastring = request.session.get('cdHitUnTruncatedAln')
    response_dict = handleCustomAlnGETRequest(fastastring)
    response_dict['cdHitReport'] = request.session['cdHitReport']
    return JsonResponse(response_dict, safe = False)

def handleCustomAlnGETRequest(fastastring):
    from alignments.Shannon import gap_adjusted_frequency
    fastastring = fastastring.replace('\n','\\n')
    concat_fasta, twc, gap_only_cols, filtered_spec_list, alignment_obj = calculateFastaProps(fastastring)
    frequency_list = list()
    for i in range(0, alignment_obj.get_alignment_length()):
        frequency_list.append(gap_adjusted_frequency(alignment_obj[:,i], IUPACData.protein_letters))
    response_dict = construct_dict_for_json_response([concat_fasta,filtered_spec_list,gap_only_cols,frequency_list,twc])

    return response_dict

def prepareCDHit (alnObj):
    seqNameDict = dict()
    cleanSeqs = MultipleSeqAlignment([])
    for i, record in enumerate(alnObj):
        seqNameDict[str(i)] = record.id
        cleanSeqs.add_sequence(str(i), str(record.seq).replace('-',''))
    return format(cleanSeqs, "fasta"), seqNameDict

def executeCDHit(fasta):
    from subprocess import Popen, PIPE
    from os import remove, path
    import datetime

    now = datetime.datetime.now()
    fileNameSuffix = "_" + str(now.year) + "_" + str(now.month) + "_" + str(now.day) + "_" + str(now.hour) + "_" + str(now.minute) + "_" + str(now.second) + "_" + str(now.microsecond)
    fastaName = f"/home/Desire-DEV/PVDev/static/cleanFastaCD{fileNameSuffix}.fa"
    cdHitOut = f"/home/Desire-DEV/PVDev/static/cleanFastaCD{fileNameSuffix}"
    cdHitClusters = f"/home/Desire-DEV/PVDev/static/cleanFastaCD{fileNameSuffix}.clstr"
    tempfiles = [fastaName, cdHitOut, cdHitClusters]
    for tempf in tempfiles:
        if path.isfile(tempf):
            remove(tempf)

    fh = open(fastaName, "w")
    fh.write(fasta)
    fh.close()

    pipe = Popen(f"cdhit -i {fastaName} -o {cdHitOut}; cat {cdHitClusters}", stdout=PIPE, shell=True)
    output = pipe.communicate()[0]

    if len(output.decode("ascii")) <= 0:
        for removeFile in tempfiles:
            remove(removeFile)
        return HttpResponseServerError("CDHIT failed!\nLikely a problem with the sequence alignment.")

    cdHitClusterOut = output.decode("ascii")
    
    for tempf in tempfiles:
        remove(tempf)
    
    clusters = parseCDHitClusters(cdHitClusterOut)
    
    return clusters, cdHitClusterOut

def parseCDHitClusters(cdHitclusterOut):
    cdHitListResults = cdHitclusterOut.split('\n>Cluster')
    resultNums = cdHitListResults[0].split('comparing sequences from')[1].split('\n')[2].split()
    if resultNums[0] == resultNums[2]:
        return False
    cdHitclusterStrings = cdHitListResults[1:]
    clusterSeqs = list()
    for clusterString in cdHitclusterStrings:
        clusterEntry = clusterString.split('>')
        for seqName in clusterEntry:
            entryNames = seqName.split('... *')
            if len(entryNames) == 2:
                clusterSeqs.append(entryNames[0])
                break
    return clusterSeqs

def truncateAlnByCDHitClusters(alnObj, clusters, seqNameDict):
    truncSeqRecordList, truncIdList = list(), list()
    for cluster in clusters:
        truncIdList.append(seqNameDict[cluster])
    for record in alnObj:
        if record.id in truncIdList:
            truncSeqRecordList.append(record)
    truncAln = MultipleSeqAlignment(truncSeqRecordList)
    return format(truncAln, "fasta")

def handleCDhit(alnObj):
    import re
    cleanFasta, seqNameDict = prepareCDHit (alnObj)
    cdHitclusters, cdHitReport = executeCDHit(cleanFasta)
    if not cdHitclusters:
        return format(alnObj, "fasta"), False
    translDict = dict()
    for num, name in seqNameDict.items():
        translDict[f'>{num}...'] = f'>{name}...'
    translCDHitReport = re.sub('({})'.format('|'.join(map(re.escape, translDict.keys()))), lambda m: translDict[m.group()], cdHitReport)
    newAlnString = truncateAlnByCDHitClusters(alnObj, cdHitclusters, seqNameDict)
    return newAlnString, translCDHitReport
