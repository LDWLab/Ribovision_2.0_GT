from django.http import JsonResponse, HttpResponseServerError
from Bio.Align import MultipleSeqAlignment

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
    fastaName = f"./static/cleanFastaCD{fileNameSuffix}.fa"
    cdHitOut = f"./static/cleanFastaCD{fileNameSuffix}"
    cdHitClusters = f"./static/cleanFastaCD{fileNameSuffix}.clstr"
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
    cdHitclusterStrings = cdHitclusterOut.split('\n>Cluster')[1:]
    clusterSeqs = list()
    for clusterString in cdHitclusterStrings:
        clusterEntry = clusterString.split('>')
        clusterSeqs.append(clusterEntry[1].split('... *')[0])
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
    translDict = dict()
    for num, name in seqNameDict.items():
        translDict[f'>{num}...'] = f'>{name}...'
    translCDHitReport = re.sub('({})'.format('|'.join(map(re.escape, translDict.keys()))), lambda m: translDict[m.group()], cdHitReport)
    newAlnString = truncateAlnByCDHitClusters(alnObj, cdHitclusters, seqNameDict)
    return newAlnString, translCDHitReport
