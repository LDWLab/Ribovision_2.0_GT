from django.http import JsonResponse, HttpResponseServerError
from Bio.Align import MultipleSeqAlignment

def prepareCDHit (alnObj):
    seqNameDict = dict()
    cleanSeqs = MultipleSeqAlignment([])
    for i, record in enumerate(alnObj):
        seqNameDict[str(i)] = record.id
        cleanSeqs.add_sequence(str(i), str(record.seq).replace('-',''))
    return format(cleanSeqs, "fasta")

def executeCDHit():
    pass
    #save fasta file
    #execute cdhit
    #read cdhit results


def parseCDHit():
    pass
    #parse cdhit results
    #cleanup results and aln file
    #store as session var
