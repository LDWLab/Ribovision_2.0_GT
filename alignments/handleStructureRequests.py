import io
from django.http import JsonResponse, HttpResponse, HttpResponseServerError

def handleCustomUploadStructure (request, strucID):
    '''We will POST all structure chains we need with uniqueIDs.
    Then when we GET them we can list the strucIDs separated by coma 
    this would mean "combine these in one CIF and return them".
    '''
    if request.method == 'POST':
        import json
        from urllib.request import urlopen
        try:
            entityIDS = request.POST.getlist("entityIDS")
        except:
            return HttpResponseServerError("POST was sent without entities to parse!")
        for entityId in entityIDS:
            ebiURL = f'https://www.ebi.ac.uk/pdbe/coordinates/{strucID.lower()}/chains?entityId={entityId}'
            try:
                data = urlopen(ebiURL)
            except:
                return HttpResponseServerError(f'Failed to fetch coordinates from PDBe for PDB {strucID} and entityID {entityId}')
            try:
                tempStrucList = str()
                for line in data:
                    tempStrucList+=line.decode('UTF-8')
                serializeData = json.dumps(tempStrucList)
            except:
                return HttpResponseServerError("Failed to parse the provided structure!")
            request.session[f'{strucID}_{entityId}'] = serializeData
        return JsonResponse("Success!", safe=False)
    
    if request.method == 'GET':
        structureDict = dict()
        for singleID in strucID.split(','):
            serializeData = request.session[singleID]
            strucObj = parse_serialized_structure(serializeData, singleID)
            structureDict[singleID] = strucObj
        if len(structureDict) == 1:
            stringStruc = strucToString(next(iter(structureDict.values())))
            return HttpResponse(stringStruc, content_type="text/plain")
        elif len(structureDict) > 1:
            #combine into single structure and return
            pass
        else:
            return HttpResponseServerError("Failed to parse structures!")

def parse_serialized_structure(serializeData, strucID):

    from Bio.PDB import MMCIFParser
    parser = MMCIFParser()
    strucFile = io.StringIO(serializeData.replace('\\n','\n')[1:-1])
    structureObj = parser.get_structure(strucID,strucFile)
    return structureObj

def strucToString(strucObj):
    
    from Bio.PDB.mmcifio import MMCIFIO
    strucFile = io.StringIO("")
    mmCIFio=MMCIFIO()
    mmCIFio.set_structure(strucObj)
    mmCIFio.save(strucFile)
    return strucFile.getvalue()