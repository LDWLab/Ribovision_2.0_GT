import io
from django.http import JsonResponse, HttpResponse, HttpResponseServerError
from Bio.PDB import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO

def handleCustomUploadStructure (request, strucID):
    '''We will POST all structure chains we need with uniqueIDs.
    Then when we GET them we can list the strucIDs separated by coma 
    this would mean "combine these in one CIF and return them".
    '''
    strucID = strucID.upper()
    if request.method == 'POST':
        from urllib.request import urlopen
        try:
            entityIDS = request.POST.getlist("entityIDS")
        except:
            return HttpResponseServerError("POST was sent without entities to parse!")
        for entityId in entityIDS:
            if request.session.get(f'{strucID}-{entityId}'):
                continue
            ebiURL = f'https://www.ebi.ac.uk/pdbe/coordinates/{strucID.lower()}/chains?entityId={entityId}&atomSitesOnly=1'
            try:
                data = urlopen(ebiURL)
            except:
                return HttpResponseServerError(f'Failed to fetch coordinates from PDBe for PDB {strucID} and entityID {entityId}')
            try:
                tempStrucStr = str()
                for line in data:
                    tempStrucStr+=line.decode('UTF-8')
            except:
                return HttpResponseServerError("Failed to parse the provided structure!")
            request.session[f'{strucID}-{entityId}'] = tempStrucStr
        return JsonResponse("Success!", safe=False)
    
    if request.method == 'GET':
        if request.session.get(strucID):
            stringStruc = request.session[strucID]
            return HttpResponse(stringStruc, content_type="text/plain")
        structureList = list()
        #try:
        chainToEntity = dict()
        for singleID in strucID.split(','):
            pdbAndEntity = singleID.split('-')
            stringData = request.session[singleID]
            strucObj = parse_string_structure(stringData, pdbAndEntity[0])
            chainToEntity[list(strucObj[0].get_chains())[0].id] = pdbAndEntity[1]
            structureList.append(strucObj)
        if len(structureList) > 1:
            for singleStruc in structureList:
                chains = list(singleStruc.get_chains())
                for removeChain in chains[1:]:
                    singleStruc[0].detach_child(removeChain.id)
            for strucToMerge in structureList[1:]:
                chain = list(strucToMerge.get_chains())
                structureList[0][0].add(chain[0])
        stringStruc = strucToString(structureList[0])
        listStruc = stringStruc.split('\n')
        tempStruc = listStruc[:21]
        for row in listStruc[21:]:
            rowList = row.split()
            if len(rowList) > 10:
                rowList[7] = chainToEntity[rowList[16]]
            tempStruc.append(' '.join(rowList))
        request.session[strucID] = '\n'.join(tempStruc)
        stringStruc = '\n'.join(tempStruc)
        #except:
        #    return HttpResponseServerError("Failed to parse structures!")
        return HttpResponse(stringStruc, content_type="text/plain")

def parse_string_structure(stringData, strucID):

    parser = MMCIFParser()
    strucFile = io.StringIO(stringData)
    structureObj = parser.get_structure(strucID,strucFile)
    return structureObj

def strucToString(strucObj):
    
    strucFile = io.StringIO("")
    mmCIFio=MMCIFIO()
    mmCIFio.set_structure(strucObj)
    mmCIFio.save(strucFile)
    return strucFile.getvalue()