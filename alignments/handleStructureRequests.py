import io, json, os
from django.http import JsonResponse, HttpResponse, HttpResponseServerError
from Bio.PDB import PDBParser, MMCIFParser, PDBIO
from Bio.PDB.mmcifio import MMCIFIO
from pdbecif.mmcif_io import CifFileReader
import datetime
from alignments.views import parse_string_structure
from alignments.mapStrucSeqToAln import constructStrucSeqMap
import alignments.config

def handleCustomUploadStructure (request, strucID):
    '''We will POST all structure chains we need with uniqueIDs.
    Then when we GET them we can list the strucIDs separated by coma 
    this would mean "combine these in one CIF and return them".
    '''
    strucID = strucID
    if request.method == 'POST':
        from urllib.request import urlopen
        try:
            entities = request.POST.get("entities")
        except:
            return HttpResponseServerError("POST was sent without entities to parse!")
        deStrEnt = json.loads(entities)
        print('json')
        if strucID == "cust":
            #### This is not dependent on topology and should return success
            strucObj, cif_file_path = parseCustomCIF(deStrEnt["stringData"], "CUST")

            strucString = strucToString(strucObj)
           
            
            fixedEntityStruc = fixEntityFieldofParsedCIF(strucString, {deStrEnt["chainID"]:deStrEnt["entityID"]})
            outStruc = fixResiFieldsofParsedCIF(fixedEntityStruc)
            request.session[f'{strucID}-{deStrEnt["entityID"]}-{deStrEnt["chainID"]}'] = outStruc
            request.session[f'PDB-{strucID}-{deStrEnt["entityID"]}-{deStrEnt["chainID"]}'] = deStrEnt["stringData"]
    
            return JsonResponse({
                "successFlag" : True,
                "cif_file_path" : cif_file_path
            }, safe=False)
        for entry in deStrEnt:
            if request.session.get(f'{strucID}-{entry["entityID"]}-{entry["chainID"]}'):
                continue
            ebiURL = f'https://coords.litemol.org/{strucID.lower()}/chains?entityId={entry["entityID"]}&authAsymId={entry["chainID"]}&atomSitesOnly=1'
            #ebiURL = f'https://www.ebi.ac.uk/pdbe/coordinates/{strucID.lower()}/chains?entityId={entityId}&atomSitesOnly=1'
            try:
                data = urlopen(ebiURL)
            except:
                return HttpResponseServerError(f'Failed to fetch coordinates from litemol for PDB {strucID} and entityID {entry["entityID"]} and chain id {entry["chainID"]}.')
            try:
                tempStrucStr = str()
                for line in data:
                    tempStrucStr+=line.decode('UTF-8')
            except:
                return HttpResponseServerError("Failed to parse the provided structure!")
            request.session[f'{strucID}-{entry["entityID"]}-{entry["chainID"]}'] = tempStrucStr
        return JsonResponse("Success!", safe=False)
    
    if request.method == 'GET':
        if request.session.get(strucID):
            stringStruc = request.session[strucID]
            return HttpResponse(stringStruc, content_type="text/plain")
        structureList = list()
        chainToEntity = dict()
        for singleID in strucID.split(','):
            pdbAndEntityAndChain = singleID.split('-')
            try:
                stringData = request.session[singleID]
            except:
                return HttpResponseServerError(f"Requested structure {singleID} was not present in the session! Wait for POST to finish.")
            strucObj = parse_string_structure(request, stringData, pdbAndEntityAndChain[0])
            chainToEntity[pdbAndEntityAndChain[2]] = pdbAndEntityAndChain[1]
            structureList.append(strucObj)
        try:
            if len(structureList) > 1:
                structureList = combineChainsInSingleStruc(structureList)
            stringStruc = strucToString(structureList[0])
            outStruc = fixEntityFieldofParsedCIF(stringStruc, chainToEntity)
            request.session[strucID] = outStruc
            stringStruc = outStruc
        except:
            return HttpResponseServerError("Failed to parse structures!")
        return HttpResponse(stringStruc, content_type="text/plain")

def handleCustomUploadStructure_CIF (request, strucID):
    '''We will POST all structure chains we need with uniqueIDs.
    Then when we GET them we can list the strucIDs separated by coma 
    this would mean "combine these in one CIF and return them".
    '''
    strucID = strucID
    if request.method == 'POST':
        from urllib.request import urlopen
        try:
            entities = request.POST.get("entities")
        except:
            return HttpResponseServerError("POST was sent without entities to parse!")
        deStrEnt = json.loads(entities)
        print('json')
        if strucID == "cust":
            #### This is not dependent on topology and should return success
            strucObj, cif_file_path = parseCustomCIF(deStrEnt["stringData"], "CUST")

            strucString = strucToString(strucObj)
           
            print('structString')
            fixedEntityStruc = fixEntityFieldofParsedCIF(strucString, {deStrEnt["chainID"]:deStrEnt["entityID"]})
            outStruc_full = fixResiFieldsofParsedCIF(fixedEntityStruc)
            outStruc = delResiofParsedPDB(outStruc_full)
            request.session[f'{strucID}-{deStrEnt["entityID"]}-{deStrEnt["chainID"]}'] = outStruc
            request.session[f'PDB-{strucID}-{deStrEnt["entityID"]}-{deStrEnt["chainID"]}'] = deStrEnt["stringData"]
    
    
            return JsonResponse({
                "successFlag" : True,
                "cif_file_path" : cif_file_path
            }, safe=False)
        for entry in deStrEnt:
            if request.session.get(f'{strucID}-{entry["entityID"]}-{entry["chainID"]}'):
                continue
            ebiURL = f'https://coords.litemol.org/{strucID.lower()}/chains?entityId={entry["entityID"]}&authAsymId={entry["chainID"]}&atomSitesOnly=1'
            #ebiURL = f'https://www.ebi.ac.uk/pdbe/coordinates/{strucID.lower()}/chains?entityId={entityId}&atomSitesOnly=1'
            try:
                data = urlopen(ebiURL)
            except:
                return HttpResponseServerError(f'Failed to fetch coordinates from litemol for PDB {strucID} and entityID {entry["entityID"]} and chain id {entry["chainID"]}.')
            try:
                tempStrucStr = str()
                for line in data:
                    tempStrucStr+=line.decode('UTF-8')
            except:
                return HttpResponseServerError("Failed to parse the provided structure!")
            request.session[f'{strucID}-{entry["entityID"]}-{entry["chainID"]}'] = tempStrucStr
        return JsonResponse("Success!", safe=False)
    
    if request.method == 'GET':
        if request.session.get(strucID):
            stringStruc = request.session[strucID]
            return HttpResponse(stringStruc, content_type="text/plain")
        structureList = list()
        chainToEntity = dict()
        for singleID in strucID.split(','):
            pdbAndEntityAndChain = singleID.split('-')
            try:
                stringData = request.session[singleID]
            except:
                return HttpResponseServerError(f"Requested structure {singleID} was not present in the session! Wait for POST to finish.")
            strucObj = parse_string_structure(request, stringData, pdbAndEntityAndChain[0])
            chainToEntity[pdbAndEntityAndChain[2]] = pdbAndEntityAndChain[1]
            structureList.append(strucObj)
        try:
            if len(structureList) > 1:
                structureList = combineChainsInSingleStruc(structureList)
            stringStruc = strucToString(structureList[0])
            outStruc = fixEntityFieldofParsedCIF(stringStruc, chainToEntity)
            request.session[strucID] = outStruc
            stringStruc = outStruc
        except:
            return HttpResponseServerError("Failed to parse structures!")
        return HttpResponse(stringStruc, content_type="text/plain")

def handleCustomUploadStructure_PDB (request, strucID):
    '''We will POST all structure chains we need with uniqueIDs.
    Then when we GET them we can list the strucIDs separated by coma 
    this would mean "combine these in one CIF and return them".
    '''
    strucID = strucID
    if request.method == 'POST':
        from urllib.request import urlopen
        try:
            entities = request.POST.get("entities")
        except:
            return HttpResponseServerError("POST was sent without entities to parse!")
        deStrEnt = json.loads(entities)
        print('json')
        if strucID == "cust":
            #### This is not dependent on topology and should return success
            strucObj = parseCustomPDB(deStrEnt["stringData"], "CUST")

            strucString = strucToString(strucObj)
           
            
            fixedEntityStruc = fixEntityFieldofParsedPDB(strucString, {deStrEnt["chainID"]:deStrEnt["entityID"]})
            outStruc_full = fixResiFieldsofParsedPDB(fixedEntityStruc)
            outStruc = delResiofParsedPDB(outStruc_full)
            request.session[f'{strucID}-{deStrEnt["entityID"]}-{deStrEnt["chainID"]}'] = outStruc
            request.session[f'PDB-{strucID}-{deStrEnt["entityID"]}-{deStrEnt["chainID"]}'] = deStrEnt["stringData"]
    
            return JsonResponse("Success!", safe=False)
        for entry in deStrEnt:
            if request.session.get(f'{strucID}-{entry["entityID"]}-{entry["chainID"]}'):
                continue
            ebiURL = f'https://coords.litemol.org/{strucID.lower()}/chains?entityId={entry["entityID"]}&authAsymId={entry["chainID"]}&atomSitesOnly=1'
            #ebiURL = f'https://www.ebi.ac.uk/pdbe/coordinates/{strucID.lower()}/chains?entityId={entityId}&atomSitesOnly=1'
            try:
                data = urlopen(ebiURL)
            except:
                return HttpResponseServerError(f'Failed to fetch coordinates from litemol for PDB {strucID} and entityID {entry["entityID"]} and chain id {entry["chainID"]}.')
            try:
                tempStrucStr = str()
                for line in data:
                    tempStrucStr+=line.decode('UTF-8')
            except:
                return HttpResponseServerError("Failed to parse the provided structure!")
            request.session[f'{strucID}-{entry["entityID"]}-{entry["chainID"]}'] = tempStrucStr
        return JsonResponse("Success!", safe=False)
    
    if request.method == 'GET':
        if request.session.get(strucID):
            stringStruc = request.session[strucID]
            return HttpResponse(stringStruc, content_type="text/plain")
        structureList = list()
        chainToEntity = dict()
        for singleID in strucID.split(','):
            pdbAndEntityAndChain = singleID.split('-')
            try:
                stringData = request.session[singleID]
            except:
                return HttpResponseServerError(f"Requested structure {singleID} was not present in the session! Wait for POST to finish.")
            strucObj = parse_string_structure(request, stringData, pdbAndEntityAndChain[0])
            chainToEntity[pdbAndEntityAndChain[2]] = pdbAndEntityAndChain[1]
            structureList.append(strucObj)
        try:
            if len(structureList) > 1:
                structureList = combineChainsInSingleStruc(structureList)
            stringStruc = strucToString(structureList[0])
            outStruc = fixEntityFieldofParsedCIF(stringStruc, chainToEntity)
            request.session[strucID] = outStruc
            stringStruc = outStruc
        except:
            return HttpResponseServerError("Failed to parse structures!")
        return HttpResponse(stringStruc, content_type="text/plain")

def fixEntityFieldofParsedCIF(stringStruc, chainToEntity):
    listStruc = stringStruc.split('\n')
    outStruc = listStruc[:21]
    for row in listStruc[21:]:
        rowList = row.split()
        if len(rowList) > 10:
           
            if rowList[16] == str(list(chainToEntity.keys())[0]):
                
                rowList[7] = chainToEntity[rowList[16]]
        outStruc.append(' '.join(rowList))
    return '\n'.join(outStruc)

def fixEntityFieldofParsedPDB(stringStruc, chainToEntity):
    listStruc = stringStruc.split('\n')
    outStruc = listStruc[:21]
    for row in listStruc[21:]:
        rowList = row.split()
        if len(rowList) > 10:    
            rowList[7] = chainToEntity[rowList[16]]         
        outStruc.append(' '.join(rowList))
    return '\n'.join(outStruc) 
   
def delResiofParsedPDB(stringStruc):
    listStruc = stringStruc.split('\n')
    outStruc = listStruc[:21]
    for row in listStruc[21:]:
        rowList = row.split()
        if len(rowList) >6:
            if str(rowList[5])=='URA':
               rowList[5]=='U'
            if str(rowList[5])=='GUA':
               rowList[5]=='G'
            if str(rowList[5])=='CYT':
               rowList[5]=='C'
            if str(rowList[5])=='ADE':
               rowList[5]=='A'          
            if str(rowList[5])=='U' or str(rowList[5])=='A' or str(rowList[5])=='G' or str(rowList[5])=='C':
                outStruc.append(' '.join(rowList))
    return '\n'.join(outStruc)   

def delResiofParsedCIF(stringStruc):
    listStruc = stringStruc.split('\n')
    outStruc = listStruc[:21]
    for row in listStruc[21:]:
        rowList = row.split()
        if len(rowList) >6:
            if str(rowList[5])=='URA':
               rowList[5]=='U'
            if str(rowList[5])=='GUA':
               rowList[5]=='G'
            if str(rowList[5])=='CYT':
               rowList[5]=='C'
            if str(rowList[5])=='ADE':
               rowList[5]=='A'          
            if str(rowList[5])=='U' or str(rowList[5])=='A' or str(rowList[5])=='G' or str(rowList[5])=='C':
                outStruc.append(' '.join(rowList))
    return '\n'.join(outStruc)   

def fixResiFieldsofParsedCIF(stringStruc):
    listStruc = stringStruc.split('\n')
    outStruc = listStruc[:21]
    
    for row in listStruc[21:]:
        rowList = row.split()
        
        if len(rowList) > 10:
            rowList[8] = rowList[15]
        outStruc.append(' '.join(rowList))
    return '\n'.join(outStruc)

def fixResiFieldsofParsedPDB(stringStruc):
    listStruc = stringStruc.split('\n')
    outStruc = listStruc[:21]
    #print(outStruc)
    for row in listStruc[21:]:
        rowList = row.split()
        if len(rowList) > 10:
            rowList[8] = rowList[15]
        
        
        outStruc.append(' '.join(rowList))
    return '\n'.join(outStruc)

def strucToString(strucObj):
    
    strucFile = io.StringIO("")
    mmCIFio=MMCIFIO()
    mmCIFio.set_structure(strucObj)
    mmCIFio.save(strucFile)    
    return strucFile.getvalue()

def strucToPDBString(strucObj):
    strucFile = io.StringIO("")
    pdbIO=PDBIO()
    chainmap = rename_chains(strucObj)
    pdbIO.set_structure(strucObj)
    pdbIO.save(strucFile)
    return strucFile.getvalue()

def combineChainsInSingleStruc(structureList):
    for singleStruc in structureList:
        chains = list(singleStruc.get_chains())
        for removeChain in chains[1:]:
            singleStruc[0].detach_child(removeChain.id)
    for strucToMerge in structureList[1:]:
        chain = list(strucToMerge.get_chains())
        structureList[0][0].add(chain[0])
    return structureList

def parseCustomPDB(stringData, id = "CUST"):
    parser = PDBParser()
    strucFile = io.StringIO(stringData)
    structureObj = parser.get_structure(id, strucFile)
    return structureObj

def parseCustomCIF(stringData, pdbid):
    #import config
    parser = MMCIFParser()
    strucFile = io.StringIO(stringData)
    strucFile1 = io.StringIO(stringData)
    ###NEED time here
    now = datetime.datetime.now()
    cif_fileNameSuffix = "_" + str(now.year) + "_" + str(now.month) + "_" + str(now.day) + "_" + str(now.hour) + "_" + str(now.minute) + "_" + str(now.second) + "_" + str(now.microsecond)
    alignments.config.cif_fileNameSuffix_share=cif_fileNameSuffix
    cif_file_path = f'/tmp/cust2{cif_fileNameSuffix}.cif'
    with open(cif_file_path, 'w') as f:
            f.write(strucFile1.read())
            f.close()
    structureObj = parser.get_structure(pdbid,strucFile)
    print('MMCIF parsed')
   
    print(cif_file_path)
    return structureObj, cif_file_path

def postTopology(request, strucID):
    if request.method == 'POST':
        try:
            topologyJSON=[]
        except:
            return HttpResponseServerError("Failed to generate topology from the provided structure!")
        request.session[f'TOPOLOGY-{strucID}'] = topologyJSON
        return JsonResponse("Success!", safe=False)
    return JsonResponse("Only for custom structure post!", safe=False)

def getTopology (request, topID):
    if topID == "EMPTY":
        return JsonResponse({}, safe=False)
    if request.session.get(topID):
        topology = request.session[topID]
        return JsonResponse(topology, safe=False)


class OutOfChainsError(Exception): pass
def rename_chains(structure):
    """Renames chains to be one-letter chains
    
    Existing one-letter chains will be kept. Multi-letter chains will be truncated
    or renamed to the next available letter of the alphabet.
    
    If more than 62 chains are present in the structure, raises an OutOfChainsError
    
    Returns a map between new and old chain IDs, as well as modifying the input structure
    """
    next_chain = 0 #
    # single-letters stay the same
    chainmap = {c.id:c.id for c in structure.get_chains() if len(c.id) == 1}
    for o in structure.get_chains():
        if len(o.id) != 1:
            if o.id[0] not in chainmap:
                chainmap[o.id[0]] = o.id
                o.id = o.id[0]
            else:
                c = int_to_chain(next_chain)
                while c in chainmap:
                    next_chain += 1
                    c = int_to_chain(next_chain)
                    if next_chain >= 62:
                        raise OutOfChainsError()
                chainmap[c] = o.id
                o.id = c
    return chainmap

def int_to_chain(i,base=62):
    """
    int_to_chain(int,int) -> str
    Converts a positive integer to a chain ID. Chain IDs include uppercase
    characters, numbers, and optionally lowercase letters.
    i = a positive integer to convert
    base = the alphabet size to include. Typically 36 or 62.
    """
    if i < 0:
        raise ValueError("positive integers only")
    if base < 0 or 62 < base:
        raise ValueError("Invalid base")

    quot = int(i)//base
    rem = i%base
    if rem < 26:
        letter = chr( ord("A") + rem)
    elif rem < 36:
        letter = str( rem-26)
    else:
        letter = chr( ord("a") + rem - 36)
    if quot == 0:
        return letter
    else:
        return int_to_chain(quot-1,base) + letter

