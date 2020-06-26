from django.shortcuts import render
from ribovision.models import Mastertable
from django.http import HttpResponse
from django.http import JsonResponse
from django.views.generic import ListView
from django.views.decorators.csrf import csrf_exempt

def fetchmasterlist(request):
    response = Mastertable.objects.values('SpeciesName', 'DataSetType', 'StructureName', 'LoadString').filter(Active = True)
    return JsonResponse(list(response), safe = False)

def fetchresidues(request):
    if request.method == "POST":
        structure_identity = request.body
        SQLStatement = 'SELECT ss.map_Index, ss.molName, resNum, X, Y, unModResName, modResName, MoleculeType, \
                    MoleculeGroup, ChainName, Domain_RN, Domain_AN, Domains_Color, Helix_Num, Helix_Color, ss.SS_Table \
                    FROM (SELECT * FROM SecondaryStructures WHERE SS_Table = '+structure_identity.decode()+') AS ss \
                    INNER JOIN MoleculeNames AS mn ON ss.molName = mn.MoleculeName \
                    INNER JOIN (SELECT MoleculeName, ChainName FROM ChainList \
                    WHERE StructureName = (SELECT DISTINCT StructureName FROM Secondary_Tertiary WHERE SS_Table = '+structure_identity.decode()+')) \
                    AS cl ON cl.MoleculeName = ss.molName \
                    LEFT JOIN (SELECT * FROM StructuralData2 WHERE SS_Table = '+structure_identity.decode()+') AS sd ON sd.map_Index = ss.map_Index'
    
    return JsonResponse(list(SQLStatement), safe = False)

def index(request):
    return render(request, 'ribovision/index.html')