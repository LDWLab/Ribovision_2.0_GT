from django.shortcuts import render
from django.db import connection
from ribovision.models import Mastertable, Secondarystructures, Species
from django.http import HttpResponse
from django.http import JsonResponse
from django.views.generic import ListView
from django.views.decorators.csrf import csrf_exempt

def fetchmasterlist(request):
    response = Mastertable.objects.values('SpeciesName', 'DataSetType', 'StructureName', 'LoadString').filter(Active = True)
    return JsonResponse(list(response), safe = False)

def get_secstr_data(secstr_name):
    secstr_data = Secondarystructures.objects.filter(name = secstr_name)
    return secstr_data

def get_species_data(species_abbreviation):
    specie_data = Species.objects.filter(abbreviation = species_abbreviation)
    return specie_data

def fetchstructurename(request):
    '''Requires:
        StructureName: "4V9D"
    '''
    if request.method == "POST":
        #response = Threedstructures.objects.values('structurename').filter(number_3d_structure_id = )
        pass

def speciestable(request):
    '''Requires :
        Circle_Radius: 1.7
        Font_Size_Canvas: 3.1
        Font_Size_SVG: 3.9
        SS_Table: "ECOLI_SSU"
        Species_Abr: "ECOLI"
        Species_Name: "Escherichia coli"
    '''
    if request.method == "POST":
        pass

def fetchresidues(request):
    '''Requires for each entry:
        ChainName: "CA"
        Domain_AN: 1
        Domain_RN: "I"
        Domains_Color: 1
        Helix_Color: 1
        Helix_Num: "H1"
        MoleculeGroup: "LSU"
        MoleculeType: "rRNA"
        SS_Table: "ECOLI_LSU"
        X: "405.009"
        Y: "484.355"
        map_Index: 1
        modResName: "G"
        molName: "23S"
        resNum: "1"
        unModResName: "G"
    '''
    if request.method == "POST":
        structure_identity = request.body
        SQLStatement = 'SELECT SS_Data.map_index, GeneSymbol, resNum, X, Y, unModResName, modResName, polymer_type, \
                        MoleculeGroup, ChainName, Domain_RN, Domain_AN, Domains_Color, Helix_Num, Helix_Color\
                        FROM (SELECT * FROM SecondaryStructures WHERE Name = '+structure_identity.decode()+') AS ss\
                        LEFT JOIN SS_Data ON ss.SecStr_id = SS_Data.ss_id\
                        LEFT JOIN Residues ON SS_Data.res_id = Residues.resi_id\
                        LEFT JOIN Polymer_Data ON Residues.PolData_id = Polymer_Data.PData_id\
                        LEFT JOIN Polymer_metadata ON Residues.PolData_id = Polymer_metadata.polymer_id\
                        LEFT JOIN ChainList ON Polymer_Data.PData_id = ChainList.polymer_id\
                        INNER JOIN StructuralData2 ON ss.SecStr_id = StructuralData2.secondary_structure_id\
                        WHERE SS_Data.map_index = StructuralData2.map_index'
        with connection.cursor() as cursor:
            cursor.execute(SQLStatement)
            response = cursor.fetchall()
        
    return JsonResponse(list(response), safe = False)

def textlabels(request):
    '''Requires for each entry:
        Fill: "#FCC595"
        Font: "MyriadPro-Regular"
        FontSize: 24
        LabelText: "0"
        SS_Table: "ECOLI_LSU"
        X: 365.759
        Y: 332.807
        id: 633
    '''
    if request.method == "POST":
        pass

def linelabels(request):
    '''Requires for each entry:
        Fill: ""
        SS_Table: "ECOLI_LSU"
        Stroke: "#231F20"
        StrokeLineJoin: "round"
        StrokeMiterLimit: 10
        StrokeWidth: 0.25
        X1: 269.233
        X2: 264.649
        Y1: 302.506
        Y2: 297.84
        id: 623
    '''
    if request.method == "POST":
        pass

def fetchinteractionsmenu(request):
    '''Requires:
        [{bp_group: "BasePairs"}, {bp_group: "Stacking"}, {bp_group: "BaseSugar"}, {bp_group: "BasePhosphate"}]
    '''
    if request.method == "POST":
        pass

#Can skip that one for now
def structdatamenu(request):
    '''Requires for each entry:
        ColName: "mean_tempFactor"
        ColorList: "Rainbow1"
        Description: "Mean B-factor per nucleotide from the crystal structure."
        ExtraArg: ""
        HelpLink: "StructuralData"
        IndexMode: "FALSE"
        StructDataName: "B factors"
        StructureName: "4V9D"
    '''
    if request.method == "POST":
        pass

def index(request):
    return render(request, 'ribovision/index.html')