from django.shortcuts import render
from django.db import connection
from ribovision.models import Mastertable, Secondarystructures, Species
from django.http import HttpResponse
from django.http import JsonResponse
from django.views.generic import ListView
from django.views.decorators.csrf import csrf_exempt
import re

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
        structure_identity = request.body.decode()
        SQLStatement = 'SELECT DISTINCT StructureName FROM ThreeDStructures WHERE 3D_structure_id = \
            (SELECT 3D_structure_id FROM Secondary_Tertiary WHERE secondary_structure_id = (\
                SELECT SecStr_id FROM SecondaryStructures WHERE Name = '+structure_identity+'))'
        with connection.cursor() as cur:
            cur.execute(SQLStatement)
            response = [dict((cur.description[i][0], value)
            for i, value in enumerate(row)) for row in cur.fetchall()]
    return JsonResponse(list(response), safe = False)

def speciestable(request):
    '''Requires :
        Circle_Radius: 1.7
        Font_Size_Canvas: 3.1
        Font_Size_SVG: 3.9
        SS_Table: "ECOLI_SSU"
        Species_Abr: "ECOLI"
        Species_Name: "Escherichia coli"
        Molecule_Names: ??? "23S:5S"
    '''
    if request.method == "POST":
        structure_identity = request.body.decode()
        SQLStatement = 'SELECT CAST(Circle_Radius AS CHAR) AS Circle_Radius, CAST(Font_Size_Canvas AS CHAR) AS Font_Size_Canvas, CAST(Font_Size_SVG AS CHAR) AS Font_Size_SVG,\
            SecondaryStructures.Name as \'SS_Table\', Abbreviation as \'Species_Abr\', Species.name as \'Species_Name\' FROM SecondaryStructures\
            INNER JOIN Species ON SecondaryStructures.strain_fk = Species.strain_id\
            WHERE SecondaryStructures.Name = '+structure_identity
        #SQLStatement_with_polymers = 'SELECT CAST(Circle_Radius AS CHAR) AS Circle_Radius, CAST(Font_Size_Canvas AS CHAR) AS Font_Size_Canvas,\
        #    CAST(Font_Size_SVG AS CHAR) AS Font_Size_SVG, CAST(GeneSymbol AS CHAR) AS Molecule_Names, ss.Name as \'SS_Table\', \
        #    Abbreviation as \'Species_Abr\', Species.name as \'Species_Name\' FROM (SELECT * FROM DESIRE.SecondaryStructures WHERE Name = '+structure_identity+') as ss \
        #    INNER JOIN Species ON ss.strain_fk = Species.strain_id\
        #    INNER JOIN Secondary_Tertiary ON ss.SecStr_id = Secondary_Tertiary.secondary_structure_id \
        #    INNER JOIN ThreeDStructures ON Secondary_Tertiary.3D_structure_id = ThreeDStructures.3D_structure_id \
        #    INNER JOIN ChainList ON ThreeDStructures.3D_structure_id = ChainList.3D_structure_id \
        #    INNER JOIN Polymer_Data on ChainList.polymer_id = Polymer_Data.PData_id\
        #    INNER JOIN Nomenclature on Polymer_Data.nomgd_id = Nomenclature.nom_id\
        #    WHERE Nomenclature.MoleculeGroup='+structure_identity.split('_')[1]
        with connection.cursor() as cur:
            cur.execute(SQLStatement)
            response = [dict()]
            for row in cur.fetchall():
                for i, value in enumerate(row):
                    if value:       #Not good cus it overwrites the Molecule_Names; have to fix it later
                        if re.match(r"\d{1}\.{1}\d{1}", value):
                            response[0][cur.description[i][0]] = float(value)
                        else:
                            response[0][cur.description[i][0]] = value
    if bool(response[0]):
        return JsonResponse(list(response), safe = False)
    else:
        return JsonResponse(list(), safe = False)

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
        structure_identity = request.body.decode()
        SQLStatement = 'SELECT SS_Data.map_index, CAST(GeneSymbol AS CHAR) AS molName, resNum, X, Y, unModResName, modResName, \
                        CAST(polymer_type AS CHAR) AS MoleculeType, MoleculeGroup, ChainName, Domain_RN, Domain_AN, Domains_Color, Helix_Num, Helix_Color\
                        FROM (SELECT * FROM SecondaryStructures WHERE Name = '+structure_identity+') AS ss\
                        LEFT JOIN SS_Data ON ss.SecStr_id = SS_Data.ss_id\
                        LEFT JOIN Residues ON SS_Data.res_id = Residues.resi_id\
                        LEFT JOIN Polymer_Data ON Residues.PolData_id = Polymer_Data.PData_id\
                        LEFT JOIN Polymer_metadata ON Residues.PolData_id = Polymer_metadata.polymer_id\
                        LEFT JOIN ChainList ON Polymer_Data.PData_id = ChainList.polymer_id\
                        INNER JOIN StructuralData2 ON ss.SecStr_id = StructuralData2.secondary_structure_id\
                        WHERE SS_Data.map_index = StructuralData2.map_index'
        SQLStatement_old = 'SELECT ss.map_Index, ss.molName, resNum, X, Y, unModResName, modResName, MoleculeType, \
                MoleculeGroup, ChainName, Domain_RN, Domain_AN, Domains_Color, Helix_Num, Helix_Color, ss.SS_Table \
                FROM (SELECT * FROM ribovision2.SecondaryStructures WHERE SS_Table = '+structure_identity+') AS ss \
                INNER JOIN ribovision2.MoleculeNames AS mn ON ss.molName = mn.MoleculeName \
                INNER JOIN (SELECT MoleculeName, ChainName FROM ribovision2.ChainList \
                WHERE StructureName = (SELECT DISTINCT StructureName FROM ribovision2.Secondary_Tertiary WHERE SS_Table = '+structure_identity+')) \
                AS cl ON cl.MoleculeName = ss.molName \
                LEFT JOIN (SELECT * FROM ribovision2.StructuralData2 WHERE SS_Table = '+structure_identity+') AS sd ON sd.map_Index = ss.map_Index'

        with connection.cursor() as cursor:
            cursor.execute(SQLStatement)
            response = [dict((cursor.description[i][0], value)
            for i, value in enumerate(row)) for row in cursor.fetchall()]
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
        structure_identity = request.body.decode()
        SQLStatement = 'SELECT Fill, Font, Font_Size as \'FontSize\', LabelText, TextLabel_id as \'id\', X, Y\
        FROM TextLabels WHERE secondary_structure_id = (SELECT SecStr_id FROM SecondaryStructures WHERE Name = '+structure_identity+')'
        with connection.cursor() as cur:
            cur.execute(SQLStatement)
            response = [dict((cur.description[i][0], value)
            for i, value in enumerate(row)) for row in cur.fetchall()]
        return JsonResponse(list(response), safe = False)

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
        structure_identity = request.body.decode()
        SQLStatement = 'SELECT * FROM LineLabels WHERE secondary_structure_id = (SELECT SecStr_id FROM SecondaryStructures WHERE Name = '+structure_identity+')'
        with connection.cursor() as cur:
            cur.execute(SQLStatement)
            response = [dict((cur.description[i][0], value)
            for i, value in enumerate(row)) for row in cur.fetchall()]
        return JsonResponse(list(response), safe = False)

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