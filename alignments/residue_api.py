from django.http import JsonResponse, Http404
from alignments.models import *
from django.db import connection

def get_taxid_from_polid(polymer):
    try:
        taxid = PolymerData.objects.filter(pdata_id = polymer)[0].strain_id
    except PolymerData.DoesNotExist:
        raise Http404("No taxid with matching polymer id"+str(polymer)+" in the database!")
    return taxid

def get_alignments_from_resid(resid):
    alns = Alignment.objects.filter(alndata__res = resid)
    if len(alns) == 0:
        raise Http404("Provided residue id does not have an alignment!")
    return alns

def get_superkingdom_id(taxid):
    sql_query = 'SELECT taxgroup_id FROM DESIRE.TaxGroups\
                WHERE groupLevel = \'superkingdom\' AND\
                taxgroup_id IN (SELECT taxgroup_id FROM DESIRE.Species_TaxGroup \
                                WHERE strain_id = '+str(taxid)+')'
    try:
        superkingdom_id = Taxgroups.objects.raw(sql_query)[0]
    except:
        raise Http404("What is this species without a superkingdom?")
    return superkingdom_id.pk

def get_anno_strain(resid, residue_alignments, superk, assoc_or_struc):
    if str(superk) == '2157':
        ad_annotated_species = '186497'
    elif str(superk) == '2759':
        ad_annotated_species = '9606'
    elif str(superk) == '2':
        if assoc_or_struc == 'assoc':
            #ad_annotated_species = '262724'
            ad_annotated_species = '511145'
            print(assoc_or_struc)
        elif assoc_or_struc == 'struc':
            ad_annotated_species = '511145'
    else:
        raise Http404("This species is not part of either of the three superkingdoms?")
    
    #annoaln = Alignment.objects.filter(alndata__res = resid, source = 'abe')
    annoaln = Alignment.objects.filter(alndata__res = resid)
    if len(annoaln) == 0:
        return None
    annoaln = annoaln[0].aln_id
    sql_query = 'SELECT * FROM Aln_Data WHERE \
                    aln_id = '+str(annoaln)+' AND \
                    aln_pos = '+str(residue_alignments[annoaln])+' AND \
                    res_id in (SELECT resi_id FROM Residues WHERE PolData_id IN\
                        (SELECT PData_id FROM Polymer_Data WHERE strain_id = '+ad_annotated_species+'))'
    print('sqlq',sql_query)
    if len(AlnData.objects.raw(sql_query)) == 0:
        return None
    return AlnData.objects.raw(sql_query)[0].res_id

def related_ad_data_annotated_resis(resid, residue_alignments, superk):
    annotated_resi = get_anno_strain(resid, residue_alignments, superk, 'assoc')
    print('ar1', annotated_resi)
    if annotated_resi is None:
        return None
    ad_filter = AssociatedData.objects.filter(adresidues__residuep = annotated_resi)
    print('ar1, af', ad_filter)
    return ad_filter

def related_struc_fold_annotated_resis(resid, residue_alignments, superk):
    annotated_resi = get_anno_strain(resid, residue_alignments, superk, 'struc')
    print('ar', annotated_resi)
    if annotated_resi is None:
        return None
    struc_filter = StructuralFolds.objects.filter(strucfoldresidues__residue_id = annotated_resi)
    return struc_filter

def recurse_get_fold_lineage(fold_id):
    '''Finish when we transition to Apollo2 (needs MYSQL 8)'''
    pass

def aln_info(request, aln_id, aln_length):
    aln_id_to_strain_id_map = {
        38 : "272569",
        39 : "511145"
    }
    ad_annotated_species = aln_id_to_strain_id_map[aln_id]
    aln_info_results = {}
    #sql_query = f'Select Type, Value, aln_pos from DESIRE.AD_Residues join DSIRE.Associated_Data on DESIRE.Associated_Data.Data_id = DESIRE.AD_Residues.AD_id JOIN (SELECT res_id, aln_pos, Aln_Data_id FROM Aln_Data AS View WHERE aln_id = {aln_id} AND aln_pos >= 1 AND aln_pos <= {aln_length} AND res_id in (SELECT resi_id FROM Residues WHERE PolData_id IN (SELECT PData_id FROM Polymer_Data WHERE strain_id = {ad_annotated_species}))) ON residueP_id = res_id'
    sql_query = f'SELECT Type, Value, aln_pos FROM DESIRE.AD_Residues JOIN DESIRE.Associated_Data ON DESIRE.Associated_Data.Data_id = DESIRE.AD_Residues.AD_id JOIN (SELECT res_id, aln_pos, Aln_Data_id FROM Aln_Data WHERE aln_id = {aln_id} AND res_id IN (SELECT resi_id FROM Residues WHERE PolData_id IN (SELECT PData_id FROM Polymer_Data WHERE strain_id = {ad_annotated_species}))) AS Subquery ON Subquery.res_id = DESIRE.AD_Residues.residueP_id'
    #print(sql_query)
    cursor = connection.cursor()
    cursor.execute(sql_query)
    results = cursor.fetchall()
    for _type, value, aln_pos in results:
        if not (aln_pos in aln_info_results):
            aln_info_results[aln_pos] = []
        aln_info_results[aln_pos].append({ "type" : _type, "value" : value })
    print("Returning now.")
    return JsonResponse(aln_info_results)

def resi_info(request, resi_id):
    try:
        residue = Residues.objects.get(pk=resi_id)
    except Residues.DoesNotExist:
        raise Http404("Residue id is not present in the database!")
    
    polymer = residue.poldata_id
    taxid = get_taxid_from_polid(polymer)
    resnum = residue.resnum
    unmodresname = residue.unmodresname
    superk = get_superkingdom_id(taxid)

    residue_alignments = dict()
    for aln in get_alignments_from_resid(resi_id):
        residue_alignments[aln.aln_id] = AlnData.objects.filter(aln = aln.aln_id, res = resi_id)[0].aln_pos

    residue_lf = list()
    residue_folds = StructuralFolds.objects.filter(strucfoldresidues__residue = resi_id)
    if len(residue_folds) == 0:
        residue_folds = related_struc_fold_annotated_resis(resi_id, residue_alignments, superk)
    if residue_folds is not None:
        for fold in residue_folds:
            residue_lf.append((fold.level,fold.name))
    
    assoc_data = list()
    ad_filter = AssociatedData.objects.filter(adresidues__residuep = resi_id)
    #print('adf',ad_filter)
    if len(ad_filter) == 0:            #In the case of no associated data check the data for aligned resi in an annotated polymer
        ad_filter = related_ad_data_annotated_resis(resi_id, residue_alignments, superk)
        #print('adf2',ad_filter)
    if ad_filter is not None:
        for ass_data in ad_filter:
            assoc_data.append((ass_data.type, ass_data.value))

    resi_data = {
        'Residue' : unmodresname,
        'Sequence position' : resnum,
        'Associated data' : assoc_data,
        'Polymer id' : polymer,
        'Taxonomy id' : taxid,
        'Superkingdom id' : superk,
        'Structural fold' : residue_lf,
        'Alignment positions' : residue_alignments,
    }

    if request == None:
        return resi_data
    else:
        return JsonResponse(resi_data, safe = False)
