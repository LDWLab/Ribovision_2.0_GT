from django.http import JsonResponse, Http404
from alignments.models import *
from alignments.residue_api import get_superkingdom_id, get_anno_strain

def pdbid_to_strainid(pdbid):
	'''Transforms PDBID to taxid'''
	structure_id = Threedstructures.objects.filter(structurename = pdbid)[0]
	secondary_id = SecondaryTertiary.objects.values("secondary_structure").filter(number_3d_structure = structure_id)[0]["secondary_structure"]
	anchor_taxid = Secondarystructures.objects.values("strain_fk").filter(secstr_id = secondary_id)[0]["strain_fk"]
	return anchor_taxid

def get_chains_polymers(struc_id):
	chain_pol_list = list()
	query = "SELECT new_name, GI, PData_id, ChainName, ChainList_id FROM DESIRE.Polymer_Data\
			INNER JOIN (SELECT * FROM DESIRE.ChainList WHERE 3D_structure_id = "+str(struc_id)+") \
				as filtered_chains \
				ON DESIRE.Polymer_Data.Pdata_id = filtered_chains.polymer_id\
			INNER JOIN Nomenclature ON DESIRE.Polymer_Data.nomgd_id = Nomenclature.nom_id"
	polymer_chains = PolymerData.objects.raw(query)
	if len(polymer_chains) == 0:
		return chain_pol_list
	return polymer_chains


def struc_info(request, struc_id):
	try:
		structure = Threedstructures.objects.get(pk=struc_id)
	except Threedstructures.DoesNotExist:
		raise Http404("Structure id "+str(struc_id)+" is not present in the database!")
	taxid = pdbid_to_strainid(structure.structurename)
	superk = get_superkingdom_id(taxid)
	ssid_name = list()
	pid_name_chid_name = list()
	for result in get_chains_polymers(struc_id):
		pid_name_chid_name.append((result.new_name, result.gi, result.pdata_id, result.ChainName, result.ChainList_id))

	struc_data = {
		'PDBID': structure.structurename,
		'Polymers and Chains': pid_name_chid_name,
		'Secondary Structures': ssid_name,
	}

	return JsonResponse(struc_data, safe = False)

