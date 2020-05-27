from django.http import JsonResponse, Http404
from alignments.models import *

'''Recursion for transition to Apollo2:
with recursive cte (struc_fold_id, name, parent, level) as (
  select     struc_fold_id,
             name,
             parent,
             level
  from       Structural_Folds
  where      parent = 2
  union all
  select     p.struc_fold_id,
             p.name,
             p.parent,
             p.level
  from       Structural_Folds p
  inner join cte
          on p.parent = cte.struc_fold_id
)
select * from cte where level = 'F';
'''

def get_lowest_level_folds(fold_id, reqest_fold):
	if reqest_fold.level == 'F':
		return [reqest_fold]
	empty_fold_list = list()
	query = "SELECT  struc_fold_id, name, parent, level \
	FROM	(SELECT * from Structural_Folds\
			ORDER BY parent, struc_fold_id) folds_sorted,\
			(select @pv := '"+str(fold_id)+"') initialisation\
	WHERE 	find_in_set(parent, @pv)\
	AND 	length(@pv := CONCAT(@pv, ',', struc_fold_id))\
	AND 	level = 'F'"
	lower_folds = StructuralFolds.objects.raw(query)
	if len(lower_folds) == 0:
		return empty_fold_list
	return lower_folds


#CONCAT(Aln_id,"_",polymer_id) AS id
def get_available_alns(f_fold_id):
	query = 'SELECT Alignment.Aln_id,Alignment.Name,Alignment.Method,Polymer_Data.PData_id FROM SEREB.Polymer_Data\
			INNER JOIN SEREB.ChainList ON Polymer_Data.PData_id = ChainList.polymer_id\
			INNER JOIN SEREB.Polymer_Alignments ON Polymer_Data.PData_id = Polymer_Alignments.PData_id\
			INNER JOIN SEREB.Alignment ON Polymer_Alignments.Aln_id = Alignment.Aln_id\
			INNER JOIN (SELECT * FROM SEREB.StrucFold_Chains WHERE strucfold_id ='+str(f_fold_id)+')\
				AS filtered_chains\
			ON SEREB.ChainList.ChainList_id = filtered_chains.chain_id'
	raw_object = PolymerData.objects.raw(query)
	output_dict = dict()
	for aln_and_pid in raw_object:
		if aln_and_pid.pdata_id not in output_dict.keys():
			output_dict[aln_and_pid.pdata_id] = []
		output_dict[aln_and_pid.pdata_id].append((aln_and_pid.Aln_id,aln_and_pid.Name,aln_and_pid.Method))
	return output_dict

def fold_info(request, fold_id):
	try:
		reqest_fold = StructuralFolds.objects.get(pk=fold_id)
	except StructuralFolds.DoesNotExist:
		raise Http404("Fold id "+str(fold_id)+" is not present in the database!")

	if reqest_fold.classification_system != 'ECOD':
		raise Http404("Fold id "+str(fold_id)+" is in a classification system ("+reqest_fold.classification_system+") not yet supported!")

	fold_pol_aln = dict()
	for f_fold in get_lowest_level_folds(fold_id, reqest_fold):
		polymers_to_alns = get_available_alns(f_fold.struc_fold_id)
		if f_fold.struc_fold_id not in fold_pol_aln.keys():
			fold_pol_aln[f_fold.struc_fold_id] = dict()
		fold_pol_aln[f_fold.struc_fold_id] = polymers_to_alns

	fold_data = {
		'Fold name': reqest_fold.name,
		'Folds to polymers to alignments': fold_pol_aln
	}

	return JsonResponse(fold_data, safe = False)

