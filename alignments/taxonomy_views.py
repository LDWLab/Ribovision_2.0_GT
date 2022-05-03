from django.http import JsonResponse
from alignments.models import *

# Recursively collect all children at level for given taxgroup_id
'''
query = '\
	with recursive cte (taxgroup_id, groupName, parent, groupLevel) as \
	(\
	select taxgroup_id, groupName, parent, groupLevel\
		from TaxGroups\
		where parent = 2157\
		union all\
		select p.taxgroup_id, p.groupName, p.parent, p.groupLevel\
		from TaxGroups p\
		inner join cte\
			on p.parent = cte.taxgroup_id\
	)\
	select taxgroup_id from cte where (groupLevel REGEXP ".*phylum")'
'''

def buildTaxonomy(request):
	taxgroups = Taxgroups.objects.raw('SELECT * FROM TaxGroups WHERE TaxGroups.groupLevel = "superkingdom";')
	taxonomy = []
	for taxgroup in taxgroups:
		subtaxonomy = {
			'label' : taxgroup.groupname,
			'children' : buildTaxonomyRecurse(taxgroup.taxgroup_id),
			'id' : taxgroup.taxgroup_id
		}
		taxonomy.append(subtaxonomy)
	tree = {
		'label' : 'Root',
		'children' : taxonomy,
		'id' : 0
	}
	return JsonResponse(tree, safe = False)

def api_showTaxonomy(request, current_tax):
	if current_tax == 0:
		toplevel_label = "Root"
		taxgroups = Taxgroups.objects.raw('SELECT * FROM TaxGroups WHERE TaxGroups.groupLevel = "superkingdom";')
		curent_parent = 0
		level = "root"
	else:
		toplevel_label = Taxgroups.objects.filter(taxgroup_id = current_tax)[0].groupname
		if Taxgroups.objects.get(pk=current_tax).grouplevel == "superkingdom":
			curent_parent = 0
		else:
			curent_parent = Taxgroups.objects.get(pk=current_tax).parent.pk
		level = Taxgroups.objects.get(pk=current_tax).grouplevel
		mySQLStr = 'SELECT * FROM TaxGroups WHERE TaxGroups.parent = "' + str(current_tax) + '";'
		taxgroups = Taxgroups.objects.raw(mySQLStr)
	nodes = list()
	for taxgroup in taxgroups:
		sql_for_children = f'SELECT * FROM TaxGroups WHERE TaxGroups.parent = "{taxgroup.taxgroup_id}"'
		# children = Taxgroups.objects.raw(sql_for_children)
		# children_for_json = list()
		# for child in children:
		# 	children_for_json.append({
		# 		'label': child.groupname,
		# 		'id': child.taxgroup_id,
		# 		'level': child.grouplevel,
		# 	})
		nodes.append({
			'label': taxgroup.groupname,
			'id': taxgroup.taxgroup_id,
			'level': taxgroup.grouplevel,
			'children': None
		})
	taxonomy = {
			'label' : toplevel_label,
			'id' : current_tax,
			'level' : level,
			'parent' : curent_parent,
			'children' : nodes
		}
	return JsonResponse(taxonomy, safe = False)

def buildTaxonomyRecurse(parentIndex):
	mySQLStr = 'SELECT * FROM TaxGroups WHERE TaxGroups.parent = "' + str(parentIndex) + '";'
	taxgroups = Taxgroups.objects.raw(mySQLStr)
	taxonomy = []
	for taxgroup in taxgroups:
		subtaxonomy = {
			'label' : taxgroup.groupname,
			'children' : buildTaxonomyRecurse(taxgroup.taxgroup_id),
			'id' : taxgroup.taxgroup_id
		}
		taxonomy.append(subtaxonomy)
	return taxonomy

def buildFoldTaxonomy(request):
	struc_groups = StructuralFolds.objects.raw('SELECT * FROM Structural_Folds WHERE Structural_Folds.Level = "Architecture";')
	taxonomy = []
	for taxgroup in struc_groups:
		subtaxonomy = {
			'label' : taxgroup.name,
			'children' : buildFoldTaxonomyRecurse(taxgroup.struc_fold_id),
			'id' : taxgroup.struc_fold_id,
			'ext_id' : taxgroup.external_id
		}
		taxonomy.append(subtaxonomy)
	tree = {
		'label' : 'Root',
		'children' : taxonomy,
		'id' : 0,
		'ext_id' : None
	}
	return JsonResponse(tree, safe = False)

def buildFoldTaxonomyRecurse(parentIndex):
	mySQLStr = 'SELECT * FROM Structural_Folds WHERE Structural_Folds.parent = "' + str(parentIndex) + '";'
	struc_groups = StructuralFolds.objects.raw(mySQLStr)
	taxonomy = []
	for curr_group in struc_groups:
		subtaxonomy = {
			'label' : curr_group.name,
			'children' : buildFoldTaxonomyRecurse(curr_group.struc_fold_id),
			'id' : curr_group.struc_fold_id,
			'ext_id' : curr_group.external_id
		}
		taxonomy.append(subtaxonomy)
	return taxonomy