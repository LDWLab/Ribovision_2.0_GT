from django.http import JsonResponse
from alignments.models import *

def buildTaxonomy(request):
	taxgroups = Taxgroups.objects.raw('SELECT * FROM SEREB.TaxGroups WHERE SEREB.TaxGroups.groupLevel = "superkingdom";')
	taxonomy = []
	for taxgroup in taxgroups:
		subtaxonomy = {
			'label' : taxgroup.groupname,
			'nodes' : buildTaxonomyRecurse(taxgroup.taxgroup_id),
			'taxID' : taxgroup.taxgroup_id
		}
		taxonomy.append(subtaxonomy)
	tree = {
		'label' : 'Root',
		'nodes' : taxonomy,
		'taxID' : 0
	}
	return JsonResponse(tree, safe = False)

def api_showTaxonomy(request, current_tax):
	if current_tax == 0:
		toplevel_label = "LUCA"
		taxgroups = Taxgroups.objects.raw('SELECT * FROM SEREB.TaxGroups WHERE SEREB.TaxGroups.groupLevel = "superkingdom";')
		curent_parent = 0
		level = "root"
	else:
		toplevel_label = Taxgroups.objects.filter(taxgroup_id = current_tax)[0].groupname
		if Taxgroups.objects.get(pk=current_tax).grouplevel == "superkingdom":
			curent_parent = 0
		else:
			curent_parent = Taxgroups.objects.get(pk=current_tax).parent.pk
		level = Taxgroups.objects.get(pk=current_tax).grouplevel
		mySQLStr = 'SELECT * FROM SEREB.TaxGroups WHERE SEREB.TaxGroups.parent = "' + str(current_tax) + '";'
		taxgroups = Taxgroups.objects.raw(mySQLStr)
	nodes = list()
	for taxgroup in taxgroups:
		nodes.append((taxgroup.groupname,taxgroup.taxgroup_id,taxgroup.grouplevel))
	taxonomy = {
			'label' : toplevel_label,
			'taxID' : current_tax,
			'level' : level,
			'parent' : curent_parent,
			'nodes' : nodes
		}
	return JsonResponse(taxonomy, safe = False)

def buildTaxonomyRecurse(parentIndex):
	mySQLStr = 'SELECT * FROM SEREB.TaxGroups WHERE SEREB.TaxGroups.parent = "' + str(parentIndex) + '";'
	taxgroups = Taxgroups.objects.raw(mySQLStr)
	taxonomy = []
	for taxgroup in taxgroups:
		subtaxonomy = {
			'label' : taxgroup.groupname,
			'nodes' : buildTaxonomyRecurse(taxgroup.taxgroup_id),
			'taxID' : taxgroup.taxgroup_id
		}
		taxonomy.append(subtaxonomy)
	return taxonomy