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

def api_showTaxonomy(request, parent):
	taxonomy = {
			'label' : Taxgroups.objects.filter(taxgroup_id = parent)[0].groupname,
			'nodes' : buildTaxonomyRecurse(parent),
			'taxID' : parent
		}
	tree = {
		'label' : 'Root',
		'nodes' : taxonomy,
		'taxID' : 0
	}
	return JsonResponse(tree, safe = False)

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