from django.urls import include, path
from rest_framework import routers
from . import views

router = routers.DefaultRouter()
router.register(r'nomenclature', views.NomenclatureViewSet)
router.register(r'old-nomenclatures', views.OldNomenclatureViewSet)
router.register(r'species', views.SpeciesViewSet)
router.register(r'polymers', views.PolymerViewSet)
router.register(r'residues', views.ResidueViewSet)
router.register(r'ss', views.SSViewSet)
router.register(r'ssdata', views.SSDataViewSet)
router.register(r'alignments', views.AlignmentViewSet)
router.register(r'adresi', views.AdResiduesViewSet)
router.register(r'ad', views.AssociatedDataViewSet)
router.register(r'residue-alignment', views.AlnDataViewSet)
router.register(r'taxonomic-groups', views.TaxGroupViewSet)
router.register(r'ECOD-domains', views.EcodDomainFilterSet)


# Wire up our API using automatic URL routing.
# Additionally, we include login URLs for the browsable API.
urlpatterns = [
    path('', include(router.urls)),
    path('api-auth/', include('rest_framework.urls', namespace='rest_framework')),
    path('filterresi/<int:resnum>/<int:strain>/<str:new_name>/<int:aln_id>/<int:parent_tx>', views.filterresi),
]

#Check this for list of values when adding multiple parents
#https://stackoverflow.com/questions/62371344/custom-django-url-path-converter-comma-separated-integers
