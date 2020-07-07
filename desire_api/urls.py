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

# Wire up our API using automatic URL routing.
# Additionally, we include login URLs for the browsable API.
urlpatterns = [
    path('', include(router.urls)),
    path('api-auth/', include('rest_framework.urls', namespace='rest_framework'))
]