from django.shortcuts import render
from rest_framework import viewsets

from .serializers import *
from .models import Residues, PolymerData, Species, Nomenclature, OldName, Secondarystructures, SsData

class SpeciesViewSet(viewsets.ModelViewSet):
    queryset = Species.objects.all().order_by('strain_id')
    serializer_class = SpeciesSerializer

class NomenclatureViewSet(viewsets.ModelViewSet):
    queryset = Nomenclature.objects.all().order_by('nom_id')
    serializer_class = NomenclatureSerializer

class OldNomenclatureViewSet(viewsets.ModelViewSet):
    queryset = OldName.objects.all().order_by('old_id')
    serializer_class = OldNameSerializer

class PolymerViewSet(viewsets.ModelViewSet):
    queryset = PolymerData.objects.all().order_by('pdata_id')
    serializer_class = PolymerSerializer

class ResidueViewSet(viewsets.ModelViewSet):
    queryset = Residues.objects.all().order_by('resi_id')
    serializer_class = ResidueSerializer

class SSViewSet(viewsets.ModelViewSet):
    queryset = Secondarystructures.objects.all().order_by('secstr_id')
    serializer_class = SecondarystructuresSerializer

class SSDataViewSet(viewsets.ModelViewSet):
    queryset = SsData.objects.all().order_by('ssd_id')
    serializer_class = SSDataSerializer