from django.shortcuts import render
from rest_framework import viewsets
from url_filter.integrations.drf import DjangoFilterBackend

from .serializers import *
from .models import *

class SpeciesViewSet(viewsets.ModelViewSet):
    queryset = Species.objects.all().order_by('strain_id')
    serializer_class = SpeciesSerializer

class NomenclatureViewSet(viewsets.ModelViewSet):
    queryset = Nomenclature.objects.all().order_by('nom_id')
    serializer_class = NomenclatureSerializer
    filter_backends = [DjangoFilterBackend]
    filter_fields = ['nom_id', 'new_name', 'occurrence', 'moleculegroup']

class OldNomenclatureViewSet(viewsets.ModelViewSet):
    queryset = OldName.objects.all().order_by('old_id')
    serializer_class = OldNameSerializer
    filter_backends = [DjangoFilterBackend]
    filter_fields = ['old_id', 'nn_fk', 'old_name', 'n_b_y_h_a']

class PolymerViewSet(viewsets.ModelViewSet):
    queryset = PolymerData.objects.all().order_by('pdata_id')
    serializer_class = PolymerSerializer
    filter_backends = [DjangoFilterBackend]
    filter_fields = ['pdata_id', 'gi', 'genesymbol', 'genedescription', 'strain', 'nomgd']

class ResidueViewSet(viewsets.ModelViewSet):
    queryset = Residues.objects.all().order_by('resi_id')
    serializer_class = ResidueSerializer
    filter_backends = [DjangoFilterBackend]
    filter_fields = ['resnum', 'poldata', 'unmodresname']

class SSViewSet(viewsets.ModelViewSet):
    queryset = Secondarystructures.objects.all().order_by('secstr_id')
    serializer_class = SecondarystructuresSerializer

class SSDataViewSet(viewsets.ModelViewSet):
    queryset = SsData.objects.all().order_by('ssd_id')
    serializer_class = SSDataSerializer

class AdResiduesViewSet(viewsets.ModelViewSet):
    queryset = AdResidues.objects.raw('SELECT CONCAT(AD_Residues.AD_id,"_",AD_Residues.residueP_id) AS id,AD_id,residueP_id FROM SEREB.AD_Residues')
    serializer_class = AdResiduesSerializer

class AlignmentViewSet(viewsets.ModelViewSet):
    queryset = Alignment.objects.all().order_by('aln_id')
    serializer_class = AlignmentSerializer

class AssociatedDataViewSet(viewsets.ModelViewSet):
    queryset = AssociatedData.objects.all().order_by('data_id')
    serializer_class = AssociatedDataSerializer

class AlnDataViewSet(viewsets.ModelViewSet):
    queryset = AlnData.objects.all()
    serializer_class = AlnDataSerializer
    filter_backends = [DjangoFilterBackend]
    filter_fields = '__all__'

class TaxGroupViewSet(viewsets.ModelViewSet):
    queryset = Taxgroups.objects.all()
    serializer_class = TaxGroupSerializer
    filter_backends = [DjangoFilterBackend]
    filter_fields = ['groupname', 'grouplevel', 'taxgroup_id']

'''
/desire-api/residues/?resnum=10&poldata__strain=74969&poldata__nomgd__new_name=uL02
/desire-api/residue-alignment/?res=404342&aln=1
/desire-api/residue-alignment/?aln_pos=61&aln=1&res__poldata__strain__in=74969,9606
'''