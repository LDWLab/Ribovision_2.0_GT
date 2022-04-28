from django.shortcuts import render
from django.http import QueryDict, JsonResponse, Http404
from django.db import connection
from rest_framework import viewsets
from rest_framework.permissions import IsAuthenticated
from url_filter.integrations.drf import DjangoFilterBackend
from url_filter.filtersets import ModelFilterSet

from .serializers import *
from .models import *
from alignments.alignment_query_and_build import construct_query, dictfetchall


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
    filter_fields = ['pdata_id', 'gi', 'genesymbol', 'genedescription', 'strain', 'nomgd', 'alns_of_polymer']

class ResidueViewSet(viewsets.ModelViewSet):
    queryset = Residues.objects.all()
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
    queryset = AdResidues.objects.raw('SELECT CONCAT(AD_Residues.AD_id,"_",AD_Residues.residueP_id) AS id,AD_id,residueP_id FROM DESIRE.AD_Residues')
    serializer_class = AdResiduesSerializer

class AlignmentViewSet(viewsets.ModelViewSet):
    queryset = Alignment.objects.all().order_by('aln_id')
    serializer_class = AlignmentSerializer
    filter_backends = [DjangoFilterBackend]
    filter_fields = ['name', 'method', 'source', 'polymers']

class AssociatedDataViewSet(viewsets.ModelViewSet):
    queryset = AssociatedData.objects.all().order_by('data_id')
    serializer_class = AssociatedDataSerializer

class AlnDataViewSet(viewsets.ModelViewSet):
    queryset = AlnData.objects.all().order_by('aln_data_id')
    serializer_class = AlnDataSerializer
    filter_backends = [DjangoFilterBackend]
    filter_fields = '__all__'

class TaxGroupViewSet(viewsets.ModelViewSet):
    queryset = Taxgroups.objects.all().order_by('taxgroup_id')
    serializer_class = TaxGroupSerializer
    filter_backends = [DjangoFilterBackend]
    filter_fields = ['groupname', 'grouplevel', 'taxgroup_id']

class ResiFilterSet(ModelFilterSet):
    class Meta(object):
        model = Residues

class ResiAlnFilterSet(ModelFilterSet):
    class Meta(object):
        model = AlnData

class EcodDomainFilterSet(viewsets.ModelViewSet):
    queryset = Ecoddomains.objects.all().order_by('uid')
    serializer_class = EcoDDomainsSerializer
    filter_backends = [DjangoFilterBackend]
    #permission_classes = [IsAuthenticated]
    filter_fields = ['pdb', 'chain', 'arch_name', 'x_name', 'h_name', 't_name', 'f_name', 'ecod_domain_id']

def get_filter_set_results(query, filterset, queryset, field):
    curr_query = QueryDict(query)
    curr_fs = filterset(data=curr_query, queryset=queryset)
    for entry in curr_fs.filter():
        yield getattr(entry, field)

def filterresi(request, resnum, strain, new_name, aln_id, parent_tx):
    #Add check for strain to be in parent_tx

    resis = list(get_filter_set_results(f'resnum={resnum}&poldata__strain={strain}&poldata__nomgd__new_name={new_name}',
                                    ResiFilterSet, Residues.objects.all(), 'pk'))
    if len(resis) == 0:
        raise Http404(f'The combination of residue number {resnum}, strain id {strain}, and polymer name {new_name} is not present in the database!')
    if len(resis) > 1:
        raise Http404(f'The combination of residue number {resnum}, strain id {strain}, and polymer name {new_name} should filter down to a single residue!')
    
    alnpositions = list(get_filter_set_results(f'res={resis[0]}&aln={aln_id}',
                                            ResiAlnFilterSet, AlnData.objects.all().order_by('aln_data_id'), 'aln_pos'))
    if len(alnpositions) == 0:
        raise Http404(f'The combination of residue id {resis[0]}, and alignment id {aln_id} is not present in the database!\n\
    Likely the polymer name {new_name} does not correspond to alignment id {aln_id}.')

    SQLStatement = construct_query(aln_id, parent_tx)+' AND DESIRE.Aln_Data.aln_pos=%s'%(str(alnpositions[0]))
    with connection.cursor() as cursor:
        cursor.execute(SQLStatement)
        raw_result = dictfetchall(cursor)
    if len(raw_result) == 0:
        raise Http404("We do not have this combination of arguments in our database.")
    residues = list()
    for row in raw_result:
        residues.append(row['unModResName'])
    data = {
        "residue_id" : resis,
        "alignment position" : alnpositions,
        "residues in column" : residues
    }
    return JsonResponse(data, safe = False)
