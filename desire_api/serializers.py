from rest_framework import serializers
from .models import Residues, PolymerData, Species, Nomenclature, OldName, Secondarystructures, SsData

class NomenclatureSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Nomenclature
        fields = '__all__'

class OldNameSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = OldName
        fields = ['old_id', 'nn_fk', 'old_name', 'n_b_y_h_a']

class SpeciesSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Species
        fields = '__all__'

class PolymerSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = PolymerData
        fields = ['pdata_id', 'gi', 'genesymbol', 'genedescription', 'strain', 'nomgd']

class ResidueSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Residues
        fields = '__all__'

class SecondarystructuresSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Secondarystructures
        fields = ['secstr_id', 'moleculegroup', 'variation', 'strain_fk', 'name']

class SSDataSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = SsData
        fields = ['ssd_id', 'res']