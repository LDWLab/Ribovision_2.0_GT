from rest_framework import serializers
from .models import *

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
        fields = '__all__'

class SSDataSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = SsData
        fields = '__all__'

class AdResiduesSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = AdResidues
        fields = ['ad', 'residuep']

class AlignmentSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Alignment
        fields = '__all__'

class AssociatedDataSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = AssociatedData
        fields = '__all__'

class AlnDataSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = AlnData
        fields = '__all__'