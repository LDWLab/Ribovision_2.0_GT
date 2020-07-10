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

class ResidueSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Residues
        fields = '__all__'

class PolResidueSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Residues
        fields = ['url']

class PolAlnSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Alignment
        fields = ['url']

class SpeciesPolymerAlignmentSerializer(serializers.HyperlinkedModelSerializer):
    alns_of_polymer = PolAlnSerializer(many=True, read_only=True)
    class Meta:
        model = PolymerData
        fields = ['url', 'alns_of_polymer']

class SpeciesSerializer(serializers.HyperlinkedModelSerializer):
    polymers_of_species = SpeciesPolymerAlignmentSerializer(many=True, read_only=True)
    class Meta:
        model = Species
        fields = '__all__'

class PolymerSerializer(serializers.HyperlinkedModelSerializer):
    residues_in_polymer = PolResidueSerializer(many=True, read_only=True)
    alns_of_polymer = PolAlnSerializer(many=True, read_only=True)
    class Meta:
        model = PolymerData
        fields = ['pdata_id', 'gi', 'genesymbol', 'genedescription', 'strain', 'nomgd', 'alns_of_polymer', 'residues_in_polymer']

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

class TaxGroupSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Taxgroups
        fields = '__all__'