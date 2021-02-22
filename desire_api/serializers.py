from rest_framework import serializers
from .models import *

class NomenclatureSerializer(serializers.HyperlinkedModelSerializer):
    #nom_id = serializers.IntegerField(read_only=True)
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
        fields = ['url', 'strain', 'name', 'abbreviation', 'polymers_of_species']

class PolymerSerializer(serializers.HyperlinkedModelSerializer):
    #residues_in_polymer = PolResidueSerializer(many=True, read_only=True)
    alns_of_polymer = PolAlnSerializer(many=True, read_only=True)
    class Meta:
        model = PolymerData
        fields = ['pdata_id', 'gi', 'genesymbol', 'genedescription', 'strain', 'nomgd', 'alns_of_polymer']

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

class GeneDescriptionSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = PolymerData
        fields = ['url', 'gi', 'strain', 'nomgd', 'genesymbol', 'genedescription']

class AlignmentSerializer(serializers.HyperlinkedModelSerializer):
    polymers = GeneDescriptionSerializer(many=True, read_only=True)
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

class AlignmentTxGrpSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Alignment
        fields = ['url', 'name', 'method', 'source']

class EcoDDomainsSerializer(serializers.HyperlinkedModelSerializer):
    class Meta:
        model = Ecoddomains
        fields = '__all__'

class TaxGroupSerializer(serializers.HyperlinkedModelSerializer):
    alignment_ids = serializers.SerializerMethodField()
    def get_alignment_ids(self, obj):
        #if obj.grouplevel == 'strain':
        #    return list(["Implement case of strain."])
        rawsql = f"\
        SELECT DISTINCT Alignment.Aln_id, Alignment.name, Alignment.Method FROM Alignment\
        INNER JOIN Polymer_Alignments ON Alignment.Aln_id = Polymer_Alignments.Aln_id\
        INNER JOIN Polymer_Data on Polymer_Alignments.PData_id = Polymer_Data.PData_id\
        WHERE Polymer_Data.strain_id IN\
        (with recursive cte (taxgroup_id, groupName, parent, groupLevel) as \
        (\
        select taxgroup_id, groupName, parent, groupLevel\
            from TaxGroups\
            where parent = {str(obj.taxgroup_id)} or taxgroup_id = {str(obj.taxgroup_id)}\
            union all\
            select p.taxgroup_id, p.groupName, p.parent, p.groupLevel\
            from TaxGroups p\
            inner join cte\
                on p.parent = cte.taxgroup_id\
        )\
        select taxgroup_id from cte where (groupLevel REGEXP 'strain')\
        )"
        aln_ids = Alignment.objects.raw(rawsql)
        outlist = [(x.aln_id, x.name, x.method) for x in aln_ids]
        return list(outlist)
    class Meta:
        model = Taxgroups
        fields = ['url', 'taxgroup_id','groupname', 'grouplevel', 'parent', 'alignment_ids']