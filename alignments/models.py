# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Make sure each ForeignKey has `on_delete` set to the desired behavior.
#   * Remove `managed = False` lines if you wish to allow Django to create, modify, and delete the table
# Feel free to rename the models, but don't rename db_table values or field names.
from django.db import models


class AdResidues(models.Model):
    ad = models.ForeignKey('AssociatedData', models.DO_NOTHING, db_column='AD_id')  # Field name made lowercase.
    residuep = models.ForeignKey('Residues', models.DO_NOTHING, db_column='residueP_id')  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'AD_Residues'
        unique_together = (('ad', 'residuep'),)


class Alignment(models.Model):
    aln_id = models.AutoField(db_column='Aln_id', primary_key=True)  # Field name made lowercase.
    name = models.CharField(db_column='Name', max_length=45)  # Make sure it's unique to make nice URLs # Field name made lowercase.
    method = models.CharField(db_column='Method', max_length=45)  # Field name made lowercase.
    source = models.CharField(db_column='Source', max_length=10)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'Alignment'


class AlnData(models.Model):
    aln_data_id = models.AutoField(db_column='Aln_Data_id', primary_key=True)  # Field name made lowercase.
    aln = models.ForeignKey(Alignment, models.DO_NOTHING)
    res = models.ForeignKey('Residues', models.DO_NOTHING)
    aln_pos = models.IntegerField()
    polymer_order = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'Aln_Data'


class AlnDomains(models.Model):
    dom_taxid = models.ForeignKey('Taxgroups', models.DO_NOTHING, db_column='dom_taxid')
    aln = models.ForeignKey(Alignment, models.DO_NOTHING)
    compartment = models.CharField(max_length=1)

    class Meta:
        managed = False
        db_table = 'Aln_Domains'
        unique_together = (('dom_taxid', 'aln'),)


class AssociatedData(models.Model):
    data_id = models.AutoField(db_column='Data_id', primary_key=True)  # Field name made lowercase.
    type = models.CharField(db_column='Type', max_length=45)  # Field name made lowercase.
    value = models.CharField(db_column='Value', max_length=45)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'Associated_Data'


class Nomenclature(models.Model):
    nom_id = models.AutoField(primary_key=True)
    new_name = models.CharField(max_length=10)
    occurrence = models.CharField(max_length=1)
    moleculegroup = models.CharField(db_column='MoleculeGroup', max_length=45, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'Nomenclature'


class OldName(models.Model):
    old_id = models.AutoField(primary_key=True)
    nn_fk = models.ForeignKey(Nomenclature, models.DO_NOTHING, db_column='nn_fk_id')
    old_name = models.CharField(max_length=45)
    n_b_y_h_a = models.CharField(db_column='N_B_Y_H_A', max_length=3)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'Old_name'


class PolymerData(models.Model):
    pdata_id = models.AutoField(db_column='PData_id', primary_key=True)  # Field name made lowercase.
    gi = models.CharField(db_column='GI', unique=True, max_length=45)  # Field name made lowercase.
    strain = models.ForeignKey('Species', models.DO_NOTHING)
    nomgd = models.ForeignKey(Nomenclature, models.DO_NOTHING, blank=True, null=True)
    genesymbol = models.CharField(db_column='GeneSymbol', max_length=45, blank=True, null=True)  # Field name made lowercase.
    genedescription = models.CharField(db_column='GeneDescription', max_length=100, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'Polymer_Data'


class PolymerMetadata(models.Model):
    polymer = models.OneToOneField(PolymerData, models.DO_NOTHING, primary_key=True)
    accession_type = models.CharField(max_length=45)
    polymer_type = models.CharField(max_length=45)
    accession = models.CharField(max_length=45, blank=True, null=True)
    fullseq = models.TextField(db_column='Fullseq', blank=True, null=True)  # Field name made lowercase.
    protein_kegg_type = models.CharField(db_column='Protein_KEGG_type', max_length=40, blank=True, null=True)  # Field name made lowercase.
    kegg_path_id = models.IntegerField(db_column='KEGG_path_id', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'Polymer_metadata'


class Residues(models.Model):
    resi_id = models.AutoField(primary_key=True)
    poldata = models.ForeignKey(PolymerData, models.DO_NOTHING, db_column='PolData_id')  # Field name made lowercase.
    resnum = models.IntegerField(db_column='resNum')  # Field name made lowercase.
    unmodresname = models.CharField(db_column='unModResName', max_length=1)  # Field name made lowercase.
    modresname = models.CharField(db_column='modResName', max_length=1, blank=True, null=True)  # Field name made lowercase.
    codon = models.CharField(db_column='codon', max_length=45, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'Residues'


class SsData(models.Model):
    ssd_id = models.AutoField(db_column='SSD_id', primary_key=True)  # Field name made lowercase.
    ss = models.ForeignKey('Secondarystructures', models.DO_NOTHING)
    res = models.ForeignKey(Residues, models.DO_NOTHING)
    map_index = models.IntegerField()
    x = models.FloatField(db_column='X', blank=True, null=True)  # Field name made lowercase.
    y = models.FloatField(db_column='Y', blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'SS_Data'


class Secondarystructures(models.Model):
    secstr_id = models.AutoField(db_column='SecStr_id', primary_key=True)  # Field name made lowercase.
    moleculegroup = models.CharField(db_column='MoleculeGroup', max_length=45)  # Field name made lowercase.
    variation = models.CharField(db_column='Variation', max_length=45)  # Field name made lowercase.
    strain_fk = models.ForeignKey('Species', models.DO_NOTHING, db_column='strain_fk', blank=True, null=True)
    name = models.CharField(db_column='Name', max_length=255, blank=True, null=True)  # Field name made lowercase.
    font_size_svg = models.DecimalField(db_column='Font_Size_SVG', max_digits=2, decimal_places=1, blank=True, null=True)  # Field name made lowercase.
    font_size_canvas = models.DecimalField(db_column='Font_Size_Canvas', max_digits=2, decimal_places=1, blank=True, null=True)  # Field name made lowercase.
    circle_radius = models.DecimalField(db_column='Circle_Radius', max_digits=2, decimal_places=1, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'SecondaryStructures'


class Species(models.Model):
    strain_id = models.IntegerField(primary_key=True)
    name = models.CharField(max_length=60, blank=True, null=True)
    strain = models.CharField(max_length=100, blank=True, null=True)
    taxid = models.IntegerField(blank=True, null=True)
    abbreviation = models.CharField(db_column='Abbreviation', max_length=45, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'Species'


class SpeciesTaxgroup(models.Model):
    strain = models.ForeignKey(Species, models.DO_NOTHING)
    taxgroup = models.ForeignKey('Taxgroups', models.DO_NOTHING)

    class Meta:
        managed = False
        db_table = 'Species_TaxGroup'
        unique_together = (('strain', 'taxgroup'),)


class Taxgroups(models.Model):
    taxgroup_id = models.IntegerField(primary_key=True)
    grouplevel = models.CharField(db_column='groupLevel', max_length=45, blank=True, null=True)  # Field name made lowercase.
    groupname = models.CharField(db_column='groupName', max_length=45, blank=True, null=True)  # Field name made lowercase.
    parent = models.ForeignKey('self', models.DO_NOTHING, db_column='parent', blank=True, null=True)
    # taxid = models.IntegerField(blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'TaxGroups'

class Threedstructures(models.Model):
    number_3d_structure_id = models.AutoField(db_column='3D_structure_id', primary_key=True)  # Field name made lowercase. Field renamed because it wasn't a valid Python identifier.
    structurename = models.CharField(db_column='StructureName', max_length=50, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'ThreeDStructures'

class Chainlist(models.Model):
    chainlist_id = models.IntegerField(db_column='ChainList_id', primary_key=True)  # Field name made lowercase.
    number_3d_structure = models.ForeignKey('Threedstructures', models.DO_NOTHING, db_column='3D_structure_id')  # Field name made lowercase. Field renamed because it wasn't a valid Python identifier.
    polymer = models.ForeignKey('PolymerData', models.DO_NOTHING)
    chainname = models.CharField(db_column='ChainName', max_length=3, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'ChainList'

class SecondaryTertiary(models.Model):
    secondary_tertiary_id = models.IntegerField(primary_key=True)
    secondary_structure = models.ForeignKey(Secondarystructures, models.DO_NOTHING)
    number_3d_structure = models.ForeignKey('Threedstructures', models.DO_NOTHING, db_column='3D_structure_id')  # Field name made lowercase. Field renamed because it wasn't a valid Python identifier.

    class Meta:
        managed = False
        db_table = 'Secondary_Tertiary'

class PolymerAlignments(models.Model):
    pdata = models.ForeignKey('PolymerData', models.DO_NOTHING, db_column='PData_id')  # Field name made lowercase.
    aln = models.ForeignKey(Alignment, models.DO_NOTHING, db_column='Aln_id')  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'Polymer_Alignments'
        unique_together = (('pdata', 'aln'),)

class StructuralFolds(models.Model):
    struc_fold_id = models.AutoField(primary_key=True)
    level = models.CharField(db_column='Level', max_length=45)  # Field name made lowercase.
    name = models.CharField(db_column='Name', max_length=45)  # Field name made lowercase.
    classification_system = models.CharField(max_length=45)
    parent = models.ForeignKey('self', models.DO_NOTHING, db_column='parent')
    external_id = models.CharField(max_length=45, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'Structural_Folds'

class StrucfoldResidues(models.Model):
    residue = models.ForeignKey(Residues, models.DO_NOTHING)
    strucfold = models.ForeignKey('StructuralFolds', models.DO_NOTHING)

    class Meta:
        managed = False
        db_table = 'StrucFold_Residues'
        unique_together = (('residue', 'strucfold'),)

class StrucfoldChains(models.Model):
    strucfold = models.ForeignKey('StructuralFolds', models.DO_NOTHING)
    chain = models.ForeignKey(Chainlist, models.DO_NOTHING)

    class Meta:
        managed = False
        db_table = 'StrucFold_Chains'
        unique_together = (('strucfold', 'chain'),)
