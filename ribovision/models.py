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
    name = models.CharField(db_column='Name', max_length=45)  # Field name made lowercase.
    method = models.CharField(db_column='Method', max_length=45)  # Field name made lowercase.
    source = models.CharField(db_column='Source', max_length=10)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'Alignment'


class AlnData(models.Model):
    alndata_id = models.AutoField(db_column='AlnData_id', primary_key=True)  # Field name made lowercase.
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


class Chainlist(models.Model):
    chainlist_id = models.IntegerField(db_column='ChainList_id', primary_key=True)  # Field name made lowercase.
    number_3d_structure = models.ForeignKey('Threedstructures', models.DO_NOTHING, db_column='3D_structure_id')  # Field name made lowercase. Field renamed because it wasn't a valid Python identifier.
    polymer = models.ForeignKey('PolymerData', models.DO_NOTHING)
    chainname = models.CharField(db_column='ChainName', max_length=3, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'ChainList'


class Default3Dstructure(models.Model):
    secondary_structure = models.ForeignKey('Secondarystructures', models.DO_NOTHING)
    number_3d_structure = models.ForeignKey('Threedstructures', models.DO_NOTHING, db_column='3D_structure_id')  # Field name made lowercase. Field renamed because it wasn't a valid Python identifier.

    class Meta:
        managed = False
        db_table = 'Default3DStructure'


class Interactions(models.Model):
    interactions_id = models.IntegerField(primary_key=True)
    residue_i = models.ForeignKey('Residues', models.DO_NOTHING, db_column='residue_i', blank=True, null=True, related_name='i_to_j')
    residue_j = models.ForeignKey('Residues', models.DO_NOTHING, db_column='residue_j', blank=True, null=True, related_name='j_to_i')
    bp_type = models.CharField(max_length=5, blank=True, null=True)
    bp_group = models.CharField(max_length=13, blank=True, null=True)
    number_3d_structure = models.ForeignKey('Threedstructures', models.DO_NOTHING, db_column='3D_structure_id', blank=True, null=True)  # Field name made lowercase. Field renamed because it wasn't a valid Python identifier.

    class Meta:
        managed = False
        db_table = 'Interactions'


class Linelabels(models.Model):
    linelabel_id = models.IntegerField(db_column='LineLabel_id', primary_key=True)  # Field name made lowercase.
    x1 = models.FloatField(db_column='X1', blank=True, null=True)  # Field name made lowercase.
    y1 = models.FloatField(db_column='Y1', blank=True, null=True)  # Field name made lowercase.
    x2 = models.FloatField(db_column='X2', blank=True, null=True)  # Field name made lowercase.
    y2 = models.FloatField(db_column='Y2', blank=True, null=True)  # Field name made lowercase.
    fill = models.CharField(db_column='Fill', max_length=7, blank=True, null=True)  # Field name made lowercase.
    stroke = models.CharField(db_column='Stroke', max_length=7, blank=True, null=True)  # Field name made lowercase.
    strokewidth = models.FloatField(db_column='StrokeWidth', blank=True, null=True)  # Field name made lowercase.
    strokelinejoin = models.CharField(db_column='StrokeLineJoin', max_length=5, blank=True, null=True)  # Field name made lowercase.
    strokemiterlimit = models.FloatField(db_column='StrokeMiterLimit', blank=True, null=True)  # Field name made lowercase.
    secondary_structure = models.ForeignKey('Secondarystructures', models.DO_NOTHING)

    class Meta:
        managed = False
        db_table = 'LineLabels'

class Mastertable(models.Model):
    master_id = models.IntegerField()
    Active = models.IntegerField(db_column='Active', blank=True, null=True)  # Field name made lowercase.
    SpeciesName = models.CharField(db_column='SpeciesName', max_length=24, blank=True, null=True)  # Field name made lowercase.
    DataSetType = models.CharField(db_column='DataSetType', max_length=33, blank=True, null=True)  # Field name made lowercase.
    StructureName = models.CharField(db_column='StructureName', max_length=14, blank=True, null=True)  # Field name made lowercase.
    LoadString = models.CharField(db_column='LoadString', max_length=29, blank=True, null=True)  # Field name made lowercase.
    Species_Abr = models.CharField(db_column='Species_Abr', max_length=5, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'MasterTable'

class Moleculenames(models.Model):
    moleculename = models.CharField(db_column='MoleculeName', primary_key=True, max_length=6)  # Field name made lowercase.
    moleculetype = models.CharField(db_column='MoleculeType', max_length=70, blank=True, null=True)  # Field name made lowercase.
    moleculegroup = models.CharField(db_column='MoleculeGroup', max_length=5, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'MoleculeNames'

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

    class Meta:
        managed = False
        db_table = 'Polymer_metadata'


class Residues(models.Model):
    resi_id = models.AutoField(primary_key=True)
    poldata = models.ForeignKey(PolymerData, models.DO_NOTHING, db_column='PolData_id')  # Field name made lowercase.
    resnum = models.IntegerField(db_column='resNum')  # Field name made lowercase.
    unmodresname = models.CharField(db_column='unModResName', max_length=1)  # Field name made lowercase.
    modresname = models.CharField(db_column='modResName', max_length=1, blank=True, null=True)  # Field name made lowercase.
    altname = models.CharField(db_column='altName', max_length=45, blank=True, null=True)  # Field name made lowercase.

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


class SecondaryTertiary(models.Model):
    secondary_tertiary_id = models.IntegerField(primary_key=True)
    secondary_structure = models.ForeignKey(Secondarystructures, models.DO_NOTHING)
    number_3d_structure = models.ForeignKey('Threedstructures', models.DO_NOTHING, db_column='3D_structure_id')  # Field name made lowercase. Field renamed because it wasn't a valid Python identifier.

    class Meta:
        managed = False
        db_table = 'Secondary_Tertiary'


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


class Structdatamenu(models.Model):
    structdatamenu_id = models.IntegerField(db_column='StructDataMenu_id', primary_key=True)  # Field name made lowercase.
    number_3d_structure = models.ForeignKey('Threedstructures', models.DO_NOTHING, db_column='3D_structure_id')  # Field name made lowercase. Field renamed because it wasn't a valid Python identifier.
    struct_data = models.ForeignKey('Structdatamenudetails', models.DO_NOTHING)

    class Meta:
        managed = False
        db_table = 'StructDataMenu'


class Structdatamenudetails(models.Model):
    struct_data_id = models.IntegerField(primary_key=True)
    structdataname = models.CharField(db_column='StructDataName', max_length=100, blank=True, null=True)  # Field name made lowercase.
    variablename = models.CharField(db_column='VariableName', max_length=15, blank=True, null=True)  # Field name made lowercase.
    colorlist = models.CharField(db_column='ColorList', max_length=12, blank=True, null=True)  # Field name made lowercase.
    indexmode = models.CharField(db_column='IndexMode', max_length=5, blank=True, null=True)  # Field name made lowercase.
    extraarg = models.CharField(db_column='ExtraArg', max_length=9, blank=True, null=True)  # Field name made lowercase.
    description = models.CharField(db_column='Description', max_length=255, blank=True, null=True)  # Field name made lowercase.
    helplink = models.CharField(db_column='HelpLink', max_length=45, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'StructDataMenuDetails'


class Structuraldata2(models.Model):
    map_index = models.IntegerField(blank=True, null=True)
    domain_rn = models.CharField(db_column='Domain_RN', max_length=4, blank=True, null=True)  # Field name made lowercase.
    domain_an = models.IntegerField(db_column='Domain_AN', blank=True, null=True)  # Field name made lowercase.
    domains_color = models.IntegerField(db_column='Domains_Color', blank=True, null=True)  # Field name made lowercase.
    helix_num = models.CharField(db_column='Helix_Num', max_length=4, blank=True, null=True)  # Field name made lowercase.
    helix_color = models.IntegerField(db_column='Helix_Color', blank=True, null=True)  # Field name made lowercase.
    secondary_structure = models.ForeignKey(Secondarystructures, models.DO_NOTHING, blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'StructuralData2'


class Structuraldata3(models.Model):
    map_index = models.IntegerField(blank=True, null=True)
    value = models.FloatField(db_column='Value', blank=True, null=True)  # Field name made lowercase.
    struct_data = models.ForeignKey(Structdatamenudetails, models.DO_NOTHING)
    number_3d_structure = models.ForeignKey('Threedstructures', models.DO_NOTHING, db_column='3D_structure_id')  # Field name made lowercase. Field renamed because it wasn't a valid Python identifier.

    class Meta:
        managed = False
        db_table = 'StructuralData3'


class Taxgroups(models.Model):
    taxgroup_id = models.IntegerField(primary_key=True)
    grouplevel = models.CharField(db_column='groupLevel', max_length=45, blank=True, null=True)  # Field name made lowercase.
    groupname = models.CharField(db_column='groupName', max_length=45, blank=True, null=True)  # Field name made lowercase.
    parent = models.ForeignKey('self', models.DO_NOTHING, db_column='parent', blank=True, null=True)

    class Meta:
        managed = False
        db_table = 'TaxGroups'


class Textlabels(models.Model):
    textlabel_id = models.IntegerField(db_column='TextLabel_id', primary_key=True)  # Field name made lowercase.
    labeltext = models.CharField(db_column='LabelText', max_length=500, blank=True, null=True)  # Field name made lowercase.
    x = models.FloatField(db_column='X', blank=True, null=True)  # Field name made lowercase.
    y = models.FloatField(db_column='Y', blank=True, null=True)  # Field name made lowercase.
    font = models.CharField(db_column='Font', max_length=100, blank=True, null=True)  # Field name made lowercase.
    font_size = models.FloatField(db_column='Font_Size', blank=True, null=True)  # Field name made lowercase.
    fill = models.CharField(db_column='Fill', max_length=50, blank=True, null=True)  # Field name made lowercase.
    secondary_structure = models.ForeignKey(Secondarystructures, models.DO_NOTHING)

    class Meta:
        managed = False
        db_table = 'TextLabels'


class Threedstructures(models.Model):
    number_3d_structure_id = models.AutoField(db_column='3D_structure_id', primary_key=True)  # Field name made lowercase. Field renamed because it wasn't a valid Python identifier.
    structurename = models.CharField(db_column='StructureName', max_length=50, blank=True, null=True)  # Field name made lowercase.

    class Meta:
        managed = False
        db_table = 'ThreeDStructures'


class AuthGroup(models.Model):
    name = models.CharField(unique=True, max_length=150)

    class Meta:
        managed = False
        db_table = 'auth_group'


class AuthGroupPermissions(models.Model):
    group = models.ForeignKey(AuthGroup, models.DO_NOTHING)
    permission = models.ForeignKey('AuthPermission', models.DO_NOTHING)

    class Meta:
        managed = False
        db_table = 'auth_group_permissions'
        unique_together = (('group', 'permission'),)


class AuthPermission(models.Model):
    name = models.CharField(max_length=255)
    content_type = models.ForeignKey('DjangoContentType', models.DO_NOTHING)
    codename = models.CharField(max_length=100)

    class Meta:
        managed = False
        db_table = 'auth_permission'
        unique_together = (('content_type', 'codename'),)


class AuthUser(models.Model):
    password = models.CharField(max_length=128)
    last_login = models.DateTimeField(blank=True, null=True)
    is_superuser = models.IntegerField()
    username = models.CharField(unique=True, max_length=150)
    first_name = models.CharField(max_length=30)
    last_name = models.CharField(max_length=150)
    email = models.CharField(max_length=254)
    is_staff = models.IntegerField()
    is_active = models.IntegerField()
    date_joined = models.DateTimeField()

    class Meta:
        managed = False
        db_table = 'auth_user'


class AuthUserGroups(models.Model):
    user = models.ForeignKey(AuthUser, models.DO_NOTHING)
    group = models.ForeignKey(AuthGroup, models.DO_NOTHING)

    class Meta:
        managed = False
        db_table = 'auth_user_groups'
        unique_together = (('user', 'group'),)


class AuthUserUserPermissions(models.Model):
    user = models.ForeignKey(AuthUser, models.DO_NOTHING)
    permission = models.ForeignKey(AuthPermission, models.DO_NOTHING)

    class Meta:
        managed = False
        db_table = 'auth_user_user_permissions'
        unique_together = (('user', 'permission'),)


class DjangoAdminLog(models.Model):
    action_time = models.DateTimeField()
    object_id = models.TextField(blank=True, null=True)
    object_repr = models.CharField(max_length=200)
    action_flag = models.SmallIntegerField()
    change_message = models.TextField()
    content_type = models.ForeignKey('DjangoContentType', models.DO_NOTHING, blank=True, null=True)
    user = models.ForeignKey(AuthUser, models.DO_NOTHING)

    class Meta:
        managed = False
        db_table = 'django_admin_log'


class DjangoContentType(models.Model):
    app_label = models.CharField(max_length=100)
    model = models.CharField(max_length=100)

    class Meta:
        managed = False
        db_table = 'django_content_type'
        unique_together = (('app_label', 'model'),)


class DjangoMigrations(models.Model):
    app = models.CharField(max_length=255)
    name = models.CharField(max_length=255)
    applied = models.DateTimeField()

    class Meta:
        managed = False
        db_table = 'django_migrations'


class DjangoSession(models.Model):
    session_key = models.CharField(primary_key=True, max_length=40)
    session_data = models.TextField()
    expire_date = models.DateTimeField()

    class Meta:
        managed = False
        db_table = 'django_session'
