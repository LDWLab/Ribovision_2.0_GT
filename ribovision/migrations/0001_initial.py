# Generated by Django 2.1.3 on 2019-10-04 16:22

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='AdResidues',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'managed': False,
                'db_table': 'AD_Residues',
            },
        ),
        migrations.CreateModel(
            name='Alignment',
            fields=[
                ('aln_id', models.AutoField(db_column='Aln_id', primary_key=True, serialize=False)),
                ('name', models.CharField(db_column='Name', max_length=45)),
                ('method', models.CharField(db_column='Method', max_length=45)),
                ('source', models.CharField(db_column='Source', max_length=10)),
            ],
            options={
                'managed': False,
                'db_table': 'Alignment',
            },
        ),
        migrations.CreateModel(
            name='AlnData',
            fields=[
                ('alndata_id', models.AutoField(db_column='AlnData_id', primary_key=True, serialize=False)),
                ('aln_pos', models.IntegerField()),
                ('polymer_order', models.IntegerField(blank=True, null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'Aln_Data',
            },
        ),
        migrations.CreateModel(
            name='AlnDomains',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('compartment', models.CharField(max_length=1)),
            ],
            options={
                'managed': False,
                'db_table': 'Aln_Domains',
            },
        ),
        migrations.CreateModel(
            name='AssociatedData',
            fields=[
                ('data_id', models.AutoField(db_column='Data_id', primary_key=True, serialize=False)),
                ('type', models.CharField(db_column='Type', max_length=45)),
                ('value', models.CharField(db_column='Value', max_length=45)),
            ],
            options={
                'managed': False,
                'db_table': 'Associated_Data',
            },
        ),
        migrations.CreateModel(
            name='AuthGroup',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=150, unique=True)),
            ],
            options={
                'managed': False,
                'db_table': 'auth_group',
            },
        ),
        migrations.CreateModel(
            name='AuthGroupPermissions',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'managed': False,
                'db_table': 'auth_group_permissions',
            },
        ),
        migrations.CreateModel(
            name='AuthPermission',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=255)),
                ('codename', models.CharField(max_length=100)),
            ],
            options={
                'managed': False,
                'db_table': 'auth_permission',
            },
        ),
        migrations.CreateModel(
            name='AuthUser',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('password', models.CharField(max_length=128)),
                ('last_login', models.DateTimeField(blank=True, null=True)),
                ('is_superuser', models.IntegerField()),
                ('username', models.CharField(max_length=150, unique=True)),
                ('first_name', models.CharField(max_length=30)),
                ('last_name', models.CharField(max_length=150)),
                ('email', models.CharField(max_length=254)),
                ('is_staff', models.IntegerField()),
                ('is_active', models.IntegerField()),
                ('date_joined', models.DateTimeField()),
            ],
            options={
                'managed': False,
                'db_table': 'auth_user',
            },
        ),
        migrations.CreateModel(
            name='AuthUserGroups',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'managed': False,
                'db_table': 'auth_user_groups',
            },
        ),
        migrations.CreateModel(
            name='AuthUserUserPermissions',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'managed': False,
                'db_table': 'auth_user_user_permissions',
            },
        ),
        migrations.CreateModel(
            name='Chainlist',
            fields=[
                ('chainlist_id', models.IntegerField(db_column='ChainList_id', primary_key=True, serialize=False)),
                ('chainname', models.CharField(blank=True, db_column='ChainName', max_length=3, null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'ChainList',
            },
        ),
        migrations.CreateModel(
            name='Default3Dstructure',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'managed': False,
                'db_table': 'Default3DStructure',
            },
        ),
        migrations.CreateModel(
            name='DjangoAdminLog',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('action_time', models.DateTimeField()),
                ('object_id', models.TextField(blank=True, null=True)),
                ('object_repr', models.CharField(max_length=200)),
                ('action_flag', models.SmallIntegerField()),
                ('change_message', models.TextField()),
            ],
            options={
                'managed': False,
                'db_table': 'django_admin_log',
            },
        ),
        migrations.CreateModel(
            name='DjangoContentType',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('app_label', models.CharField(max_length=100)),
                ('model', models.CharField(max_length=100)),
            ],
            options={
                'managed': False,
                'db_table': 'django_content_type',
            },
        ),
        migrations.CreateModel(
            name='DjangoMigrations',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('app', models.CharField(max_length=255)),
                ('name', models.CharField(max_length=255)),
                ('applied', models.DateTimeField()),
            ],
            options={
                'managed': False,
                'db_table': 'django_migrations',
            },
        ),
        migrations.CreateModel(
            name='DjangoSession',
            fields=[
                ('session_key', models.CharField(max_length=40, primary_key=True, serialize=False)),
                ('session_data', models.TextField()),
                ('expire_date', models.DateTimeField()),
            ],
            options={
                'managed': False,
                'db_table': 'django_session',
            },
        ),
        migrations.CreateModel(
            name='Interactions',
            fields=[
                ('interactions_id', models.IntegerField(primary_key=True, serialize=False)),
                ('bp_type', models.CharField(blank=True, max_length=5, null=True)),
                ('bp_group', models.CharField(blank=True, max_length=13, null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'Interactions',
            },
        ),
        migrations.CreateModel(
            name='Linelabels',
            fields=[
                ('linelabel_id', models.IntegerField(db_column='LineLabel_id', primary_key=True, serialize=False)),
                ('x1', models.FloatField(blank=True, db_column='X1', null=True)),
                ('y1', models.FloatField(blank=True, db_column='Y1', null=True)),
                ('x2', models.FloatField(blank=True, db_column='X2', null=True)),
                ('y2', models.FloatField(blank=True, db_column='Y2', null=True)),
                ('fill', models.CharField(blank=True, db_column='Fill', max_length=7, null=True)),
                ('stroke', models.CharField(blank=True, db_column='Stroke', max_length=7, null=True)),
                ('strokewidth', models.FloatField(blank=True, db_column='StrokeWidth', null=True)),
                ('strokelinejoin', models.CharField(blank=True, db_column='StrokeLineJoin', max_length=5, null=True)),
                ('strokemiterlimit', models.FloatField(blank=True, db_column='StrokeMiterLimit', null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'LineLabels',
            },
        ),
        migrations.CreateModel(
            name='Nomenclature',
            fields=[
                ('nom_id', models.AutoField(primary_key=True, serialize=False)),
                ('new_name', models.CharField(max_length=10)),
                ('occurrence', models.CharField(max_length=1)),
                ('moleculegroup', models.CharField(blank=True, db_column='MoleculeGroup', max_length=45, null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'Nomenclature',
            },
        ),
        migrations.CreateModel(
            name='OldName',
            fields=[
                ('old_id', models.AutoField(primary_key=True, serialize=False)),
                ('old_name', models.CharField(max_length=45)),
                ('n_b_y_h_a', models.CharField(db_column='N_B_Y_H_A', max_length=3)),
            ],
            options={
                'managed': False,
                'db_table': 'Old_name',
            },
        ),
        migrations.CreateModel(
            name='PolymerData',
            fields=[
                ('pdata_id', models.AutoField(db_column='PData_id', primary_key=True, serialize=False)),
                ('gi', models.CharField(db_column='GI', max_length=45, unique=True)),
                ('genesymbol', models.CharField(blank=True, db_column='GeneSymbol', max_length=45, null=True)),
                ('genedescription', models.CharField(blank=True, db_column='GeneDescription', max_length=100, null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'Polymer_Data',
            },
        ),
        migrations.CreateModel(
            name='Residues',
            fields=[
                ('resi_id', models.AutoField(primary_key=True, serialize=False)),
                ('resnum', models.IntegerField(db_column='resNum')),
                ('unmodresname', models.CharField(db_column='unModResName', max_length=1)),
                ('modresname', models.CharField(blank=True, db_column='modResName', max_length=1, null=True)),
                ('altname', models.CharField(blank=True, db_column='altName', max_length=45, null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'Residues',
            },
        ),
        migrations.CreateModel(
            name='Secondarystructures',
            fields=[
                ('secstr_id', models.AutoField(db_column='SecStr_id', primary_key=True, serialize=False)),
                ('moleculegroup', models.CharField(db_column='MoleculeGroup', max_length=45)),
                ('variation', models.CharField(db_column='Variation', max_length=45)),
                ('name', models.CharField(blank=True, db_column='Name', max_length=255, null=True)),
                ('font_size_svg', models.DecimalField(blank=True, db_column='Font_Size_SVG', decimal_places=1, max_digits=2, null=True)),
                ('font_size_canvas', models.DecimalField(blank=True, db_column='Font_Size_Canvas', decimal_places=1, max_digits=2, null=True)),
                ('circle_radius', models.DecimalField(blank=True, db_column='Circle_Radius', decimal_places=1, max_digits=2, null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'SecondaryStructures',
            },
        ),
        migrations.CreateModel(
            name='SecondaryTertiary',
            fields=[
                ('secondary_tertiary_id', models.IntegerField(primary_key=True, serialize=False)),
            ],
            options={
                'managed': False,
                'db_table': 'Secondary_Tertiary',
            },
        ),
        migrations.CreateModel(
            name='Species',
            fields=[
                ('strain_id', models.IntegerField(primary_key=True, serialize=False)),
                ('name', models.CharField(blank=True, max_length=60, null=True)),
                ('strain', models.CharField(blank=True, max_length=100, null=True)),
                ('taxid', models.IntegerField(blank=True, null=True)),
                ('abbreviation', models.CharField(blank=True, db_column='Abbreviation', max_length=45, null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'Species',
            },
        ),
        migrations.CreateModel(
            name='SpeciesTaxgroup',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
            ],
            options={
                'managed': False,
                'db_table': 'Species_TaxGroup',
            },
        ),
        migrations.CreateModel(
            name='SsData',
            fields=[
                ('ssd_id', models.AutoField(db_column='SSD_id', primary_key=True, serialize=False)),
                ('map_index', models.IntegerField()),
                ('x', models.FloatField(blank=True, db_column='X', null=True)),
                ('y', models.FloatField(blank=True, db_column='Y', null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'SS_Data',
            },
        ),
        migrations.CreateModel(
            name='Structdatamenu',
            fields=[
                ('structdatamenu_id', models.IntegerField(db_column='StructDataMenu_id', primary_key=True, serialize=False)),
            ],
            options={
                'managed': False,
                'db_table': 'StructDataMenu',
            },
        ),
        migrations.CreateModel(
            name='Structdatamenudetails',
            fields=[
                ('struct_data_id', models.IntegerField(primary_key=True, serialize=False)),
                ('structdataname', models.CharField(blank=True, db_column='StructDataName', max_length=100, null=True)),
                ('variablename', models.CharField(blank=True, db_column='VariableName', max_length=15, null=True)),
                ('colorlist', models.CharField(blank=True, db_column='ColorList', max_length=12, null=True)),
                ('indexmode', models.CharField(blank=True, db_column='IndexMode', max_length=5, null=True)),
                ('extraarg', models.CharField(blank=True, db_column='ExtraArg', max_length=9, null=True)),
                ('description', models.CharField(blank=True, db_column='Description', max_length=255, null=True)),
                ('helplink', models.CharField(blank=True, db_column='HelpLink', max_length=45, null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'StructDataMenuDetails',
            },
        ),
        migrations.CreateModel(
            name='Structuraldata2',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('map_index', models.IntegerField(blank=True, null=True)),
                ('domain_rn', models.CharField(blank=True, db_column='Domain_RN', max_length=4, null=True)),
                ('domain_an', models.IntegerField(blank=True, db_column='Domain_AN', null=True)),
                ('domains_color', models.IntegerField(blank=True, db_column='Domains_Color', null=True)),
                ('helix_num', models.CharField(blank=True, db_column='Helix_Num', max_length=4, null=True)),
                ('helix_color', models.IntegerField(blank=True, db_column='Helix_Color', null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'StructuralData2',
            },
        ),
        migrations.CreateModel(
            name='Structuraldata3',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('map_index', models.IntegerField(blank=True, null=True)),
                ('value', models.FloatField(blank=True, db_column='Value', null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'StructuralData3',
            },
        ),
        migrations.CreateModel(
            name='Taxgroups',
            fields=[
                ('taxgroup_id', models.IntegerField(primary_key=True, serialize=False)),
                ('grouplevel', models.CharField(blank=True, db_column='groupLevel', max_length=45, null=True)),
                ('groupname', models.CharField(blank=True, db_column='groupName', max_length=45, null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'TaxGroups',
            },
        ),
        migrations.CreateModel(
            name='Textlabels',
            fields=[
                ('textlabel_id', models.IntegerField(db_column='TextLabel_id', primary_key=True, serialize=False)),
                ('labeltext', models.CharField(blank=True, db_column='LabelText', max_length=500, null=True)),
                ('x', models.FloatField(blank=True, db_column='X', null=True)),
                ('y', models.FloatField(blank=True, db_column='Y', null=True)),
                ('font', models.CharField(blank=True, db_column='Font', max_length=100, null=True)),
                ('font_size', models.FloatField(blank=True, db_column='Font_Size', null=True)),
                ('fill', models.CharField(blank=True, db_column='Fill', max_length=50, null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'TextLabels',
            },
        ),
        migrations.CreateModel(
            name='Threedstructures',
            fields=[
                ('number_3d_structure_id', models.AutoField(db_column='3D_structure_id', primary_key=True, serialize=False)),
                ('structurename', models.CharField(blank=True, db_column='StructureName', max_length=50, null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'ThreeDStructures',
            },
        ),
        migrations.CreateModel(
            name='PolymerMetadata',
            fields=[
                ('polymer', models.ForeignKey(on_delete=django.db.models.deletion.DO_NOTHING, primary_key=True, serialize=False, to='ribovision.PolymerData')),
                ('accession_type', models.CharField(max_length=45)),
                ('polymer_type', models.CharField(max_length=45)),
                ('accession', models.CharField(blank=True, max_length=45, null=True)),
                ('fullseq', models.TextField(blank=True, db_column='Fullseq', null=True)),
            ],
            options={
                'managed': False,
                'db_table': 'Polymer_metadata',
            },
        ),
    ]