from django.urls import include, path

from . import views

app_name = 'alignments'

urlpatterns = [
	path('', views.index_orthologs, name='index'),
	path('orthologs/', views.index_orthologs, name='orthologs'),
	path('node_test/', views.index, name='node_test'),
	path('paralogs/', views.paralog_entry_form, name='paralogs'),
	path('orthologs/rRNA/<str:align_name>/<int:tax_group>', views.rRNA, name='rRNA'),
	path('orthologs/rProtein/<str:align_name>/<int:tax_group>', views.rProtein, name='rProtein'),
	path('orthologs/Visualizer/<str:align_name>/<int:tax_group1>/<int:tax_group2>/<str:anchor_structure>', views.visualizer),
	path('rRNA/<str:align_name>/<int:tax_group>', views.rRNA, name='rRNA'),
	path('rProtein/<str:align_name>/<int:tax_group>', views.rProtein, name='rProtein'),
	path('showTaxonomy', views.buildTaxonomy, name='showTaxonomy'),
	path('showStrucTaxonomy', views.buildFoldTaxonomy, name='showStrucTaxonomy'),
	path('showTaxonomy-api/<int:current_tax>', views.api_showTaxonomy, name='api_showTaxonomy'),
	path('entropy/<str:align_name>/<int:tax_group>/<str:anchor_structure>', views.entropy, name='entropy'),
	path('orthologs/twincons/<str:anchor_structure>/<str:chain>/<str:align_name>/<int:tax_group1>/<int:tax_group2>/<int:minIndex>/<int:maxIndex>', views.twincons_handler, name='twincons'),
	path('orthologs/twincons/<str:anchor_structure>/<str:chain>/<str:align_name>/<int:tax_group1>/<int:tax_group2>', views.twincons_handler, name='twincons'),
	path('twincons/<str:anchor_structure>/<str:chain>/<str:align_name>/<int:tax_group1>/<int:tax_group2>', views.twincons_handler, name='twincons'),
	path('entropy-api/<str:align_name>/<int:tax_group>/<str:anchor_structure>', views.api_entropy, name='api_entropy'),
	path('twc-api/<str:align_name>/<int:tax_group1>/<int:tax_group2>/<str:anchor_structure>', views.api_twc, name='api_twc'),
	path('twc-api/', views.api_twc_parameterless, name='api_twc'),
	path('mapSeqAln/', views.make_map_from_alnix_to_sequenceix, name='mapping_aln_to_seq'),
	path('twc-api/<str:align_name>/<int:tax_group1>/<int:tax_group2>', views.api_twc, name='api_twc_no_struc'),
	path('resi-api/<int:resi_id>', views.resi_info),
	path('struc-api/<int:struc_id>', views.struc_info),
	path('fold-api/<int:fold_id>', views.fold_info),
	path('paralog-aln-api/<int:aln_id>', views.para_aln),
	path('ortholog-aln-api/<int:aln_id>/<int:tax_group>', views.simple_fasta),
	path('ortholog-aln-api/<int:aln_id>/<str:tax_group>', views.simple_fasta, name ='ortholog_aln_api'),
	path('orthologs/twc-api/<str:align_name>/<int:tax_group1>/<int:tax_group2>/<str:anchor_structure>', views.api_twc, name='api_twc'),
	path('orthologs/upload/twincons/<str:anchor_structure>/<str:chain>/<int:minIndex>/<int:maxIndex>', views.twincons_handler, name='twc_with_upload'),
	path('orthologs/upload/twincons/<str:anchor_structure>/<str:chain>', views.twincons_handler, name='twc_with_upload'),
	path('upload/twc-api/<str:anchor_structure>', views.api_twc_with_upload, name='api_twc_with_upload'),
	path('upload-custom-csv/<str:anchor_structure>/<str:chain>', views.twincons_handler, name='custom_csv_data_viewer'),
	path('custom-csv-data', views.upload_custom_data_for_mapping, name='custom_csv_data_handler'),
	path('upload-custom-csv/', views.upload_custom_data, name='upload_custom_data'),
    path('propensity-data/<int:aln_id>/<int:tax_group>', views.propensity_data, name = 'propensity_data'),
    path('propensities/<int:aln_id>/<int:tax_group>', views.propensities, name = 'propensities'),
]