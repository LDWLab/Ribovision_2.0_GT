from alignments import handleStructureRequests, mapStrucSeqToAln
from django.urls import include, path
from django.views.generic.base import RedirectView
from django.contrib.staticfiles.storage import staticfiles_storage
from . import views

app_name = 'alignments'

urlpatterns = [
    path('', views.index, name='index'),
    path('favicon.ico', RedirectView.as_view(url=staticfiles_storage.url('img/favicon128.png'))),
    path('showTaxonomy', views.buildTaxonomy, name='showTaxonomy'),
    path('flush-session', views.flushSession, name='flushSession'),
    path('showStrucTaxonomy', views.buildFoldTaxonomy, name='showStrucTaxonomy'),
    path('showTaxonomy-api/<int:current_tax>', views.api_showTaxonomy, name='api_showTaxonomy'),
    #path('entropy/<str:align_name>/<int:tax_group>/<str:anchor_structure>', views.entropy, name='entropy'),
    path('orthologs/twincons/<str:anchor_structure>/<str:chain>/<str:align_name>/<int:tax_group1>/<int:tax_group2>/<int:minIndex>/<int:maxIndex>', views.twincons_handler, name='twincons'),
    path('orthologs/twincons/<str:anchor_structure>/<str:chain>/<str:align_name>/<int:tax_group1>/<int:tax_group2>', views.twincons_handler, name='twincons'),
    path('twincons/<str:anchor_structure>/<str:chain>/<str:align_name>/<int:tax_group1>/<int:tax_group2>', views.twincons_handler, name='twincons'),
    path('entropy-api/<str:align_name>/<int:tax_group>/<str:anchor_structure>', views.api_entropy, name='api_entropy'),
    path('twc-api/<str:align_name>/<int:tax_group1>/<int:tax_group2>/<str:anchor_structure>', views.api_twc, name='api_twc'),
    path('twc-api/', views.api_twc_parameterless, name='api_twc'),
    path('mapSeqAln/', mapStrucSeqToAln.make_map_from_alnix_to_sequenceix_new, name='mapping_aln_to_seq'),
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
    path('custom-aln-data', views.handle_custom_upload_alignment, name='custom_aln_data_handler'),
    path('upload-custom-csv/', views.upload_custom_data, name='upload_custom_data'),
    path('propensity-data/<int:aln_id>/<int:tax_group>', views.propensity_data, name = 'propensity_data'),
    path('propensity-data/<int:aln_id>/<str:tax_group>', views.propensity_data, name = 'propensity_data'),
    path('propensity-data-custom/', views.propensity_data_custom, name = 'propensity_data_custom'),
    path('propensities/<str:align_name>/<int:tax_group>', views.propensities, name = 'propensities'),
    path('propensities/<str:align_name>/<str:tax_group>', views.propensities, name = 'propensities'),
    path('custom-struc-data/<str:strucID>', handleStructureRequests.handleCustomUploadStructure, name = 'custom_structure'),
    path('authEcodQuery', views.ecodPassThroughQuery, name = 'ecodQuery'),
]