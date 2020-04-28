 
from django.urls import include, path

from . import views

app_name = 'alignments'

urlpatterns = [
	path('', views.index, name='index'),
	path('orthologs/', views.index_orthologs, name='orthologs'),
	path('orthologs/rRNA/<str:align_name>/<int:tax_group>', views.rRNA, name='rRNA'),
	path('orthologs/rProtein/<str:align_name>/<int:tax_group>', views.rProtein, name='rProtein'),
	path('upload/', views.upload, name='upload'),
	path('rRNA/<str:align_name>/<int:tax_group>', views.rRNA, name='rRNA'),
	path('rProtein/<str:align_name>/<int:tax_group>', views.rProtein, name='rProtein'),
	path('showTaxonomy', views.buildTaxonomy, name='showTaxonomy'),
	path('showTaxonomy-api/<int:parent>', views.api_showTaxonomy, name='api_showTaxonomy'),
	path('entropy/<str:align_name>/<int:tax_group>/<str:anchor_structure>', views.entropy, name='entropy'),
	path('twincons/<str:align_name>/<int:tax_group1>/<int:tax_group2>/<str:anchor_structure>', views.twincons, name='twincons'),
	path('orthologs/twincons/<str:align_name>/<int:tax_group1>/<int:tax_group2>/<str:anchor_structure>', views.twincons, name='twincons'),
	path('entropy-api/<str:align_name>/<int:tax_group>/<str:anchor_structure>', views.api_entropy, name='api_entropy'),
	path('twc-api/<str:align_name>/<int:tax_group1>/<int:tax_group2>/<str:anchor_structure>', views.api_twc, name='api_twc'),
	path('resi-api/<int:resi_id>', views.resi_info),
	path('struc-api/<int:struc_id>', views.struc_info),
	path('orthologs/twc-api/<str:align_name>/<int:tax_group1>/<int:tax_group2>/<str:anchor_structure>', views.api_twc, name='api_twc'),
	path('upload/submit', views.submit, name = 'submit'),
	path('upload/submitAlignment', views.submitAlignment, name='submitAlignment'),
	path('upload/submitAlignmentText', views.submitAlignmentText, name='submitAlignmentText'),
	path('upload/twincons/<str:anchor_structure>', views.twincons_with_upload, name='api_twc_with_upload'),
	path('upload/twc-api/<str:anchor_structure>', views.api_twc_with_upload, name='api_twc_with_upload')
]