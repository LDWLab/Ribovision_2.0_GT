 
from django.urls import include, path

from . import views

app_name = 'alignments'

urlpatterns = [
	path('', views.index, name='index'),
	path('rRNA/<str:align_name>/<int:tax_group>', views.rRNA, name='rRNA'),
	path('rProtein/<str:align_name>/<int:tax_group>', views.rProtein, name='rProtein'),
	path('showTaxonomy', views.buildTaxonomy, name='showTaxonomy'),
	path('entropy/<str:align_name>/<int:tax_group>/<int:taxid>', views.entropy, name='entropy'),
]