 
from django.urls import include, path

from . import views

app_name = 'alignments'

urlpatterns = [
	path('', views.index, name='index'),
	path('rRNA/<str:name>/', views.rRNA, name='rRNA'),
	path('rProtein/<str:align_name>/', views.detail, name='detail'),
	path('showTaxonomy', views.buildTaxonomy, name='showTaxonomy'),
]