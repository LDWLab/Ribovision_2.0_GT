 
from django.urls import include, path

from . import views

app_name = 'alignments'

urlpatterns = [
	path('', views.index, name='index'),
	path('rRNA/<str:name>/', views.rRNA, name='rRNA'),
	path('rProtein/<str:align_name>/', views.detail, name='detail'),
	path('', views.TaxgroupListView.as_view(), name='taxgroups_changelist'),
	path('add/', views.TaxgroupCreateView.as_view(), name ='taxgroups_add'),
	path('<int:pk>/', views.TaxgroupUpdateView.as_view(), name = 'taxgroups_change'),
	path('showTaxonomy', views.buildTaxonomy)
]