from django.urls import path

from . import views

app_name = 'alignments'

urlpatterns = [
	path('', views.index, name='index'),
	path('rRNA/<str:name>/', views.rRNA, name='rRNA'),
	path('<str:align_name>/', views.detail, name='detail'),
	# path('', views.TaxgroupsListView.as_view(), name='taxgroups_changelist'),
	# path('add/', views.TaxgroupsCreateView.as_view(), name='taxgroups_add'),
	# path('<int:pk>/', views.TaxgroupsUpdateView.as_view(), name='taxgroups_change'),
	# path('ajax/load-phyla/', views.load_phyla, name='ajax_load_phyla')
]