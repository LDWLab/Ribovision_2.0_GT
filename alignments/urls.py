from django.urls import path

from . import views

app_name = 'alignments'

urlpatterns = [
	path('', views.index, name='index'),
	path('rRNA/<str:name>/', views.rRNA, name='rRNA'),
	path('<str:align_name>/', views.detail, name='detail'),
]