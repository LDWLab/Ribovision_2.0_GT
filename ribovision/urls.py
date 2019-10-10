from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('api/RiboVision/v1.0/fetchMasterList', views.fetchmasterlist),
]