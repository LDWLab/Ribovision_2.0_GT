from django.urls import path

from . import views

urlpatterns = [
    path('', views.index, name='index'),
    path('api/RiboVision/v1.0/fetchMasterList', views.fetchmasterlist),
    path('api/RiboVision/v1.0/fetchResidues', views.fetchresidues),
    path('api/RiboVision/v1.0/speciesTable', views.speciestable),
    path('api/RiboVision/v1.0/fetchStructureName', views.fetchstructurename),
    path('api/RiboVision/v1.0/textLabels', views.textlabels),
    path('api/RiboVision/v1.0/lineLabels', views.linelabels),
    path('api/RiboVision/v1.0/fetchInteractionsMenu', views.fetchinteractionsmenu),
    path('api/RiboVision/v1.0/structdatamenu', views.structdatamenu),   #Skip for now
    #path('api/RiboVision/v1.0/fullTable', views.fulltable),
    #path('api/RiboVision/v1.0/basePairs', views.basepairs),
    #path('api/RiboVision/v1.0/fetchStructData', views.fetchstructdata),
    #path('api/RiboVision/v1.0/fetchInteractions', views.fetchinteractions),
    #path('api/RiboVision/v1.0/savepml', views.savepml),
    #path('api/RiboVision/v1.0/save1D', views.save1D),
    #path('api/RiboVision/v1.0/save2D', views.save2D),
]