from django.shortcuts import render
from ribovision.models import Mastertable
from django.http import HttpResponse
from django.http import JsonResponse
from django.views.generic import ListView

def fetchmasterlist(request):
    response = Mastertable.objects.values('SpeciesName', 'DataSetType', 'StructureName', 'LoadString').filter(Active = True)
    return JsonResponse(list(response), safe = False)

def index(request):
    return render(request, 'ribovision/index.html')