from django.shortcuts import render
from django.http import HttpResponse
from django.http import JsonResponse
from django.views.generic import ListView

class FetchMasterList(ListView):
    def __init__(self,**kwargs):
        self.db = kwargs['db']
        self.cnx = self.db.raw_connection()

    def get(self):
        try:
            cur= self.cnx.cursor()
            SQLStatement = 'SELECT SpeciesName, DataSetType, StructureName, LoadString FROM MasterTable WHERE Active=True'
            cur.execute(SQLStatement)
            r = [dict((cur.description[i][0], value) \
               for i, value in enumerate(row)) for row in cur.fetchall()]
            cur.close()
            return JsonResponse(r)  

        except Exception as e:
            return {'error': str(e)}

def index(request):
    return render(request, 'ribovision/index.html')