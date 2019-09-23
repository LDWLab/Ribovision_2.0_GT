from django import forms
from .models import Taxgroups

class TaxgroupForm(forms.ModelForm):
    class Meta:
        model = Taxgroups
        fields = ('superkingdom', 'phyla', 'alignment')

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fields['phyla'].queryset = Taxgroups.objects.none()
        self.fields['alignment'].queryset = Taxgroups.objects.none()