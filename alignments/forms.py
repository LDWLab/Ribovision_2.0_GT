# from django import forms
# from .models import Taxgroups

# class TaxgroupsForm(forms.ModelForm):
#     class Meta:
#         model = Taxgroups
#         fields = ('superkingdom', 'phyla', 'alignment')

#     def __init__(self, *args, **kwargs):
#         super().__init__(*args, **kwargs)
#         self.fields['superkingdom'].queryset = Taxgroups.objects.none()