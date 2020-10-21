"""DESIRE URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.urls import path, include

urlpatterns = [
    path('', include('alignments.urls')),
    path('alignments/', include('alignments.urls')),
    path('ribovision/', include('ribovision.urls')),
    path('upload/', include('ribovision.urls')),
    path('desire-api/', include('desire_api.urls')),
    path('admin/', admin.site.urls),
]

urlpatterns += staticfiles_urlpatterns()
