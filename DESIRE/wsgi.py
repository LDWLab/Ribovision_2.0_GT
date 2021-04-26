"""
WSGI config for DESIRE project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/2.1/howto/deployment/wsgi/
"""

import os, sys

from django.core.wsgi import get_wsgi_application
path = "/home/Desire-DEV/PVDev/env/lib64/python3.6/site-packages/"
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'DESIRE.settings')
if path not in sys.path: sys.path.append(path)
application = get_wsgi_application()
