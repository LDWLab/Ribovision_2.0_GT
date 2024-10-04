#!/usr/bin/env python
import os
import sys
import logging
import subprocess
from utils import LoggerSetup

if __name__ == '__main__':
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'DESIRE.settings')    
    try:
        from django.core.management import execute_from_command_line
    except ImportError as exc:
        raise ImportError(
            "Couldn't import Django. Are you sure it's installed and "
            "available on your PYTHONPATH environment variable? Did you "
            "forget to activate a virtual environment?"
        ) from exc
    execute_from_command_line(sys.argv)
    
    # res = subprocess.run(["which", "python3"], capture_output=True, text=True)
    res = subprocess.run(["which", "python3"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    
    try:
        LoggerSetup("./.logs")
        logger = logging.getLogger("ribovision3-logger")
        
        if res.stdout:
            logger.info(res.stdout)
        if res.stderr:
            logger.error(res.stderr)
        
    except Exception as e:
        print("Failed to setup logging!")
        print(str(e))
    
    
