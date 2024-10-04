import logging
from logging.handlers import TimedRotatingFileHandler
import os

class LoggerSetup:
    def __init__(self, log_directory):
        self.log_directory = log_directory
        self.setup_logger()

    def setup_logger(self):
        logger = logging.getLogger("ribovision3-logger")
        logger.setLevel(logging.DEBUG)
        
        # Ensure log directory exists
        os.makedirs(self.log_directory, exist_ok=True)
        
        # Setup TimedRotatingFileHandler
        handler = TimedRotatingFileHandler(
            os.path.join(self.log_directory, "latest.log"),
            when="midnight",
            interval=1
        )
        
        # Custom namer function to rename files after rotation
        def namer(name):
            base_filename, ext, date = name.rsplit(".", 2)
            return f"{base_filename}.{date}.log"
        
        handler.namer = namer
        
        # Set formatter
        formatter = logging.Formatter(
            "%(asctime)s [%(filename)30s:%(lineno)4s - %(funcName)20s()] %(levelname)10s :  %(message)s"
        )
        handler.setFormatter(formatter)
        
        # Add handler to logger
        logger.addHandler(handler)
        
        # Optionally, add console handler
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        logger.addHandler(console_handler)

    def get_logger(self):
        return logging.getLogger("ribovision3-logger")
    

# def pre_setup():
#     logger = LoggerSetup("/home/.logs/ribovision_live")
    
    
    