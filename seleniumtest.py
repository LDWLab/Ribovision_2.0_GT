from selenium import webdriver  
import time  
from selenium.webdriver.common.keys import Keys 
from selenium.webdriver.chrome.service import Service 
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
import filecmp
import datetime
import os

print("RiboVision UI test started")  
webdriver_service = Service(f"/home/hmccann3/chromedriver/stable/chromedriver")
chrome_options = Options()
chrome_options.add_argument("--headless") # Ensure GUI is off
chrome_options.add_argument("--no-sandbox")
driver = webdriver.Chrome(service=webdriver_service, options=chrome_options)  
driver.maximize_window()  
#navigate to the url  
driver.get("http://127.0.0.1:8000/")  
ele = driver.find_element(By.CLASS_NAME, "vue-treeselect__multi-value")
ele.click()
time.sleep(3)  
button = driver.find_element(By.XPATH, "//*[text()='Bacteria']")
button.click()
print("Taxonomy browser loaded from database")
ele = driver.find_element(By.CLASS_NAME, "alignment_section")
ele.click()
time.sleep(3) 
ele = driver.find_element(By.ID, "select_protein_type")
ele.click()
ele = driver.find_element(By.XPATH, "//*[text()='LSU-rRNA']")
ele.click()
ele = driver.find_element(By.ID, "selectaln")
ele.click()
time.sleep(3) 
ele = driver.find_element(By.XPATH, "//*[text()='5S']")
ele.click()
time.sleep(3) 
ele = driver.find_element(By.CLASS_NAME, "input-group-text")
ele.send_keys("7k00")
ele.send_keys(Keys.ENTER)
time.sleep(10)  
ele = driver.find_element(By.ID, "downloadAlnImageBtn")
ele.click()
ele = driver.find_element(By.XPATH, "//*[text()='Full alignment']")
ele.click()
time.sleep(3) 
now = datetime.datetime.now()
curr_date = str(now.month) + '-' + str(now.day) + '-' + str(now.year)
assert(filecmp.cmp('./PValignment-full-' + curr_date + '.png', './test_files/PValignment-full-3-12-2023.png'))
print("Alignment viewer successfully loaded")
ele = driver.find_element(By.XPATH, "//*[text()='5S ribosomal RNA']")
print("RNA chains successfully loaded")
ele.click()
time.sleep(90)
img = driver.find_element(By.CLASS_NAME, "saveSVG")
img.click()
time.sleep(3)   
driver.close()  
assert(filecmp.cmp('./rv3Topology.svg', './test_files/rv3Topology.svg'))
print("RNA topology viewer successfully loaded")
print("All test cases successfully completed")  
os.remove('./PValignment-full-' + curr_date + '.png')
os.remove('./rv3Topology.svg')

