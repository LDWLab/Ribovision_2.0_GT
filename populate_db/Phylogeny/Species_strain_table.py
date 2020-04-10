#species and strain and species ID
'''everything is in taxonomy, not everything is in current species list'''
import csv
import xlrd
myFile = open('./data/Current_species_list_edit.csv','r')
listofLines = myFile.readlines()
aList = []
for string in listofLines:
    stringline = string.split(',')
    aList.append(stringline)

fileTwo = open('Taxonomy.csv','r')
listofLines2 = fileTwo.readlines()
bList = []
for string in listofLines2:
    stringline = string.split(',')
    bList.append(stringline)

speciesStrain = {}
speciesID = {}

speciesList = []
for x in bList:
    if x[2] == 'species':
        speciesList.append(x[1])

speciesUnderscoreList = []

for x in speciesList:
    x = x.replace(' ','_')
    speciesUnderscoreList.append(x)

strainList = []
strainID = []
strainSpecies = []
z = 0
for x in aList:
    if x[1] == 'Bacteria' or x[1] == 'Archaea' or x[1] == 'Eukaryota':
            strainList.append(x[0])
            strainID.append(x[2])
            strainSpecies.append(x[3])

impStrain = []


matchingList = []
matchingStrainID = []
matchingStrain = []

for x in strainSpecies:
    if x in speciesUnderscoreList:
        matchingList.append(x)    
    else:
        impStrain.append(x)
        strainSpecies.remove(x)

for x in aList: #every line in the current species file
    if x[3] in matchingList:
        x[3].replace('_',' ')
        anotherList = []
        anotherList.append(x[3].replace('_',' '))
        speciesStrain[x[0]] = anotherList

for value in speciesStrain.values():
    for x in bList:
        if x[1] in value:
            value.append(int(x[3]))

for x in bList:
    if x[1].replace(' ','_') in matchingList:
        speciesID[x[1]] = [int(x[3])]


with open('./data/species_and_strain.csv',mode = 'w', newline = '') as speciesStrainFile_file:
    speciesStrainFile_writer = csv.writer(speciesStrainFile_file, delimiter = ',', quotechar = '"')
    speciesStrainFile_writer.writerow(['Species','Strain','Strain ID'])
    for key in speciesStrain:
        speciesStrainFile_writer.writerow([speciesStrain[key][0],key,speciesStrain[key][1]])
    