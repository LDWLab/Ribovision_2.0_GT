#!/usr/bin/env python3
import re, csv, sys, getopt, getpass, mysql.connector

def usage():
    print (\
    "USAGE:\n./upload_accession.py -c [csv_file_path]-h\n\
    -c: defines path to csv file with txids, accessions, database, protein name, description, and sequence.\tREQUIRED\n\
    -h: prints this\
")

try:
    opts, args = getopt.getopt(sys.argv[1:], 'c:h', ['alignment=', 'help'])
except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
    if opt in ('-h', '--help'):
        usage()
        sys.exit(2)
    elif opt in ('-c', '--alignment'):
        csv_path = arg
    else:
        usage()
        sys.exit(2)

uname = input("User name: ")
pw = getpass.getpass("Password: ")
cnx = mysql.connector.connect(user=uname, password=pw, host='130.207.36.76', database='DESIRE')
cursor = cnx.cursor()

def read_csv(csv_path):
    with open(csv_path, 'r') as csv_file:
        reader = csv.reader(csv_file)
        csv_list = list(reader)
    return csv_list

def superkingdom_info(ID):
    '''
    Gets the superkingdom for a strain ID
    '''
    #print(ID)
    cursor.execute("SELECT TaxGroups.groupName FROM Species_TaxGroup\
        INNER JOIN TaxGroups ON Species_TaxGroup.taxgroup_id=TaxGroups.taxgroup_id\
        INNER JOIN Species ON Species_TaxGroup.strain_id=Species.strain_id\
        WHERE TaxGroups.groupLevel = 'superkingdom' AND Species.strain_id = '"+ID+"'")
    results = cursor.fetchall()
    #print(ID,results)
    try:
        superkingdom=(results[0][0])
    except:
        raise ValueError ("No result for specie "+str(ID)+" in the MYSQL query")
    return superkingdom

def check_nomo_id(occur, name):
    '''
    Gets nom_id for new name and superkingdom
    '''
    #cursor.execute("SELECT Nomenclature.nom_id FROM Nomenclature\
    #    INNER JOIN Old_name ON Nomenclature.nom_id=Old_name.nomo_id\
    #    WHERE Old_name.old_name = '"+name+"' AND Old_name.N_B_Y_H_A = 'BAN' AND Nomenclature.occurrence = '"+occur+"'")
    cursor.execute("SELECT Nomenclature.nom_id FROM Nomenclature\
        WHERE Nomenclature.new_name = '"+name+"' AND Nomenclature.PhylogeneticOccurrence = '"+occur+"'")
    result = cursor.fetchall()
    # print ("result: " + str(result))
    #nom_id=result[0][0]
    try:
        nom_id=result[0][0]
    except:
        raise ValueError ("No result for new_name "+name+" and occurrence "+occur+" in the MYSQL query")
    return nom_id

def upload_resi(poldata_id, fullseq):
    i = 1
    for resi in fullseq:
        query = "INSERT INTO `Residues`(`PolData_id`,`resNum`,`unModResName`) VALUES('"+poldata_id+"','"+str(i)+"','"+resi+"')"
        cursor.execute(query)
        #print(query)
        i+=1
    return True

def main():
    csv_list = read_csv(csv_path)
    for entry in csv_list:
        if re.match(r'^#', entry[0]):
            continue
        superK = superkingdom_info(entry[0])
        nom_id = check_nomo_id(superK[0], entry[3])
        strain_id = str(entry[0])
        gi = entry[1]
        query = "INSERT INTO `Polymer_Data`(`GI`,`strain_ID`,`nomgd_id`, `GeneDescription`, `GI_type`) \
                        VALUES('"+gi+"','"+strain_id+"','"+str(nom_id)+"','"+entry[4].rstrip()+"','"+entry[2]+"')"
        print(query)
        cursor.execute(query)

        lastrow_id = str(cursor.lastrowid)
        query = "INSERT INTO `Polymer_metadata`(`polymer_id`,`accession_type`,`polymer_type`,`encoding_location`,`classification`,`Fullseq`) \
                                            VALUES('"+str(lastrow_id)+"','LDW-prot','protein','"+entry[5]+"','"+entry[6]+"','"+entry[7]+"')"
        cursor.execute(query)

        # pdata_id = cursor.lastrow_id


        query = "INSERT INTO `Species_Polymer`(`strain_id`, `nomgd_id`, `GI`) VALUES("+strain_id+","+str(nom_id)+",'" + gi + "')"
        print(query)
        cursor.execute(query)

        upload_resi(str(lastrow_id), entry[5])
    
    cnx.commit()
    cursor.close()
    cnx.close()
    print("Success!")

