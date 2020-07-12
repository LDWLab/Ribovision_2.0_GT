#!/usr/bin/env python3
import csv, sys, getopt, getpass, mysql.connector

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

#uname = input("User name: ")
pw = getpass.getpass("Password: ")
cnx = mysql.connector.connect(user='ppenev', password=pw, host='130.207.36.76', database='SEREB')
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
    cursor.execute("SELECT SEREB.TaxGroups.groupName FROM SEREB.Species_TaxGroup\
        INNER JOIN SEREB.TaxGroups ON SEREB.Species_TaxGroup.taxgroup_id=SEREB.TaxGroups.taxgroup_id\
        INNER JOIN SEREB.Species ON SEREB.Species_TaxGroup.strain_id=SEREB.Species.strain_id\
        WHERE SEREB.TaxGroups.groupLevel = 'superkingdom' AND SEREB.Species.strain_id = '"+ID+"'")
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
    #cursor.execute("SELECT SEREB.Nomenclature.nom_id FROM SEREB.Nomenclature\
    #    INNER JOIN SEREB.Old_name ON SEREB.Nomenclature.nom_id=SEREB.Old_name.nomo_id\
    #    WHERE SEREB.Old_name.old_name = '"+name+"' AND SEREB.Old_name.N_B_Y_H_A = 'BAN' AND SEREB.Nomenclature.occurrence = '"+occur+"'")
    cursor.execute("SELECT SEREB.Nomenclature.nom_id FROM SEREB.Nomenclature\
        WHERE SEREB.Nomenclature.new_name = '"+name+"' AND SEREB.Nomenclature.occurrence = '"+occur+"'")
    result = cursor.fetchall()
    #nom_id=result[0][0]
    try:
        nom_id=result[0][0]
    except:
        raise ValueError ("No result for nom_id "+name+" and occurrence "+occur+" in the MYSQL query")
    return nom_id

def upload_resi(poldata_id, fullseq):
    i = 1
    for resi in fullseq:
        query = "INSERT INTO `SEREB`.`Residues`(`PolData_id`,`resNum`,`unModResName`) VALUES('"+poldata_id+"','"+str(i)+"','"+resi+"')"
        cursor.execute(query)
        #print(query)
        i+=1
    return True

def main():
    csv_list = read_csv(csv_path)
    for entry in csv_list:
        superK = superkingdom_info(entry[0])
        nom_id = check_nomo_id(superK[0], entry[3])
        query = "INSERT INTO `SEREB`.`Polymer_Data`(`GI`,`strain_ID`,`nomgd_id`, `GeneDescription`) \
                        VALUES('"+entry[1]+"','"+str(entry[0])+"','"+str(nom_id)+"','"+entry[4].rstrip()+"')"
        print(query)
        cursor.execute(query)
        lastrow_id = str(cursor.lastrowid)
        query = "INSERT INTO `SEREB`.`Polymer_metadata`(`polymer_id`,`accession_type`,`polymer_type`, `Fullseq`) \
                                            VALUES('"+str(lastrow_id)+"','LDW-prot','protein','"+entry[5]+"')"
        cursor.execute(query)
        #print(query)
        upload_resi(str(lastrow_id), entry[5])
    

if __name__ == "__main__":
    main()

cnx.commit()
cursor.close()
cnx.close()
print("Success!")