import csv, mysql.connector, getpass
with open("../Data_tables/Interactions.csv",'r') as csv_file:
    reader = csv.reader(csv_file)
    csv_list =list(reader)

    uname = input("User name: ")
pw = getpass.getpass("Password: ")
cnx = mysql.connector.connect(user=uname, password=pw, host='130.207.36.75', database='DESIRE')
cursor = cnx.cursor()

def get_resi_j_fk(residue_j):
    cursor.execute("SELECT * FROM DESIRE.Residues\
    INNER JOIN Polymer_Data on Residues.PolData_id\
    where Polymer_Data.PData_id=Residues.PolData_id\
    and where res_num=%'"+str(residue_j)+"'%;")
    results = cursor.fetchall()
    return results[0][0]

def get_resi_i_fk(residue_i):
    cursor.execute("Select * FROM Residues\
    INNER JOIN Polymer_Data on Residues.PolData_id\
    where Polymer_Data.PData_id=Residues.PolData_id\
    and where res_num=%'"+str(residue_i)+"'%;")
    results = cursor.fetchall()
    return results[0][0]

def get_3D_structure_id(structure_name):
    cursor.execute("SELECT 3D_structure_id from 3DStructures\
    where StructureName='"+str(structure_name)+"'")
    results = cursor.fetchall()
    return results[0][0]


for entry in csv_list[1:]:
    interactions_id = entry[0]
    residue_i= get_resi_i_fk(entry[1])
    residue_j= get_resi_j_fk(entry[2])
    bp_type = entry[3]
    bp_group = entry[4]
    _3D_structure_id = get_3D_structure_id(entry[1])
            #query = "INSERT INTO `DESIRE`.`Interactions`(`interactions_id`,`residue_i`,`residue_j`,`bp_type`,`bp_group`,`3D_structure_id`) VALUES('"+entry[0]+"','"+str(resi_id)+"','"+str(resi_id)+"','"+(entry[3])+"','"+entry[4]+"','"+str(_3D_structure_id)+"')"
            #print(query)
            #cursor.execute(query)
cnx.commit()
cursor.close()
cnx.close()