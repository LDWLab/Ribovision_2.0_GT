import csv, mysql.connector, getpass
with open("../Data_tables/3DStructures.csv",'r') as csv_file:
    reader = csv.reader(csv_file)
    csv_list =list(reader)


uname = input("User name: ")
pw = getpass.getpass("Password: ")
cnx = mysql.connector.connect(user=uname, password=pw, host='130.207.36.75', database='DESIRE')
cursor = cnx.cursor()

for entry in csv_list:
    3D_structure_id = (entry[0])
    StructureName =(entry[1])
        query = "INSERT INTO `DESIRE`.`3DStructures`(`3D_structure_id`,`StructureName`) VALUES('"+entry[0]+"','"+entry[0]+"')"
        print(entry[10], secondary_structure_id)
        #cursor.execute(query)

cnx.commit()
cursor.close()
cnx.close()