import csv, mysql.connector, getpass
with open("../Data_tables/StructuralData2.csv",'r', encoding="utf8") as csv_file:
    reader = csv.reader(csv_file)
    csv_list =list(reader)


uname = input("User name: ")
pw = getpass.getpass("Password: ")
cnx = mysql.connector.connect(user=uname, password=pw, host='130.207.36.75', database='DESIRE')
cursor = cnx.cursor()

def get_SecStr_pk(SecStr_name):
    cursor.execute("SELECT SecStr_id FROM DESIRE.SecondaryStructures\
                WHERE Name = '"+str(SecStr_name)+"'")
    results = cursor.fetchall()
    try:
        return results[0][0]
    except IndexError:
        return 0

for entry in csv_list[1:]:
    map_index = (entry[0])
    Domain_RN = (entry[1])
    Domain_AN = (entry[2])
    Domains_Color = (entry[3])
    Helix_Num = (entry[4])
    Helix_Color = (entry[5])
    secondary_structure_id = get_SecStr_pk(entry[6])
    if secondary_structure_id != 0:
        query = "INSERT INTO `DESIRE`.`StructuralData2`(`map_index`,`Domain_RN`,`Domain_AN`,`Domains_Color`,`Helix_Num`,`Helix_Color`,`secondary_structure_id`) VALUES('"+entry[0]+"','"+entry[1]+"','"+entry[2]+"','"+entry[3]+"','"+entry[4]+"','"+entry[5]+"','"+str(secondary_structure_id)+"')"
        print(query)
        cursor.execute(query)

cnx.commit()
cursor.close()
cnx.close()