import csv, mysql.connector, getpass
with open("../Data_tables/TextLabels.csv",'r', encoding="utf8") as csv_file:
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
    TextLabel_id = (entry[0])
    LabelText = (entry[1])
    X = (entry[2])
    Y = (entry[3])
    Font = (entry[4])
    Font_Size = (entry[5])
    Fill = (entry[6])
    secondary_structure_id = get_SecStr_pk(entry[7])
    if secondary_structure_id != 0:
        query = "INSERT INTO `DESIRE`.`TextLabels`(`TextLabel_id`,`LabelText`,`X`,`Y`,`Font`,`Font_Size`,`Fill`,`secondary_structure_id`) VALUES('"+entry[0]+"','"+entry[1]+"','"+entry[2]+"','"+entry[3]+"','"+entry[4]+"','"+entry[5]+"','"+entry[6]+"','"+str(secondary_structure_id)+"')"
        print(query)
        cursor.execute(query)

cnx.commit()
cursor.close()
cnx.close()