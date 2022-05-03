import csv, mysql.connector, getpass
with open("../Data_tables/LineLabels.csv",'r') as csv_file:
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
    LineLabel_id = (entry[0])
    X1 = (entry[1])
    Y1 = (entry[2])
    X2 = (entry[3])
    Y2 = (entry[4])
    Fill = (entry[5])
    Stroke = (entry[6])
    StrokeWidth = (entry[7])
    StrokeLineJoin = (entry[8])
    StrokeMiterLimit = (entry[9])
    secondary_structure_id = get_SecStr_pk(entry[10])
    if secondary_structure_id != 0:
        query = "INSERT INTO `DESIRE`.`LineLabels`(`LineLabel_id`,`X1`,`Y1`,`X2`,`Y2`,`Fill`,`Stroke`,`StrokeWidth`,`StrokeLineJoin`,`StrokeMiterLimit`,`secondary_structure_id`) VALUES('"+entry[0]+"','"+entry[1]+"','"+entry[2]+"','"+entry[3]+"','"+entry[4]+"','"+entry[5]+"','"+entry[6]+"','"+entry[7]+"','"+entry[8]+"','"+entry[9]+"','"+str(secondary_structure_id)+"')"
        print(entry[10], secondary_structure_id)
        cursor.execute(query)

cnx.commit()
cursor.close()
cnx.close()