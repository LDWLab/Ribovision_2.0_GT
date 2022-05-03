import csv, mysql.connector, getpass, sys
from ..Alignments.upload_aln import check_nomo_id, check_polymer, superkingdom_info
with open(sys.argv[1],'r') as csv_file:
    reader = csv.reader(csv_file)
    csv_list =list(reader)

uname = input("User name: ")
pw = getpass.getpass("Password: ")
cnx = mysql.connector.connect(user=uname, password=pw, host='130.207.36.76', database='DESIRE')
cursor = cnx.cursor()

def get_polymer_fk(strain_id,new_name):
    cursor.execute("Select * from Polymer_Data\
    INNER JOIN Nomenclature on Polymer_Data.nomgd_id\
    where Polymer_Data.nomgd_id=Nomenclature.nom_id\
    and strain_id='"+str(strain_id)+"' and new_name='"+str(new_name)+"';")
    results = cursor.fetchall()
    try:
        return results[0][0]
    except IndexError:
        return 'NOVAL'

def get_3D_structure_id(structure_name):
    cursor.execute("SELECT 3D_structure_id from ThreeDStructures\
    where StructureName='"+str(structure_name)+"'")
    results = cursor.fetchall()
    return results[0][0]

def check_chainListID(pk):
    cursor.execute(f"SELECT ChainList_id FROM DESIRE.ChainList WHERE ChainList_id = {pk}")
    results = cursor.fetchall()
    if len(results) == 0:
        return True
    else:
        return False

pdb_to_strain = {
    '1VY4':262724,
    '3J6B':559292,
    '4V6U':186497,
    '4V6W':7227,
    '4V88':559292,
    '4V9D':511145,
    '4V9F':272569
}
for entry in csv_list[1:]:
    ChainList_id = entry[0]
    occur = superkingdom_info(cursor, str(pdb_to_strain[entry[1]]))
    try:
        nom_id = check_nomo_id(cursor, occur[0], entry[2])
        polymer_id = check_polymer(cursor, str(pdb_to_strain[entry[1]]), str(nom_id))
    except ValueError:
        polymer_id = get_polymer_fk(pdb_to_strain[entry[1]],entry[2])
    if polymer_id != 'NOVAL':
        _3D_structure_id = get_3D_structure_id(entry[1])
        if _3D_structure_id != 0:
            ChainName = entry[3]
            if check_chainListID(entry[0]):
                query = "INSERT INTO `DESIRE`.`ChainList`(`ChainList_id`,`3D_structure_id`,`polymer_id`,`ChainName`) VALUES('"+entry[0]+"','"+str(_3D_structure_id)+"','"+str(polymer_id)+"','"+(entry[3])+"')"
                print(query)
                cursor.execute(query)

cnx.commit()
cursor.close()
cnx.close()

#Run from root dir with:
#python3 -m populate_db.Secondary_Structures.Insert_ChainList ./populate_db/Secondary_Structures/DATA/ChainList.csv