# #!/usr/bin/env python3

# from Bio import Entrez
# from xml.dom.minidom import parse, parseString
# import xml.dom.minidom as minidom
# import subprocess
# import re
# import mysql.connector
# import getpass

# fileRead = './phase-definitions.txt'
# #uname = input("User name: ")
# #pw = getpass.getpass("Password: ")
# #cnx = mysql.connector.connect(user=uname, password=pw, host='130.207.36.75', database='DESIRE')
# #cursor = cnx.cursor()

# def get_id(file):
#     #for info in file:
#     #    info = info.split(',')
#     #    prot_info = info[1].split('_')
#     #    value = info[2].split()
#     #    values = value.split(';')        
#     file_name = 'phase-definitions.txt'
#     open_file = open(file_name, 'r')
#     read_lines = open_file.readlines()
#     open_file.close()
#     excelFile = open('phase_definitions.csv','w')
#     for line in read_lines:
#         line_parts_list = line.split(',')
#         if ';' in line_parts_list[2]:
#             newsplit = line_parts_list[2].split(';')
#             print(newsplit)
#             line_parts_list.remove(line_parts_list[2])
#             line_parts_list.append(newsplit[0])
#             #print(line_parts_list)
#             line_parts_list.append(newsplit[1].replace('\n',''))
#             if len(newsplit)==3:
#                 line_parts_list.append(newsplit[2].replace('\n',''))

#             print('yas'+str(line_parts_list))
#         else:
#             line_parts_list[2] = line_parts_list[2].replace('\n','')
#             print(line_parts_list)
#         nums = line_parts_list[2].split('-')
#         #print(line_parts_list)
#         if len(line_parts_list) == 3:
#             nums2 = line_parts_list[2].split('-')
#             #print('this is nums2'+' '+str(nums2))
#         if len(line_parts_list) == 4:
#             nums3 = line_parts_list

#         if len(line_parts_list)==3:
#             for x in range(int(nums[1])-int(nums[0])+1):    
#                 query = ("SELECT resi_id FROM DESIRE.Residues\
#                 INNER JOIN DESIRE.Polymer_Data ON DESIRE.Residues.PolData_id = DESIRE.Polymer_Data.PData_id\
#                 INNER JOIN DESIRE.Nomenclature ON DESIRE.Polymer_Data.nomgd_id = DESIRE.Nomenclature.nom_id\
#                 WHERE DESIRE.Polymer_Data.strain_id = '+str(prot_info[1])+' AND DESIRE.Nomenclature.new_name = '+str(prot_info[0][:4])+' AND DESIRE.Residues.resNum = '"+str((int(nums[0])+x))+"';")
#                 #print('norm')
#                 #print(str((int(nums[0])+x)))
#                 #cursor.execute(query)
#                 #print(query)
#         elif len(line_parts_list)==4:
#             jlo = []
#             for x in range(int(nums[1])-int(nums[0])+1):
#                 #print(int(nums[0])+x)
#                 jlo.append(int(nums[0])+x)
#             for x in range(int(nums2[1])-int(nums2[0])+1):
#                 jlo.append(int(nums2[0])+(x))
#                 #print(jlo)
#             for x in jlo:
#                 #print('excepton')
#                 #print(str((int(nums[0])+int(x))))
#                 query = ("SELECT resi_id FROM DESIRE.Residues\
#                 INNER JOIN DESIRE.Polymer_Data ON DESIRE.Residues.PolData_id = DESIRE.Polymer_Data.PData_id\
#                 INNER JOIN DESIRE.Nomenclature ON DESIRE.Polymer_Data.nomgd_id = DESIRE.Nomenclature.nom_id\
#                 WHERE DESIRE.Polymer_Data.strain_id = '' AND DESIRE.Nomenclature.new_name = '' AND DESIRE.Residues.resNum = '"+str((int(nums[0])+int(x)))+"';")
#                 #print(query)
#         elif len(line_parts_list) ==5:
#             newn = []
            
#         #nums2 = line_parts_list[3].split('-')
#         #if len(line_parts_list)==3:
#         #    for x in range(int(nums[1])-int(nums[0])):    
#         #        cursor.execute("SELECT resi_id FROM DESIRE.Residues\
#         #        INNER JOIN DESIRE.Polymer_Data ON DESIRE.Residues.PolData_id = DESIRE.Polymer_Data.PData_id\
#         #        INNER JOIN DESIRE.Nomenclature ON DESIRE.Polymer_Data.nomgd_id = DESIRE.Nomenclature.nom_id\
#         #        WHERE DESIRE.Polymer_Data.strain_id = '"+(prot_info[1])+"' AND DESIRE.Nomenclature.new_name = '"+(prot_info[0][:4])+"' AND DESIRE.Residues.resNum = '"+str((int(nums[0])+x))+"';")
#         #elif len(line_parts_list)==4:
#         #    for x in range((int(nums[1])-int(nums[0]))+(int(nums2[1])-int(nums2[0]))):
                
#     ########################results = cursor.fetchall()
#     # print(results)
    
#     #query = 'SELECT Residues.resi_id, Associated_Data.Data_id FROM Residues INNER JOIN Associated_Data ON Residues'
#     #cursor.execute(query)
    
    
    
    
    
#     # cursor.execute("SELECT DESIRE.Associated_Data.Data_Id FROM DESIRE.Associated_Data")
#     # result_AD = cursor.fetchall()
#     # cursor.execute("SELECT DESIRE.Residues.resi_id FROM DESIRE.Residues")
#     # result_res = cursor.fetchall()
#     # for AD_id in result_AD:
#         # for residueP_id in result_res:
#             # query = ("INSERT INTO `DESIRE`.`AD_Residues`(`AD_id`,`residueP_id`) VALUES('"+str(AD_id)+"','"+str(residueP_id)+"')")
#             # print(query)
#                 # # cursor.execute(query)

# get_id('phase_definitions.txt')
# '''
# with open(fileRead) as f:
#     content = f.readlines()
#     get_id(content)
# cursor.close()
# cnx.close()
# '''
                
#!/usr/bin/env python3

from Bio import Entrez
from xml.dom.minidom import parse, parseString
import xml.dom.minidom as minidom
import subprocess
import re
import mysql.connector
import getpass

fileRead = './Helix_defs_LSUb.txt'
uname = input("User name: ")
pw = getpass.getpass("Password: ")
cnx = mysql.connector.connect(user=uname, password=pw, host='130.207.36.76', database='DESIRE')
cursor = cnx.cursor(buffered=True)
#print('got here')
def get_id(file):
    for info in file:
        info = info.split(',')
        print(info[0])
        prot_info = info[1].split('_')
        #print(prot_info)
    file_name = './Helix_defs_LSUb.txt'
    open_file = open(file_name, 'r')
    read_lines = open_file.readlines()
    open_file.close()
    excelFile = open('phase_definitions.csv','w')
    dicto = {}
    for line in read_lines:
        #print(line)
        line_parts_list = line.split(',')
        print(line_parts_list[0])
        if ';' in line_parts_list[2]:
            newsplit = line_parts_list[2].split(';')
            #print(newsplit)
            line_parts_list.remove(line_parts_list[2])
            for x in newsplit:
                line_parts_list.append(x.replace('\n',''))
        else:
            line_parts_list[2] = line_parts_list[2].replace('\n','')
            #print(line_parts_list)
        nums = line_parts_list[2].split('-')
        nums_range = int(nums[1])-int(nums[0])+1
        #print(nums_range)
        #print('this is nums '+str(nums))
        if len(line_parts_list) > 3:
            nums2 = line_parts_list[3].split('-')
            nums2_range = int(nums2[1])-int(nums2[0])+1
            #print('this is nums2 '+str(nums2))
            if len(line_parts_list) == 5:
                nums3 = line_parts_list[4].split('-')
                #print('this is nums3 '+str(nums3))
                nums3_range = int(nums3[1])-int(nums3[0])+1
        
        length = len(line_parts_list)
        print(line_parts_list, length)
        if length > 2:
            for num in range(nums_range):
                query = ("SELECT resi_id FROM DESIRE.Residues\
                INNER JOIN DESIRE.Polymer_Data ON DESIRE.Residues.PolData_id = DESIRE.Polymer_Data.PData_id\
                INNER JOIN DESIRE.Nomenclature ON DESIRE.Polymer_Data.nomgd_id = DESIRE.Nomenclature.nom_id\
                WHERE DESIRE.Polymer_Data.strain_id = '"+str(prot_info[1])+"' AND DESIRE.Nomenclature.new_name = '"+str(prot_info[0][:4])+"' AND DESIRE.Residues.resNum = '"+str(int(nums[0])+num)+"';")
                print(query)
                cursor.execute(query)
                results = cursor.fetchall()
                query1 = ("SELECT DESIRE.Associated_Data.Data_id FROM DESIRE.Associated_Data WHERE DESIRE.Associated_Data.Type = 'Helix' AND DESIRE.Associated_Data.Value = '"+str(line_parts_list[0])+"';")
                print(query1)
                cursor.execute(query1)
                results1 = cursor.fetchall()
                print('result1', results1)
                query2 = ("INSERT INTO `DESIRE`.`AD_Residues`(`AD_id`,`residueP_id`) VALUES('"+str(results1[0][0])+"','"+str(results[0][0])+"')")
                print(query2)
                try:
                    cursor.execute(query2)
                except Exception as e:
                    print(f"Error: {str(e)}")
                
               
                
            if length > 3:
                for num in range(nums2_range):
                    #print('nums', nums2_range, str(int(nums2[0])), str(num))
                    query = ("SELECT resi_id FROM DESIRE.Residues\
                    INNER JOIN DESIRE.Polymer_Data ON DESIRE.Residues.PolData_id = DESIRE.Polymer_Data.PData_id\
                    INNER JOIN DESIRE.Nomenclature ON DESIRE.Polymer_Data.nomgd_id = DESIRE.Nomenclature.nom_id\
                    WHERE DESIRE.Polymer_Data.strain_id = '"+str(prot_info[1])+"' AND DESIRE.Nomenclature.new_name = '"+str(prot_info[0][:4])+"' AND DESIRE.Residues.resNum = '"+str(int(nums2[0])+num))+"';"
                    cursor.execute(query)
                    results = cursor.fetchall()
                    print(query)
                    print('results3', results)
                    
                    query1 = ("SELECT DESIRE.Associated_Data.Data_id FROM DESIRE.Associated_Data WHERE DESIRE.Associated_Data.Type = 'Helix' AND DESIRE.Associated_Data.Value = '"+str(line_parts_list[0])+"';")
                    print(query1)
                    
                    cursor.execute(query1)
                    results1 = cursor.fetchall()
                    print('result1', results1)
                    query2 = ("INSERT INTO `DESIRE`.`AD_Residues`(`AD_id`,`residueP_id`) VALUES('"+str(results1[0][0])+"','"+str(results[0][0])+"')")
                    print(query2)
                    try:
                        cursor.execute(query2)
                    except Exception as e:
                        print(f"Error: {str(e)}")
                    
                    #view_query = f"SELECT * FROM DESIRE.AD_Residues WHERE residueP_id = '{results[0][0]}';"
                    #cursor.execute(view_query)
                    #view_res = cursor.fetchall()
                    #print('view_res', view_res)
                    
                    
                if length > 4:
                    for num in range(nums3_range):
                        print(nums3_range)
                        query = ("SELECT resi_id FROM DESIRE.Residues\
                        INNER JOIN DESIRE.Polymer_Data ON DESIRE.Residues.PolData_id = DESIRE.Polymer_Data.PData_id\
                        INNER JOIN DESIRE.Nomenclature ON DESIRE.Polymer_Data.nomgd_id = DESIRE.Nomenclature.nom_id\
                        WHERE DESIRE.Polymer_Data.strain_id = '"+str(prot_info[1])+"' AND DESIRE.Nomenclature.new_name = '"+str(prot_info[0][:4])+"' AND DESIRE.Residues.resNum = '"+str(int(nums3[0])+num)+"';")
                        cursor.execute(query)
                        results = cursor.fetchall()
                        
                        print(query)
                        print('results4', results)
                        
            print(info[0])
            #value = info[2].split(';')
            #values = [s.rstrip() for s in value]
            #values = [s.rstrip('-') for s in value]
            #res = values[0].split(';')
            #res = res[0].split('-')
            #x=int(res[0])
            #y=int(res[1])
            #for resNum in range(x,y+1):
                #print('resnum', resNum)
                #cursor.execute("SELECT resi_id FROM DESIRE.Residues\
                     #INNER JOIN DESIRE.Polymer_Data ON DESIRE.Residues.PolData_id = DESIRE.Polymer_Data.PData_id\
                 #INNER JOIN DESIRE.Nomenclature ON DESIRE.Polymer_Data.nomgd_id = DESIRE.Nomenclature.nom_id\
                 #WHERE DESIRE.Polymer_Data.strain_id = '"+str(prot_info[1])+"' AND DESIRE.Nomenclature.new_name = '"+str(prot_info[0][:4])+"' AND DESIRE.Residues.resNum = '"+str(resNum)+"';")
            #results = cursor.fetchall()
            #print('results_out', results)
            #try:
                #residue_id = (results[0][0])
            #except:
                #raise ValueError("No result for residue "+str(prot_info[1]), str(prot_info[0][:4]), str(resNum)+"in the MYSQL query!")
            #query = ("SELECT DESIRE.Associated_Data.Data_id FROM DESIRE.Associated_Data WHERE DESIRE.Associated_Data.Type = 'Helix' AND DESIRE.Associated_Data.Value = '"+str(line_parts_list[0])+"';")
            #print(query)
            #query1 = ("INSERT INTO `DESIRE`.`AD_Residues`(`AD_id`,`residueP_id`) VALUES('"+str(line_parts_list[0])+"','"+str(results[0][0])+"')")
            #print(query1)
    #cursor.execute(query)




with open(fileRead) as f:
    content = f.readlines()
    get_id(content)
    
#cnx.commit()
cursor.close()
cnx.close()
            

    