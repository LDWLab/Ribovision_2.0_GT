#!/usr/bin/env python3
import re, sys, getpass, mysql.connector, argparse
from Bio import AlignIO

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description='Upload data to alignment related tables in DESIRE from a fasta.\
                                        \nMake sure upload_accession has updated the tables with Residue data.')
    parser.add_argument('alignment_file', help='Path to alignment file', type=str)
    parser.add_argument('source', help='Defines superkingdom source (e.g. abe)', type=str)
    parser.add_argument('-aln_method','--alignment_method', help='Alignment method used (default: PROMALS3D)', type=str, default='PROMALS3D')
    parser.add_argument('-host','--db_host', help='Defines database host (default: 130.207.36.76)', type=str, default='130.207.36.76')
    parser.add_argument('-schema','--db_schema', help='Defines schema to use (default: DESIRE)', type=str, default='DESIRE')
    parser.add_argument('-user_name','--uname', help='Defines user name to use (default: ppenev)', type=str, default='ppenev')
    parser.add_argument('-pw','--password', help='Defines user password to use', type=str)
    parser.add_argument('-aln_id','--alignment_id', help='Defines alignment id to add entries to. If not specified makes a new alignment entry.', type=int)
    parser.add_argument('-commit','--commit_changes', help='Commit the changes to the DB', action="store_true")
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def read_align(aln_name):
    '''
    Reads the fasta file and gets the sequences.
    '''
    alignments = AlignIO.read(open(aln_name), "fasta")
    return alignments

def superkingdom_info(cursor, ID):
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
        raise ValueError ("No result for species "+str(ID)+" in the MYSQL query!")
    return superkingdom

def check_nomo_id(cursor, occur, name):
    '''
    Gets nom_id for new name and superkingdom
    '''
    occur = occur.capitalize()
    cursor.execute("SELECT Nomenclature.nom_id FROM Nomenclature\
        WHERE Nomenclature.new_name = '"+name+"' AND Nomenclature.PhylogeneticOccurrence = '"+occur+"'")
    result = cursor.fetchall()
    try:
        nom_id=result[0][0]
    #If no result maybe alignment is using BAN nomenclature
    except:
        cursor.execute("SELECT Nomenclature.nom_id FROM Nomenclature\
            INNER JOIN Old_name ON Nomenclature.nom_id=Old_name.nn_fk_id\
            WHERE Old_name.old_name = '"+name+"' AND Old_name.N_B_Y_H_A = 'BAN'\
            AND Nomenclature.PhylogeneticOccurrence = '"+occur+"'")
        result = cursor.fetchall()
        try:
            nom_id=result[0][0]
        except:
            raise ValueError ("No result for name "+name+" and phylogenetic occurrence "+occur+" in the MYSQL query!")
    return nom_id

def check_polymer(cursor, taxid, nomid, gi):
    '''
    Gets polymer id for a given taxid and nomid (LDW-prot requirement)
    '''
    cursor.execute("SELECT Polymer_Data.PData_id FROM Polymer_Data\
                    INNER JOIN Polymer_metadata ON Polymer_Data.PData_id = Polymer_metadata.polymer_id WHERE \
                    Polymer_metadata.accession_type = 'LDW-prot' AND \
                    Polymer_Data.nomgd_id = "+nomid+" AND \
                    Polymer_Data.strain_id = "+taxid+" AND \
                    Polymer_Data.GI = '"+gi+"'")
    result = cursor.fetchall()
    try:
        pol_id=result[0][0]
    except:
        pol_id = 'NOVAL'
        print("No result for nomgd_id "+nomid+" and taxid "+taxid+" and gi "+gi+" in the MYSQL query!")
        #raise ValueError ("No result for nom_id "+nomid+" and taxid "+taxid+" in the MYSQL query!")
    return pol_id

def upaln_getid(cursor, aln_name, source_string, method_name):
    '''
    Uploads alignment name and method, then returns its primary key from the DB.
    '''

    query = "INSERT INTO `Alignment`(`Name`,`Method`,`Source`) VALUES('"+aln_name+"','"+method_name+"','"+source_string+"')"
    print(query)
    cursor.execute(query)
    lastrow_id = str(cursor.lastrowid)
    return lastrow_id

def upload_pol_aln(cursor, pol_id, aln_id):
    '''Populate table Polymer_Alignments'''
    cursor.execute("SELECT PData_id,Aln_id FROM Polymer_Alignments WHERE\
                    PData_id = '"+pol_id+"' AND\
                    Aln_id = '"+aln_id+"'")
    result = cursor.fetchall()
    if len(result) == 0:
        query = "INSERT INTO `Polymer_Alignments`(`PData_id`, `Aln_id`) VALUES('"+pol_id+"','"+aln_id+"')"
        cursor.execute(query)
        return True
    if len(result) == 1:
        return True
    if len(result) > 1:
        raise ValueError("Bad primary key combination in Polymer_alignments! PData_id:"+pol_id+" Aln_id: "+aln_id)

def upload_aln_data(cursor, entry, seq_aln_pos, aln_id, polymer_id):
    '''Uploads aln_data table'''
    #print (entry.seq[seq_aln_pos[1]-1],end='')
    # if str(seq_aln_pos[0]) == '1' and str(polymer_id) == '11090' and str(entry.seq[seq_aln_pos[1]-1]) == 'T':
    #     flag = True
    resi_id = check_resi_id(cursor, str(seq_aln_pos[0]), str(polymer_id), str(entry.seq[seq_aln_pos[1]-1]))
    #print(resi_id)
    cursor.execute("SELECT aln_id,res_id FROM Aln_Data WHERE\
            res_id = '"+str(resi_id)+"' AND\
            aln_id = '"+str(aln_id)+"'")
    result = cursor.fetchall()
    if len(result) == 0:
        query = "INSERT INTO `Aln_Data`(`aln_id`,`res_id`,`aln_pos`) VALUES('"+str(aln_id)+"','"+str(resi_id)+"','"+str(seq_aln_pos[1])+"')"
        #print(query)
        cursor.execute(query)
        return True
    if len(result) == 1:
        return False

def increment_by_range(sequence, increment_range):
    '''Cretes the mapping for a set of sequence ranges'''
    aln_pos = 0
    seq_pos = 0
    seq_order = list()
    resi_mapping=[]
    for one_range in increment_range:
        for ix in range(one_range[0], one_range[1]+1):
            seq_order.append(ix)
    for seq in sequence:
        if seq == '-':
            aln_pos+=1
            continue
        resi_mapping.append((seq_order[seq_pos],aln_pos+1))
        seq_pos+=1
        aln_pos+=1
    return resi_mapping

def create_aln_mapping(entry):
    '''
    Creates list of pairs mapping sequence to alignment positions.
    '''
    sequence = entry.seq
    increment_ranges = list()
    if len(entry.id.split("_")) == 3:
        increment_ranges.append((1,len(str(sequence).replace("-",""))))
    if len(entry.id.split("_")) > 3:
        string_ranges = re.sub(r'\/.*','', entry.id.split("_")[3])
        for one_range in string_ranges.split(","):
            increment_ranges.append((int(one_range.split('-')[0]), int(one_range.split('-')[1])))

    resi_mapping = increment_by_range(sequence, increment_ranges)

    return resi_mapping

def check_resi_id(cursor, seqnum, polid, resname):
    '''
    Gets residue id for a given polymer and a sequence number.
    Also checks if the data (resname) is correct.
    '''
    query = "SELECT Residues.resi_id, Residues.unModResName FROM Residues WHERE \
            Residues.PolData_id = "+polid+" AND \
            Residues.resNum = "+seqnum
    # print (query)
    cursor.execute(query)
    result = cursor.fetchall()
    try:
        res_id=result[0][0]
    except:
        raise ValueError ("No result for pol_id "+polid+" and sequence number "+seqnum+" in the MYSQL query!")
    if result[0][1] == resname:
        return res_id
    else:
        raise ValueError ("Wrong residue name for polymer id "+polid+" and residue "+seqnum+". Is there something wrong with the numbering?")

def fix_old_taxid(taxid):
    if taxid == '1148':
        return '1080228'
    if taxid == '269483':
       return '482957'
    if taxid == '83333':
       return '511145'
    if taxid == '45157':
        return '280699'
    if taxid == '1936271':
       return '1841599'
    if taxid == '999953':
       return '185431'
    if taxid == '5660':
       return '420245'
    if taxid == '35128':
       return '296543'
    if taxid == '44689':
        return '352472'
    return taxid

def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)
    aln_path = comm_args.alignment_file
    source_string = comm_args.source
    pw = comm_args.password
    if pw is None:
        pw = getpass.getpass("Password: ")

    cnx = mysql.connector.connect(user=comm_args.uname, password=pw, host=comm_args.db_host, database=comm_args.db_schema)
    cursor = cnx.cursor()

    alns = read_align(aln_path)
    aln_name = aln_path.split("/")[-1]\
        .replace('.taxid_tagged.fa', '')\
        .replace('_txid_tagged_nucl.fas', '')\
        .replace('_txid_tagged.fas', '')\
        .replace('_new.fas', '')\
        .replace('.fas', '')\
        .replace('.fa', '')    #Fix that with re
    if comm_args.alignment_id:
        aln_id = comm_args.alignment_id
    else:
        aln_id = upaln_getid(cursor, aln_name, source_string, comm_args.alignment_method)
    for entry in alns:
        entry_id_split = entry.id.split('_')
        taxid = fix_old_taxid(entry_id_split[1])
        superK = superkingdom_info(cursor, taxid)
        entry_id_split_2 = entry_id_split[2]
        gi = entry_id_split_2[entry_id_split_2.index('|') + 1:]
        # nom_id = check_nomo_id(cursor, superK[0], entry.id.split('_')[0][:4])
        print ('entry.id: ' + str(entry.id))
        nom_id = check_nomo_id(cursor, superK[0], entry.id.split('_')[0])
        polymer_id = check_polymer(cursor, str(taxid),str(nom_id), gi)
        if polymer_id == 'NOVAL':
            continue
        upload_pol_aln(cursor, str(polymer_id), str(aln_id))
        mapped_resis = create_aln_mapping(entry)
        for seq_aln_pos in mapped_resis:
            upload_aln_data(cursor, entry, seq_aln_pos, aln_id, polymer_id)
        #print()
    if comm_args.commit_changes:
        cnx.commit()
    cursor.close()
    cnx.close()
    print("Success!")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

