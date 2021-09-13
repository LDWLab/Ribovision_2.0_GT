#!/usr/bin/env python3
import re, sys, json, getpass, mysql.connector, argparse

from Bio import SeqIO, AlignIO
from subprocess import Popen, PIPE
from os import error, remove, path, mkdir
from urllib.request import urlopen
from io import StringIO
import ssl
ssl._create_default_https_context = ssl._create_unverified_context

from upload_aln import main as uploadMain

def create_and_parse_argument_options(argument_list):
    parser = argparse.ArgumentParser(description='Upload data to alignment related tables in DESIRE from a fasta sequence file.\
                                        \nMake sure upload_accession has updated the tables with Residue data.')
    parser.add_argument('seq_file', help='Path to sequence file', type=str)
    parser.add_argument('source', help='Defines superkingdom source (e.g. abe)', type=str)
    parser.add_argument('-aln_method','--alignment_method', help='Alignment method used (default: PROMALS3D)', type=str, default='PROMALS3D')
    parser.add_argument('-host','--db_host', help='Defines database host (default: 130.207.36.76)', type=str, default='130.207.36.76')
    parser.add_argument('-schema','--db_schema', help='Defines schema to use (default: DESIRE)', type=str, default='DESIRE')
    parser.add_argument('-user_name','--uname', help='Defines user name to use (default: ppenev)', type=str, default='ppenev')
    parser.add_argument('-pw','--password', help='Defines user password to use', type=str)
    parser.add_argument('-aln_id','--alignment_id', help='Defines alignment id to add entries to. If not specified uses the name of the file.', type=int)
    parser.add_argument('-commit','--commit_changes', help='Commit the changes to the DB', action="store_true")
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def getAlnId(cursor, alnName):
    cursor.execute(f"SELECT Aln_id from Alignment WHERE Name = '{alnName}'")
    result = cursor.fetchall()
    if len(result) == 0:
        print (f"No alignment with name {alnName}! Skipping!")
        return False
    if len(result) == 1:
        return result[0][0]
    if len(result) > 1:
        print (f"More than one alignment with name {alnName}! Skipping!")
        return False

def fetchAlnOnline(alnID, source):
    url = f'https://proteovision.chemistry.gatech.edu/ortholog-aln-api/{alnID}/{source}'
    data = urlopen(url)
    tempStr = str()
    outDat = str()
    try:
        for line in data:
            tempStr+=line.decode('UTF-8')
        parsedData = json.loads(tempStr)
        outDat = parsedData["Alignment"]
    except:
        outDat = False
    return outDat

def fixSeqID(seqID, alnName, cursor):
    cursor.execute(f"SELECT strain_id FROM Species WHERE strain = '{seqID}'")
    result = cursor.fetchall()
    return f'{alnName}_{str(result[0][0])}_A'

def deleteAlnData(cursor, alnID):
    cursor.execute(f"DELETE FROM `Aln_Data` WHERE (`aln_id` = '{alnID}')")

def closeAndExitWithErr(errMsg, cursor, cnx):
    cursor.close()
    cnx.close()
    sys.exit(errMsg)

def main(commandline_arguments):
    comm_args = create_and_parse_argument_options(commandline_arguments)
    seq_path = comm_args.seq_file
    ###Use source_string to generalize fetchAlnOnline
    source_string = comm_args.source
    translator = {'b': 2, 'a': 2157, 'e':2759}
    sourceDigits = ''
    for letter in source_string:
        sourceDigits += f"{translator[letter]},"
    sourceDigits = sourceDigits[:-1]
    pw = comm_args.password
    if pw is None:
        pw = getpass.getpass("Password: ")
    cnx = mysql.connector.connect(user=comm_args.uname, password=pw, host=comm_args.db_host, database=comm_args.db_schema)
    cursor = cnx.cursor()

    dirPath = f"{path.dirname(seq_path)}/temp/"
    if not path.isdir(dirPath):
        mkdir(dirPath)
    sequences = list(SeqIO.parse(seq_path, "fasta"))
    for sequence in sequences:
        sequence.id = re.sub('\|.*','',sequence.id)
        sequence.description = ''
    cleanSeqPath = f'{dirPath}clean_{path.split(seq_path)[-1]}'
    with open(cleanSeqPath, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

    alnID = comm_args.alignment_id
    aln_name = path.split(seq_path)[-1]\
        .replace('_txid_tagged_nucl.fas', '')\
        .replace('_txid_tagged.fas', '')\
        .replace('e_new.fa', '')\
        .replace('_new.fas', '')\
        .replace('.fas', '')\
        .replace('.fa', '')
    if not comm_args.alignment_id:
        alnID = getAlnId(cursor, aln_name)
        if alnID == False:
            closeAndExitWithErr('Couldn\'t get alignment ID!', cursor, cnx)
    
    alnString = fetchAlnOnline(alnID, sourceDigits)
    if alnString == False:
        closeAndExitWithErr(f'Couldn\'t get alignment id {alnID} from ProteoVision!', cursor, cnx)

    alnStringFile = StringIO(alnString)
    origAln = AlignIO.read(alnStringFile, "fasta")

    truncName = aln_name
    if len(aln_name) > 4:
        truncName=aln_name[1:]

    for seq in origAln:
        seq.id = fixSeqID(' '.join(seq.id.split('_')[1:]), truncName, cursor)
        seq.description = ''
    
    tempAlnPath = f'{dirPath}temp_{path.split(seq_path)[-1]}'
    fh = open(tempAlnPath, "w")
    fh.write(format(origAln, 'fasta'))
    fh.close()

    newAlnPath = f'{dirPath}out_{path.split(seq_path)[-1]}'
    pipe = Popen(f"mafft --quiet --auto --add {cleanSeqPath} {tempAlnPath} > {newAlnPath}", stdout=PIPE, shell=True)
    output = pipe.communicate()[0]

    newAln = AlignIO.read(newAlnPath, "fasta")

    if (origAln.get_alignment_length() != newAln.get_alignment_length()):
        deleteAlnData(cursor, alnID)
        upAln = newAlnPath
    else:
        truncSeqsAlnPath = f'{dirPath}trunc_{path.split(seq_path)[-1]}'
        outSeqs = list()
        seqDict = SeqIO.to_dict(sequences)
        for entry in newAln:
            if entry.id in seqDict.keys():
                outSeqs.append(entry)
        with open(truncSeqsAlnPath, "w") as output_handle:
            SeqIO.write(outSeqs, output_handle, "fasta")
        upAln = truncSeqsAlnPath

    if comm_args.commit_changes:
        cnx.commit()
    cursor.close()
    cnx.close()
    
    if comm_args.commit_changes:
        try:
            print("Running upload alignment with arguments:\n",upAln,source_string,"-schema",comm_args.db_schema, "-aln_id",f"{alnID}","-commit")
            uploadMain([upAln,source_string,"-schema",comm_args.db_schema,"-pw",pw, "-aln_id",f"{alnID}","-commit"])
        except(error):
            print(f"Uploading alignment failed with error\n{error}")
            sys.exit()

    print("Success!")

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

