from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

# takes in a SeqRecord
def translate(sequence):
    reduced_alphabet = {
        'C' : 'C',
        'F' : 'F', 'W' : 'F', 'Y' : 'F',
        'G' : 'G',
        'P' : 'P',
        'D' : 'D', 'E' : 'D',
        'N' : 'N', 'Q' : 'N', 'S' : 'N', 'T' : 'N',
        'R' : 'K', 'K' : 'K', 'H' : 'K',
        'A' : 'V', 'V' : 'V', 'L' : 'V', 'I' : 'V', 'M' : 'V'
    }

    translated_string = ''

    for aa in sequence:
        if aa in reduced_alphabet:
            translated_string += reduced_alphabet[aa]
        else:
            translated_string += aa
    
    translated = SeqRecord(
        Seq(translated_string),
        id = sequence.id,
        name = sequence.name,
        description = sequence.description) 

    return translated

# translates an alignment file from modern to reduced alphabet
def translate_fasta(input_file):
    records = []
    for seq in SeqIO.parse(input_file, "fasta"):
        records.append(translate(seq))

    # length adjustment
    # longest_length = max(len(s) for s in records)
    # records = [s + ((longest_length - len(s)) * '-') for s in records]
    # records = MultipleSeqAlignment(records)
    # SeqIO.write(records, ra_output, "fasta")
    
    return records

def aa_composition(file_path, reduced = True):
    # reads an alignment file in the reduced alphabet and computes the amino acid composition 
    # for that alignment file
    # reduced parameter - is the alignment file in the reduced alphabet?
    # output: a pandas dataframe with the amino acid composition for each polypeptide in the
    # alignment file
    # the mean of the dataframe is quite different from the values in alignment_dictionary

    # alignment = AlignIO.read(file_path, "fasta")
    if reduced:
        # call the translate_fasta function here
        aa_counts = {'C' : 0, 'F' : 0, 'G' : 0, 'P' : 0, 'D' : 0, 'N' : 0, 'K' : 0, 'V' : 0}
        records = translate_fasta(file_path)
    else:
        aa_counts = {'A' : 0, 'C' : 0, 'D' : 0, 'E' : 0, 'F' : 0, 
        'G' : 0, 'H' : 0, 'I' : 0, 'K' : 0, 'L' : 0, 
        'M' : 0, 'N' : 0, 'P' : 0, 'Q' : 0, 'R' : 0, 
        'S' : 0, 'T' : 0, 'V': 0, 'W': 0, 'Y': 0}
        records = SeqIO.parse(file_path, "fasta")
    length = 0
    aa_dict = {}
    
    for record in records:
        seq = str(record.seq)
        length += (len(seq) - seq.count('-'))
        for aa in aa_counts:
            aa_counts[aa] += seq.count(aa)
            species_dict = {}
            for aa in aa_counts:
                species_dict[aa] = aa_counts[aa] / length
            name = record.name.split('|')[0]
            aa_dict[name] = species_dict

    return aa_dict