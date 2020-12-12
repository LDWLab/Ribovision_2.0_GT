from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# from Bio.Align import MultipleSeqAlignment

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

# trims fasta by a list of indices
def trim_fasta_by_index(input_file, indices):
    align = AlignIO.read(input_file, "fasta")
    trimmed_align = align[:,indices[0]:indices[0]+1] # initialize align object
    for i in indices[1:]:
        trimmed_align += align[:,i:i+1]
    return trimmed_align

# translates an alignment file from modern to reduced alphabet
def translate_fasta(input_file, indices = None):
    records = []
    if indices is not None:
        for seq in trim_fasta_by_index(input_file, indices):
            records.append(translate(seq))
    else:
        for seq in AlignIO.read(input_file, "fasta"):
            records.append(translate(seq))
        

    # length adjustment
    # longest_length = max(len(s) for s in records)
    # records = [s + ((longest_length - len(s)) * '-') for s in records]
    # records = MultipleSeqAlignment(records)
    
    return records

def aa_composition(file_path, reduced = True, indices = None):
    # reads an alignment file in the reduced alphabet and computes the amino acid composition 
    # for that alignment file
    # reduced parameter - is the alignment file in the reduced alphabet?
    # output: a pandas dataframe with the amino acid composition for each polypeptide in the
    # alignment file
    # the mean of the dataframe is quite different from the values in alignment_dictionary

    if reduced:
        # call the translate_fasta function here
        aa_counts = {'C' : 0, 'F' : 0, 'G' : 0, 'P' : 0, 'D' : 0, 'N' : 0, 'K' : 0, 'V' : 0}
        records = translate_fasta(file_path)

    else: # not reduced
        aa_counts = {'A' : 0, 'C' : 0, 'D' : 0, 'E' : 0, 'F' : 0, 
        'G' : 0, 'H' : 0, 'I' : 0, 'K' : 0, 'L' : 0, 
        'M' : 0, 'N' : 0, 'P' : 0, 'Q' : 0, 'R' : 0, 
        'S' : 0, 'T' : 0, 'V': 0, 'W': 0, 'Y': 0}
        if indices is not None:
            records = trim_fasta_by_index(file_path, indices)
        else:
            records = AlignIO.read(file_path, "fasta")

    length = 0
    aa_dict = {}
    
    for record in records:
        seq = str(record.seq)
        length += (len(seq) - seq.count('-'))
        for aa in aa_counts:
            aa_counts[aa] += seq.count(aa)
            name = record.name.split('|')[0]
            species_dict = {'name' : name}
            for aa in aa_counts:
                if length > 0:
                    species_dict[aa] = aa_counts[aa] / length
                else:
                    species_dict[aa] = 0
            aa_dict[name] = species_dict

    return aa_dict