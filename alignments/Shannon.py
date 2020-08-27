# This script will calculate Shannon entropy from a MSA.

# Dependencies:

# Biopython, Matplotlib [optionally], Math

"""
Shannon's entropy equation (latex format):
    H=-\sum_{i=1}^{M} P_i\,log_2\,P_i
    Entropy is a measure of the uncertainty of a probability distribution (p1, ..... , pM)
    https://stepic.org/lesson/Scoring-Motifs-157/step/7?course=Bioinformatics-Algorithms&unit=436
    Where, Pi is the fraction of nuleotide bases of nuleotide base type i,
    and M is the number of nuleotide base types (A, T, G or C)
    H ranges from 0 (only one base/residue in present at that position) to 4.322 (all 20 residues are equally
    represented in that position).
    Typically, positions with H >2.0 are considerered variable, whereas those with H < 2 are consider conserved.
    Highly conserved positions are those with H <1.0 (Litwin and Jores, 1992).
    A minimum number of sequences is however required (~100) for H to describe the diversity of a protein family.
"""
import os
import sys
import math
import warnings
import traceback

__author__ = "Joe R. J. Healey"
__version__ = "1.0.0"
__title__ = "ShannonMSA"
__license__ = "GPLv3"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"


def parseArgs(argument_list):
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            description='Compute per base/residue Shannon entropy of a Multiple Sequence Alignment.')

        parser.add_argument('-a',
                            '--alignment',
                            action='store',
                            required=True,
                            help='The multiple sequence alignment (MSA) in any of the formats supported by Biopython\'s AlignIO.')
        parser.add_argument('-s',
                            '--anchor_sequence',
                            action='store',
                            required=False,
                            help='Sequence id for which to trim down the alignment colmn results.')
        parser.add_argument('-f',
                            '--alnformat',
                            action='store',
                            default='fasta',
                            help='Specify the format of the input MSA to be passed in to AlignIO.')
        parser.add_argument('-t',
                            '--alntype',
                            action='store',
                            choices=['aa', 'nucl'],
                            default='aa',
                            help='Specify the type of the input MSA (default = aa).')
        parser.add_argument('-v',
                            '--verbose',
                            action='count',
                            default=0,
                            help='Verbose behaviour, printing parameters of the script.')
        parser.add_argument('-m',
                            '--runningmean',
                            action='store',
                            type=int,
                            default=0,
                            help='Return the running mean (a.k.a moving average) of the MSAs Shannon Entropy. Makes for slightly smoother plots. Providing the number of points to average over switches this on.')
        parser.add_argument('--return_within',
                            action='store_true',
                            help='Return results within another script.')
    except:
        print ("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()
    commandline_args = parser.parse_args(argument_list)
    return commandline_args

def parsefastastring(fastastring):
    '''Parse a fastastring into AlignIO object'''

    from Bio.Align import MultipleSeqAlignment
    align = MultipleSeqAlignment([])
    for entry in fastastring.split('>'):
        if entry == '':
            continue
        align.add_sequence(entry.split('\\n')[0], entry.split('\\n')[1])
    return align

def species_index_to_aln_index(alignment_obj, species_id):
    '''Given a species id and an alignment object,
    returns a index mapping of sequence id -> alignment id.
    Also returns the numerical index for the anchor sequence
    in the alignment object.
    '''
    #print(species_id)
    #for aln in alignment_obj:
    #    print('___' + aln)
    aln_anchor_index=0
    for aln in alignment_obj:
        if species_id in aln.id:
            aln_anchor_map=dict()
            aln_i = 0
            seq_i = 0
            for letter in aln:
                if letter != '-':
                    seq_i+=1
                    aln_anchor_map[seq_i] = aln_i
                aln_i +=1
            break
        aln_anchor_index+=1
    return aln_anchor_map, aln_anchor_index

def truncate_aln(alignment_obj, index_positions, *args, **kwargs):
    aln_anchor_map = kwargs.get('aln_anchor_map', None)
    truncated_aln = alignment_obj[:,1:1]
    if aln_anchor_map is not None:
        for index in index_positions:
            truncated_aln+=alignment_obj[:,aln_anchor_map[index]:aln_anchor_map[index]+1]
    else:
        for index in index_positions:
            truncated_aln+=alignment_obj[:,index-1:index]
    return truncated_aln

def parseMSA(msa, alnformat, verbose, *args, **kwargs):
    """Parse in the MSA file using Biopython's AlignIO
    It can trim down the alignment given an index_sequence id
    and return the index sequence itself."""
    
    index_sequence_id = kwargs.get('index_sequence_id', None)
    index_sequence = ''
    from Bio import AlignIO

    if alnformat == 'fastastring':
        alignment = parsefastastring(msa)
    else:
        alignment = AlignIO.read(msa, alnformat)

    if index_sequence_id is not None:
        aln_anchor_map, aln_anchor_index = species_index_to_aln_index(alignment, index_sequence_id)
        alignment = truncate_aln(alignment, list(aln_anchor_map.keys()), aln_anchor_map=aln_anchor_map)
        index_sequence = str(alignment[aln_anchor_index].seq)
    # Do a little sanity checking:
    seq_lengths_list = []
    for record in alignment:
       seq_lengths_list.append(len(record))

    seq_lengths = set(seq_lengths_list)

    if verbose > 0: print("Alignment length is:" + str(list(seq_lengths)))

    if len(seq_lengths) != 1:
        sys.stderr.write("Your alignment lengths aren't equal. Check your alignment file.")
        sys.exit(1)

    index = range(1, list(seq_lengths)[0]+1)

    return alignment, list(seq_lengths), index, index_sequence

##################################################################
# Function to calcuate the Shannon's entropy per alignment column
# H=-\sum_{i=1}^{M} P_i\,log_2\,P_i (http://imed.med.ucm.es/Tools/svs_help.html)
# Gaps and N's are included in the calculation
##################################################################

def gap_adjusted_frequency(column_list, all_residues):
    
    #Still doesn't handle ambiguous letters well
    if len(all_residues) >= 20:
        abs_length = 20
        adjsuted_column_list = ['-' if resi=='X' else resi for resi in column_list]
        all_residues = all_residues.replace('X', '')
    else:
        abs_length = 4
        adjsuted_column_list = ['-' if resi=='N' else resi for resi in column_list]
        all_residues.replace('N', '')
    
    #Gap adjustment
    num_gaps = adjsuted_column_list.count('-')
    gap_freq = num_gaps/abs_length

    M   =  len(adjsuted_column_list)
    frequency_list = list()
    # Number of residues in column
    for base in all_residues:
        n_i = adjsuted_column_list.count(base) # Number of residues of type i
        n_i += gap_freq
        P_i = n_i/float(M) # n_i(Number of residues of type i) / M(Number of residues in column)
        frequency_list.append(P_i)
    return frequency_list

def shannon_entropy(list_input, all_residues):
    """Calculate Shannon's Entropy per column of the alignment (H=-\sum_{i=1}^{M} P_i\,log_2\,P_i)"""
    entropy_list = []
    frequency_list = gap_adjusted_frequency(list_input, all_residues)
    for P_i in frequency_list:
        if P_i == 0:
            continue
        entropy_i = P_i*(math.log(P_i,2))
        entropy_list.append(entropy_i)
    sh_entropy = -(sum(entropy_list))
    return sh_entropy

def shannon_entropy_list_msa(alignment, alntype):
    """Calculate Shannon Entropy across the whole MSA"""
    from Bio.SeqUtils import IUPACData
    if alntype == 'aa':
        all_residues = IUPACData.protein_letters
    if alntype == 'nucl':
        all_residues = IUPACData.unambiguous_rna_letters
    shannon_entropy_list = []
    for col_no in range(len(list(alignment[0]))):
        list_input = list(alignment[:, col_no])
        shannon_entropy_list.append(shannon_entropy(list_input, all_residues))

    return shannon_entropy_list

def running_mean(l, N):
    sum = 0
    result = list(0 for x in l)

    for i in range( 0, N ):
        sum = sum + l[i]
        result[i] = sum / (i+1)

    for i in range( N, len(l) ):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N

    return result

def main(commandline_arguments):
    """Compute Shannon Entropy from a provided MSA."""

    # Parse arguments
    args = parseArgs(commandline_arguments)

    # Convert object elements to standard variables for functions
    msa = args.alignment
    alnformat = args.alnformat
    verbose = args.verbose
    return_within = args.return_within
    runningmean = args.runningmean
    anchorseq = args.anchor_sequence.replace(" ", "_")
    alntype = args.alntype

# Start calling functions to do the heavy lifting

    #Trim or no trim
    if args.anchor_sequence:
        alignment, seq_lengths, index, ix_seq = parseMSA(msa, alnformat, verbose, index_sequence_id=anchorseq)

    else:
        alignment, seq_lengths, index, ix_seq = parseMSA(msa, alnformat, verbose)
    sel = shannon_entropy_list_msa(alignment, alntype)

    if runningmean > 0:
        sel = running_mean(sel, runningmean)

    if return_within is True:
        return_list = []
        for c1, c2, c3 in zip(index, sel, ix_seq):
            return_list.append((c1,c2,c3))
        return return_list

    if verbose > 0: print("Index" + '\t' + "Entropy")
    for c1, c2 in zip(index, sel):
        print(str(c1) + '\t' + str(c2))



if __name__ == '__main__':
    main(main(sys.argv[1:]))
