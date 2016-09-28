from Bio.Alphabet import IUPAC
import os
import os.path as op
from Bio import SeqIO
from Bio import Alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ssbio import utils
date = utils.Date()

# TODO: what if i want to write one fasta file for multiple sequences?

def write_fasta_file(seq_str, ident, description='{}-ssbioSeq'.format(date.short_date),
                     extension='faa', outpath=None, overwrite=False, ignore_alphabet=False):
    '''
    This writes a fasta file for a single sequence string.
    It also checks if the file exists already and returns the filename.
    You can overwrite existing fasta files by setting overwrite=True

    Input:  seq_str (str) - amino acid string
            ident (str) - ID of the sequence
            outpath (str) - desired PATH of file output
            overwrite (bool) - if you want to overwrite existing files
            ignore_alphabet (boolean): OPTIONAL check the alphabet to see if it contains valid amino acids.
    Output: Filename of fasta file
    '''
    seq = load_seq_string(seq_str, ident=ident, desc=description, ignore_alphabet=ignore_alphabet)
    outfile = "{}.{}".format(ident, extension)

    if outpath:
        if not op.exists(outpath):
            os.mkdir(outpath)
        outfile = op.join(outpath, outfile)

    if op.isfile(outfile) and not overwrite:
        return outfile
    else:
        SeqIO.write(seq, outfile, "fasta")
        return outfile


def load_seq_string(seq_string, ident, name='', desc='', ignore_alphabet=False):
    """Load an amino acid sequence string.

    Args:
        seq_string (str): A protein's amino acid sequence
        ident (str): Database identifier (ie. UniProt or PDB ID)
        name (str): OPTIONAL protein name (ie. gene name)
        desc (str): OPTIONAL description of this sequence (ie. catabolizes G6P)
        ignore_alphabet (boolean): OPTIONAL check the alphabet to see if it contains valid amino acids.

    Returns:
        my_seq (SeqRecord): Biopython SeqRecord object

    """
    my_seq = Seq(seq_string, Alphabet.IUPAC.extended_protein)

    # NOTE: using private _verify_alphabet method in Bio.Alphabet module
    if not Alphabet._verify_alphabet(my_seq) and not ignore_alphabet:
        raise ValueError('Sequence contains invalid characters')

    my_seq_record = SeqRecord(my_seq, id=ident, name=name, description=desc)

    return my_seq_record

def load_fasta_file(filename):
    """Load a FASTA file and returns the sequences as a list of SeqRecords

    Args:
        filename (str): Path to the FASTA file to load

    Returns:
        records (list): list of all sequences in the FASTA file as
        Biopython SeqRecord objects
    """
    if os.path.exists(filename):
        handle = open(filename, "r")
        records = list(SeqIO.parse(handle, "fasta"))
        handle.close()
        return records
    else:
        raise IOError('File does not exist.')




if __name__ == '__main__':
    print(write_fasta_file('ALALLALAL', ident='mypdb', outpath='/tmp/', overwrite=True))
