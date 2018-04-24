from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import argparse
import sys

def createparser():
    parser = argparse.ArgumentParser(
             prog='merge_fasta_contigs.py',
             usage='\n%(prog)s <report_file> <input_file> <output_file> [options]',
             description='''This script artificially merges contigs in the assembly''',
             epilog='(c) Aliaksandr Damienikan, 2018.')

    parser.add_argument('input_file',
                        help='path to input Genbank file.')
    parser.add_argument('output_file', help='path to output Genbank file.')
    parser.add_argument('id', help='identificator for merged record')
    parser.add_argument('-N', '--nucleotide', 
                        type=int,
                        help='specify how many N nucleotides will split each contig')

    return parser


args = createparser()
enter = args.parse_args()
arguments = sys.argv[1:0]

input_handle = open(enter.input_file, 'r')
if enter.input_file == enter.output_file:
    sys.exit('Sorry, but we can\'t edit input file. Plese give another name \
              to output file!')
try:
    output_handle = open(enter.output_file, 'w')
except IOError:
    sys.exit('Open error! Please check your genbank output path!')

records = SeqIO.parse(input_handle, 'fasta')
merged_seq = ''
for record in records:
    print len(str(record.seq))
    merged_seq += str(record.seq)
    merged_seq +='N'*enter.nucleotide
merged_contigs = SeqRecord(Seq(merged_seq), description=enter.id)
SeqIO.write(merged_contigs, output_handle, 'fasta')
input_handle.close()
output_handle.close()
print 'Looks fine!'
