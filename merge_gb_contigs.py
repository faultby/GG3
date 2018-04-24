from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
import argparse
import sys

def createparser():
    parser = argparse.ArgumentParser(
             prog='merge_gb_contigs.py',
             usage='\n%(prog)s <report_file> <input_file> <output_file> [options]',
             description='''This script artificially merges contigs in the assembly''',
             epilog='(c) Aliaksandr Damienikan, 2018.')

    parser.add_argument('input_file',
                        help='path to input Genbank file.')
    parser.add_argument('output_file', help='path to output Genbank file.')
    parser.add_argument('-N', '--nucleotides',
                        default = 0,
                        type=int,
                        help='number of N between merged contigs')

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

records = SeqIO.parse(input_handle, 'genbank')
merged_seq = Seq('', generic_dna)
merged_record = SeqRecord(merged_seq, features=[])
location_offset = 0
for record in records:
    merged_record.seq += record.seq + 'N'*enter.nucleotides
    for feature in record.features:
        my_feature = SeqFeature(
                         location=FeatureLocation(feature.location.start+location_offset, feature.location.end+location_offset),
                         type=feature.type,
                         strand=feature.strand,
                         qualifiers=feature.qualifiers)
        merged_record.features.append(my_feature)
    location_offset += len(record.seq) + enter.nucleotides
print merged_record
        
SeqIO.write(merged_record, output_handle, 'genbank')
input_handle.close()
output_handle.close()
print 'Looks fine!'
