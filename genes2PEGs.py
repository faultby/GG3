from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
from Bio.SeqFeature import SeqFeature
import argparse
import sys

def createparser():
    parser = argparse.ArgumentParser(
             prog='genes2PEGs.py',
             usage='\n%(prog)s <report_file> <input_file> <output_file> [options]',
             description='''This script adds gene names according to RAST PEG numbers''',
             epilog='(c) Aliaksandr Damienikan, 2018.')

    parser.add_argument('input_file',
                        help='path to input Genbank file.')
    parser.add_argument('output_file', help='path to output Genbank file.')
    parser.add_argument('-P', '--prefix',
                        help='prefix before PEG numbers')

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
for record in records:
    for feature in record.features:
        if feature.type == "CDS" and \
           'db_xref' in feature.qualifiers.keys() and \
            any(enter.prefix in line for line in feature.qualifiers['db_xref']):
            peg_number = feature.qualifiers['db_xref'][0].replace(enter.prefix, '')
            feature.qualifiers['gene'] = ['PEG '+ peg_number]
            my_feature = SeqFeature(
                         location=feature.location,
                         type='gene',
                         strand=feature.strand,
                         qualifiers={'gene':['PEG '+ peg_number]})
            record.features.append(my_feature)

SeqIO.write(record, output_handle, 'genbank')
input_handle.close()
output_handle.close()
print 'Looks fine!'
