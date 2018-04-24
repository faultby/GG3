from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein
import os
import argparse
import sys
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

def createparser():
    parser = argparse.ArgumentParser(
             prog='blast_search.py',
             usage='\n%(prog)s <input_peg_list> <input_query> <input_subjec> <output_folder> [options]',
             description='''All-to-all protein BLAST searches''',
             epilog='(c) Aliaksandr Damenikan, 2018.')

    parser.add_argument('input_peg_list',
                        help='list of pegs to take to balst')
    parser.add_argument('input_query',
                        help='path to input fasta aminoacids for proteins in query genome')
    parser.add_argument('input_subject',
                        help='path to input database fatsa aminoacids.')
    parser.add_argument('output_folder', help='path to output folder.')
    parser.add_argument('-Pq', '--prefix_query',
                        help='prefix before PEG numbers in query database')
    parser.add_argument('-Ps', '--prefix_subject',
                        help='prefix before PEG numbers in subject database')
    parser.add_argument('-Qc', '--query_coverage',
                        default=False,
                        type=float,
                        help='sets threshold for maximal query cover hits to view')
    parser.add_argument('-E', '--evalue',
                        type=float,
                        default=0.001,
                        help='e-value cut-off for BLAST searches, default is 0.001')

    return parser


args = createparser()
enter = args.parse_args()
arguments = sys.argv[1:0]

input_peg_numbers = open(enter.input_peg_list, 'r')
input_query = open(enter.input_query, 'r')
input_subject = open(enter.input_subject, 'r')
output_query = open(enter.output_folder+"/out.faa", 'w')

#creating blast database
make_db_cline = 'makeblastdb -in \' \"%s\" \' -parse_seqids -dbtype prot -out  newdb' % (enter.input_subject)
os.system(make_db_cline)

def peg_number(string, prefix):
    return string.replace(prefix, "").replace("\n", "")

def output_builder(blastrecords,
                   blastrecord,
                   query_pegs_list, 
                   subject_pegs_list,
                   best_coverage_alignment, 
                   best_score_alignment):
    #out = 'Query_PEG,Query_product,Subject_PEG,Subject_product,Query_coverage,Score\n'
    
    query_peg = query_pegs_list[blastrecords.index(blastrecord)][0]
    query_product = query_pegs_list[blastrecords.index(blastrecord)][1]
    print query_product
    if "hypothetical protein" in query_product or \
       "Hypothetical protein" in query_product or \
       "Mobile element protein" in query_product:
        print 'hah'
        main_out = ''
        return main_out
    else:
        pass
    subject_peg = ''
    subject_product = ''
    query_coverage = ''
    subject_coverage = ''
    identities = ''
    score = ''
    if enter.query_coverage != False and best_coverage_alignment[0] != "No record" and \
       best_coverage_alignment[2] > enter.query_coverage and \
       best_score_alignment[2] > enter.query_coverage:
        main_out = ''
        return main_out
    else:
        pass
    if best_coverage_alignment[0] != "No record":
        for peg in subject_pegs_list:
            if peg[0] == best_coverage_alignment[0].title.split(' ')[0]: #peg[0] is subject record id
                #print "WOW"
                subject_peg += '\"1)'+' '+peg[0]+'\n'
                subject_product += '\"1)'+' '+peg[1]+'\n'
                query_coverage += '\"1)'+' '+str(best_coverage_alignment[2])+'\n'
                if best_coverage_alignment[1] == 's':
                    score += '\"1)'+' '+str(best_coverage_alignment[3])+'\n'
                    subject_coverage += '\"1) '+str(blastrecord.query_length*best_coverage_alignment[2]/best_coverage_alignment[0].length)+'\n'
                    identities += '\"1) '+str(float(best_coverage_alignment[0].hsps[0].identities)/float(best_coverage_alignment[0].hsps[0].align_length))+'\n'
                    break
                elif best_coverage_alignment[1] == 'm':
                    max_score = 0
                    max_scoverage = 0
                    max_identity = 0
                    for hsp in best_coverage_alignment[0].hsps:
                        scoverage = blastrecord.query_length*best_coverage_alignment[2]/best_coverage_alignment[0].length
                        if hsp.bits > max_score:
                            max_score = hsp.bits
                        if scoverage > max_scoverage:
                            max_scoverage = scoverage
                        if hsp.identities > max_identity:
                            max_identity = hsp.identities
                            ali_length = float(hsp.align_length)
                    score = '\"1) '+str(max_score)+'\n'
                    subject_coverage += '\"1) '+str(float(max_scoverage)/ali_length)+'\n'
                    identities += '\"1) '+str(max_identity)+'\n'
                    break 
    else:
        subject_peg += '\"1)'+' '+'No record'+'\n'
        subject_product += '\"1)'+' '+'No record'+'\n'
        query_coverage += '\"1)'+' '+'-'+'\n'
        subject_coverage += '\"1)'+' '+'-'+'\n'
        identities += '\"1)'+' '+'-'+'\n'
        score += '\"1)'+' '+'-'+'\n'

    if best_score_alignment[0] != "No record":
        for peg in subject_pegs_list:
            if peg[0] == best_score_alignment[0].title.split(' ')[0]: #peg[0] is subject record id
                subject_peg += '2)'+' '+peg[0]+' \"'
                subject_product += "2)"+' '+peg[1]+' \"'
                query_coverage += '2)'+' '+str(best_score_alignment[2])+' \"'
                subject_coverage += '2) '+str(blastrecord.query_length*best_score_alignment[2]/best_score_alignment[0].length)+' \"'
                score += '2)'+' '+str(best_score_alignment[3])+' \"'
                if best_score_alignment[1] == 's':
                    identities += '2) '+str(float(best_score_alignment[0].hsps[0].identities)/float(best_score_alignment[0].hsps[0].align_length))+' \"'
                elif best_score_alignment[1] == 'm':
                    max_identity = 0
                    for hsp in best_score_alignment[0].hsps:
                        if hsp.identities > max_identity:
                            max_identity = hsp.identities
                            ali_length = float(hsp.align_length)
                    identities += '2) '+str(float(max_identity)/ali_length)+' \"'
                break
    else:
        subject_peg += '2)'+' '+'No record'+' \"'
        subject_product += '2)'+' '+'No record'+' \"'
        query_coverage += '2)'+' '+'-'+' \"'
        subject_coverage += '2)'+' '+'-'+' \"'
        identities += '2)'+' '+'-'+' \"'
        score += '2)'+' '+'-'+' \"'
    
    out = '%s,%s,%s,%s,%s,%s,%s,%s\n' % (query_peg, query_product, subject_peg, subject_product, query_coverage, identities, score, subject_coverage)
    return out
    
    
    

query_peg_list = []
for line in input_peg_numbers.readlines():
    if line.startswith(enter.prefix_query):
        #print 'OK'
        pegs_in_line = line.split(', ')
        for peg in pegs_in_line:
            query_peg_list.append(peg_number(peg, enter.prefix_query))

subject_pegs_description = []
subject_records = SeqIO.parse(enter.input_subject, 'fasta')
for record in subject_records:
    subject_pegs_description.append([record.id, record.description.replace(',', '')])
    #print record.id, record.description.replace(',', '')

records = list(SeqIO.parse(input_query, 'fasta'))
records2blast = []
query_pegs_description = []
for peg in query_peg_list:
    for record in records:
        peg_num = peg_number(record.id, enter.prefix_query)
        if peg_num == peg:
            records2blast.append(record) 
            query_pegs_description.append([peg_num, record.description.replace(',', '')])
        
SeqIO.write(records2blast, output_query, "fasta")
output_query.close()
blastp_cline = NcbiblastpCommandline(query=enter.output_folder+"/out.faa", db="~/newdb", evalue=enter.evalue,
                                      outfmt=5, out=enter.output_folder+"/out.xml")
stdout, stderr = blastp_cline()
result_handle = open(enter.output_folder+"/out.xml")
blast_records = list(NCBIXML.parse(result_handle))
main_out='Query_PEG,Query_product,Subject_PEG,Subject_product,Query_coverage(max),Identities(max),Max_score,Subject_coverage(max)\n'
for record in blast_records:
    best_coverage = float(0)
    bestc_alignment = ["No record", 'n', '-', '-']    # 'n'-no record, 's'- single hsp, 'm'-multiple hsp
    best_score = float(0)
    bests_alignment = ["No record", 'n', '-', '-']
        
    for alignment in record.alignments:
        if len(alignment.hsps) == 1:
            #print "I AM HERE", alignment.hsps[0].align_length, record.query_length
            query_coverage = float(alignment.hsps[0].align_length)/float(record.query_length)
            #print query_coverage
            if query_coverage > best_coverage:
                #print 'OK'
                best_coverage = query_coverage
                bestc_alignment[0], bestc_alignment[1] = alignment, 's'
                bestc_alignment[2], bestc_alignment[3] = best_coverage, alignment.hsps[0].bits
            if alignment.hsps[0].bits > best_score:
                best_score = alignment.hsps[0].bits
                bests_alignment[0], bests_alignment[1] = alignment, 's'
                bests_alignment[2], bests_alignment[3] = query_coverage, best_score 
        elif len(alignment.hsps) > 1:
            #print alignment.title
            query_coverage = 0
            max_score = 0
            for hsp in alignment.hsps:
                if float(hsp.align_length)/float(record.query_length) > query_coverage:
                    query_coverage = float(hsp.align_length)/float(record.query_length) #best query coverage from all hsps
                if hsp.bits > max_score:
                    max_score = hsp.bits
            if query_coverage > best_coverage:
                best_coverage = query_coverage
                bestc_alignment[0], bestc_alignment[1] = alignment, 'm'
                bestc_alignment[2], bestc_alignment[3] = best_coverage, max_score
            if max_score > best_score:
                #print max_score, best_score
                best_score = max_score
                bests_alignment[0], bests_alignment[1] = alignment, 'm'
                bests_alignment[2], bests_alignment[3] = query_coverage, best_score
            
    #print bestc_alignment[0].title, bestc_alignment[0].length 
    main_out += output_builder(blast_records,
                                  record,
                                  query_pegs_description,
                                  subject_pegs_description,
                                  bestc_alignment,
                                  bests_alignment)
main_output_file = open(enter.output_folder+"/main_out.csv", 'w')
main_output_file.writelines(main_out)
main_output_file.close()
input_query.close()
input_subject.close()
         
print "Done"
