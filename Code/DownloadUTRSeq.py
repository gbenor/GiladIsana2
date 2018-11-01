import pandas as pd
import numpy as np
from RichDuplex import *
from SeedFeatures import *
#from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os
#Blast parameter: https://www.ncbi.nlm.nih.gov/books/NBK279684/


# Generate fasta file with all miRNA
# -----------------------------------
fasta_filename = "Data/Human/Parsed/clash_mrna.fasta"
human_clash_data_utr3 = pd.read_csv("Data/Human/Parsed/human_clash_data_utr3.csv")
sequences=[]
for index, row in human_clash_data_utr3.iterrows():
    mRNA_seq = row.mRNA_seq_extended
    mRNA_name = row.mRNA_name
    record = SeqRecord(Seq(mRNA_seq), description=mRNA_name)
    sequences.append (record)
SeqIO.write(sequences, fasta_filename, "fasta")


# Run BLAST on the mrna seq
# -----------------------------------
blast_output_filname = "Data/Human/Parsed/clash_blast.xml"
E_VALUE_THRESH = 1e-20
cline = NcbiblastnCommandline(query=fasta_filename, db="refseq_mrna", strand="plus", evalue=0.001, out=blast_output_filname, outfmt=5, max_hsps=1)
print cline
os.system(str(cline))


#
#
#
#
# i=0
# if b.alignments : #skip queries with no matches
#     print "QUERY: %s" % b.query[:100]
#     for align in b.alignments:
#         for hsp in align.hsps:
#             if hsp.expect < E_VALUE_THRESH:
#                 i+=1
#                 print "MATCH: %s " % align.title[:1000]
#                 print hsp.expect
#             if i>5 :
#                 raise Exception ("gilad")
#
# from Bio.Blast import NCBIXML
# E_VALUE_THRESH = 1e-20
# i=0
# if b.alignments : #skip queries with no matches
#     print "QUERY: %s" % b.query[:100]
#     for align in b.alignments:
#         for hsp in align.hsps:
#             if hsp.expect < E_VALUE_THRESH:
#                 i+=1
#                 print "MATCH: %s " % align.title[:1000]
#                 print hsp.expect
#             if i>5 :
#                 raise Exception ("gilad")