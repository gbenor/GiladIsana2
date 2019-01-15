from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

import os
import numpy as np


def create_blast_db(db_fasta, db_title):
    print ("create_blast_db")
    cmd = "makeblastdb -in {fasta} -parse_seqids -dbtype nucl -out {out}".format(fasta=db_fasta,
                                                                                 out=db_title)
    print (cmd)
    os.system(cmd)

def run_blastn(mrna, db_title, blast_output_filname):
    mRNA_seq = mrna
    mRNA_name = "mrna_to_find"
    filename = "mrna_to_find.fasta"
    record = SeqRecord(Seq(mRNA_seq), description=mRNA_name)
    SeqIO.write(record, filename, "fasta")
    #think about add strand

    cline = NcbiblastnCommandline(query=filename, db=db_title, evalue=0.001,
                                  out=blast_output_filname, outfmt=5, max_hsps=1,
                                  max_target_seqs=1)  # it should return the best result
    cmd = str(cline)
    print (cmd)
    os.system(cmd)


def parse_blast(blast_result, query_seq):
    # assuming only the best result stored in blast_result
    E_VALUE_THRESH = 10
    blast_xml_handle = open(blast_result, "r")
    records = NCBIXML.parse(blast_xml_handle)
    for b in records:
        if b.alignments:  # skip queries with no matches
            for align in b.alignments:
                for hsp in align.hsps:
                    if hsp.expect < E_VALUE_THRESH:
                        identities = hsp.identities
                        full_match = identities == len(query_seq)
                        sbjct_start = hsp.sbjct_start
                        sbjct_end = hsp.sbjct_end
                        title = align.title
    try:
        return full_match, title, sbjct_start, sbjct_end, identities
    except UnboundLocalError:
        return np.nan, np.nan, np.nan, np.nan, np.nan
