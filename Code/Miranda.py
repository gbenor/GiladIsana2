import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import itertools as it

class Miranda(object):

    def __init__(self, mir, mrna, tmpdir):
        self.run_miranda(mir, mrna, tmpdir)
        self.duplex = self.parse_miranda_output(self.miranda_out_file)




    def run_miranda (self, mir, mrna, tmpdir):
        mrna_fasta_filename = tmpdir + "/tmp_mrna.fasta"
        mirna_fasta_filename = tmpdir + "/tmp_mirna.fasta"
        miranda_out = tmpdir + "/miranda_out.txt"

        mRNA_record = SeqRecord(Seq(mrna), description="mRNA")
        miRNA_record = SeqRecord(Seq(mir), description="miRNA")
        SeqIO.write(mRNA_record, mrna_fasta_filename, "fasta")
        SeqIO.write(miRNA_record, mirna_fasta_filename, "fasta")

        miranda_cmd = "miranda {mir} {mrna} -out {out} -en 10000 -sc 60 ".format(mir=mirna_fasta_filename, mrna=mrna_fasta_filename, out=miranda_out)
        os.system (miranda_cmd)
        self.miranda_out_file = miranda_out

    def parse_miranda_output (self, miranda_file):
        fle = open (miranda_file,"r")
        mr = fle.readlines()
        fle.close()

        output=""
        line_tuples = it.tee(mr, 3)
        line_tuples[1].next()
        line_tuples[2].next()
        line_tuples[2].next()

        for l1, l2, l3 in it.izip(*line_tuples):
            if l1.strip().startswith('Query:'):
                a = [l3.strip("\n"), l2.strip("\n"),l1.strip("\n")]
                maxlen = len(max(a, key=len))
                a = [x + '^' * (maxlen - len(x))  for x in a]
                a = [x[::-1] for x in a]
                output += "{}\n{}\n{}\n".format(a[0],a[1],a[2])

                break
        return output

    def __str__(self):
        return self.duplex



