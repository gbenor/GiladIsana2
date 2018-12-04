import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
from ViennaRNADuplex import *

from InteractionRichPresentation import *

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import itertools as it
from itertools import islice


class MirandaDuplex(object):

    def __init__(self, mir, mrna, tmpdir):
        self.run_miranda(mir, mrna, tmpdir)
        self.IRP = self.parse_miranda_output(self.miranda_out_file)




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

        def extract_seq (s):
            return s.split("'")[1].strip()[0:-2]

        fle = open (miranda_file,"r")
        mr = fle.readlines()
        fle.close()

        output=""
        line_tuples = it.tee(mr, 5)
        line_tuples[1].next()
        line_tuples[2].next()
        line_tuples[2].next()
        next(islice(line_tuples[3], 3, 3), None)
        next(islice(line_tuples[4], 4, 4), None)

        found = False
        for l1, l2, l3, l4, l5 in it.izip(*line_tuples):
            if l1.strip().startswith('Query:'):
                query = extract_seq(l1.strip("\n"))
                ref = extract_seq(l3.strip("\n"))
                interaction = l2[l1.find(query):l1.find(query)+len(query)]
                found = True
                break

                # a = [l3.strip("\n"), l2.strip("\n"),l1.strip("\n")]
                # maxlen = len(max(a, key=len))
                # a = [x + '^' * (maxlen - len(x))  for x in a]
                # a = [x[::-1] for x in a]
                # output += "{}\n{}\n{}\n".format(a[0],a[1],a[2])
                #
                # break

        mrna_bulge = ""
        mrna_inter = ""
        mir_inter = ""
        mir_bulge = ""
        if found :
            query = query[::-1].upper()
            ref = ref[::-1].upper()
            interaction = interaction[::-1]

            self.num_of_pairs = 0

            for i in range (len(interaction)) :
                if interaction[i] == " " :
                    mrna_inter += " "
                    mir_inter += " "
                    mrna_bulge += ref[i]
                    mir_bulge += query[i]
                else :
                    self.num_of_pairs+=1
                    mrna_bulge += " "
                    mrna_inter += ref[i]
                    mir_inter += query[i]
                    mir_bulge += " "
            mrna_bulge = mrna_bulge.replace ("-", " ")
            mrna_inter = mrna_inter.replace ("-", " ")
            mir_inter = mir_inter.replace ("-", " ")
            mir_bulge = mir_bulge.replace ("-", " ")

            self.miranda_presentation = ""
            self.miranda_presentation = self.miranda_presentation + "miranda" +"\n"
            self.miranda_presentation = self.miranda_presentation + "--------"  +"\n"
            self.miranda_presentation = self.miranda_presentation + ref +"\n"
            self.miranda_presentation = self.miranda_presentation + interaction +"\n"
            self.miranda_presentation = self.miranda_presentation + query +"\n"
            self.miranda_presentation = self.miranda_presentation + l5.strip() +"\n"
            # mir_for_vienna = query.replace ("-", "")
            # mrna_for_vienna = ref.replace ("-", "")[::-1]
            # print "Vienna on miranda result"
            # print "-------------------------"
            # print ViennaRNADuplex (mir_for_vienna, mrna_for_vienna)
        else:
            self.num_of_pairs = -1 #No hits
            raise NoMirandaHits

        return InteractionRichPresentation (mrna_bulge, mrna_inter, mir_inter, mir_bulge)

    def __str__(self):
        return self.miranda_presentation + "\n" + str(self.IRP)



class NoMirandaHits(Exception):
    pass

