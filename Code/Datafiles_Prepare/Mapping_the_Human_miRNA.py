import BlastUtils
from pathlib import Path
import JsonLog


import pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqIO

from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML

import datetime
import os
import pandas as pd
import numpy as np

class Mapping_the_Human_miRNA(object):

    def __init__(self, input_file, data_dir):
        organism = "human"
        print("Mapping the Human miRNA Interactome by CLASH Reveals Frequent Noncanonical Binding")
        print ("https://www.sciencedirect.com/science/article/pii/S009286741300439X")
        print("#############################################")
        print("Organism: {}".format(organism))
        print("#############################################")

        JsonLog.add_to_json('paper', "Mapping the Human miRNA Interactome by CLASH Reveals Frequent Noncanonical Binding")
        JsonLog.add_to_json('Organism', organism)

        self.organism = organism
        self.input_file = input_file
        self.data_dir = data_dir

        self.Mapping_the_Human_miRNA_Files()
        self.read_paper_data()


    def Mapping_the_Human_miRNA_Files (self):
        self.biomart_file = str(self.data_dir / "biomart_3utr.csv")
        self.final_output = str("Datafiles_Prepare/CSV"/ Path (self.organism + "_Mapping_the_Human_miRNA_Data.csv"))


    def read_paper_data(self):
        self.inter_df = pd.read_csv(self.input_file, sep="\t", skiprows=30)
        JsonLog.add_to_json('Num of samples', self.inter_df.shape[0])

    def filter_non_UTR (self):
        # filtter out all the rows without 3'UTR
        self.inter_df = self.inter_df[self.inter_df['3\'UTR'].notnull()]
        JsonLog.add_to_json('UTR3 samples', self.inter_df.shape[0])

    def filter_CDS (self):
        self.inter_df = self.inter_df[self.inter_df['CDS']!=1]
        JsonLog.add_to_json('No CDS samples', self.inter_df.shape[0])

    def extract_ensg_enst (self):
        def ENSG(name, sep):
            a = name.split(sep)
            return a[0]

        def ENST(name, sep):
            a = name.split(sep)
            return a[1]

        self.inter_df['ensg'] = self.inter_df.apply(lambda row: ENSG(row['mRNA_name'], "_"), axis=1)
        self.inter_df['enst'] = self.inter_df.apply(lambda row: ENST(row['mRNA_name'], "_"), axis=1)

    def BIOMART_MERGE_FULL_MATCH(self):
        print ("Start merging CLASH<-BIOMART")
        biomart_df = pd.read_csv(self.biomart_file)

        self.inter_df['full_mrna_seq'] = np.nan
        self.inter_df['full_mrna_seq_match_start'] = np.nan
        self.inter_df['full_mrna_source'] = np.nan


        for index, row in self.inter_df.iterrows():
            # find candidates with the same ensg
            a = biomart_df['ensg'] == row.ensg
            if sum(a) > 0:
                full_mrna_seq_candidate = biomart_df[a]
                # Choose the candidates that contains the mRNA_seq_extended
                f_result = []
                for candidate_index, candidate_row in full_mrna_seq_candidate.iterrows():
                    try:
                        f = candidate_row.sequence.find(row.mRNA_seq_extended)
                    except Exception:
                        f = -1
                    f_result.append(f)
                    # Insert the candidates to a dedicated dataframe
                f_result = np.array(f_result)
                b = f_result != -1
                if sum(b) > 0:
                    final_candidates = full_mrna_seq_candidate[b]
                    # find the candidate with the longest full mrna seq
                    sequence_len = final_candidates.sequence.str.len()
                    row_with_the_longest_sequence = final_candidates[sequence_len == sequence_len.max()].iloc[0]
                    self.inter_df.loc[index, 'full_mrna_seq'] = row_with_the_longest_sequence.sequence
                    self.inter_df.loc[
                        index, 'full_mrna_seq_match_start'] = row_with_the_longest_sequence.sequence.find(
                        row.mRNA_seq_extended)
                    self.inter_df.loc[index, 'full_mrna_source'] = "BIOMART"

    def filter_no_biomart_match (self):
        # Take only the valid rows: full_mrna_source == "BIOMART"
        valid_rows = self.inter_df['full_mrna_source'] == "BIOMART"
        JsonLog.add_to_json("valid biomart results", sum(valid_rows))
        JsonLog.add_to_json("invalid biomart results", (self.inter_df.shape[0] - sum(valid_rows)))
        self.inter_df = self.inter_df[valid_rows]



    def file_formatting (self):
        in_df = self.inter_df
        in_df['Source'] = "Mapping_the_Human_miRNA"
        in_df['Organism'] = self.organism
        in_df['Creation_time'] = JsonLog.get_creation_time()

        # Choose the necessary columns
        in_df_filter = in_df.filter(
            ['Source', 'Organism', 'microRNA_name', 'miRNA_seq', 'mRNA_seq_extended', 'chimeras_decompressed',
             'mRNA_name', 'mRNA_start', 'mRNA_end_extended', 'full_mrna_seq'], axis=1)

        in_df_filter.rename(columns={'mRNA_end_extended': 'mRNA_end', 'miRNA_seq':'miRNA sequence',
                                     'mRNA_seq_extended': 'target sequence', 'chimeras_decompressed':'number of reads',
                                     'full_mrna_seq':'full_mrna'}, inplace=True)


        # reset the index
        in_df_filter.reset_index(drop=True, inplace=True)

        in_df_filter.to_csv(self.final_output)

    def run (self):
        self.filter_non_UTR()
        self.filter_CDS()
        self.extract_ensg_enst()
        self.BIOMART_MERGE_FULL_MATCH()
        self.filter_no_biomart_match()
        self.file_formatting()


        #
        #
        # #####################################################
        # # build the blast DB
        # #####################################################
        # self.create_blast_db()
        #
        # #####################################################
        # # Run the blast against biomart
        # #####################################################
        # self.blast_mRNA()
        #
        # #####################################################
        # # Add the full mRNA seq and report the status
        # #####################################################
        # self.add_full_mRNA()
        #


def main ():
    interaction_file = str(Path("Raw/1-s2.0-S009286741300439X-mmc1.txt"))
    log_dir = "Datafiles_Prepare/Logs/"

    #####################################################
    # Human
    #####################################################
    p_dir = Path ("Human/Raw")
    JsonLog.set_filename(Path(log_dir) / "Mapping_the_Human_miRNA_human.json")
    JsonLog.add_to_json('file name', interaction_file)
    hu = Mapping_the_Human_miRNA (interaction_file, p_dir)
    hu.run()
  #  hu.file_formatting()
    JsonLog.json_print()




if __name__ == "__main__":
    main()


