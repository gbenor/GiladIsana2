import JsonLog

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from pathlib import Path
from functools import lru_cache
import datetime
import os
import pandas as pd
import numpy as np
from Bio.SeqFeature import SeqFeature, FeatureLocation
import mirBaseUtils as MBU
from Unambiguous_Identification import Unambiguous_Identification

class Pairing_Beyond_Seed(object):
################################################################################
# Strategy:
# Tranform the input file into the form of Unambiguous_Identification.
# Ihe steps are:
# * Extract the target from the genome
# * Extract the miRNA from mature mrna
# When finishing, we can run Unambiguous_Identification and get the desire datafile.


    def __init__(self, input_file, data_dir):
        organism = "elegans"
        print ("Pairing beyond the Seed Supports MicroRNA Targeting Specificity")
        print ("https://www.sciencedirect.com/science/article/pii/S1097276516305214#mmc3")
        print("#############################################")
        print("Organism: {}".format(organism))
        print("#############################################")

        self.organism = organism
        self.input_file = input_file
        self.data_dir = data_dir

        self.Pairing_Beyond_Seed_Files()
        self.read_paper_data()

    def Pairing_Beyond_Seed_Files (self):
        self.genome_file = "Celegans/Raw/ce10/Caenorhabditis_elegans/UCSC/ce10/Sequence/WholeGenomeFasta/genome.fa"
        self.elegans_with_target = "Celegans/Raw/elegans_with_target.csv"
        self.mirbase_fasta = "all_mature_miRNA.fa"
        self.pairing_beyond_to_unambiguous_identification = "Celegans/Raw/pairing_beyond_to_unambiguous_identification.csv"
        self.unambiguous_output = str(self.data_dir / Path (self.organism + "_Pairing_Beyond_Seed_unambiguous_output.csv"))
        self.final_output = str(self.data_dir / Path (self.organism + "_Pairing_Beyond_Seed_Data.csv"))


    def read_paper_data(self):
        self.sheets = {"elegans" : "Table S2"}

        current_sheet = self.sheets[self.organism]

        # Read the interaction file
        inter_df = pd.read_excel(self.input_file, sheet_name=current_sheet, header=None)
        assert "Table S2".lower() in inter_df.iloc[0, 0].lower(), "Read the wrong sheet. no {} in the first cell".format("Table S2")
        self.inter_df = pd.read_excel(self.input_file, sheet_name=current_sheet, skiprows=3)


    def read_genome (self):
        gen = SeqIO.parse(self.genome_file, "fasta")
        self.genome_list = list (gen)

    @lru_cache(maxsize=None)
    def get_genome_by_id(self, id):
        return [g for g in self.genome_list if g.id == id][0]

    def extract_target (self):
        inter_df = self.inter_df
        #fasta_out, utr_max_len=10000):
        for index, row in inter_df.iterrows():
            chr_id = row['chr']
            chr = self.get_genome_by_id(chr_id)
            strand = row['strand']
            assert strand=='+' or strand=='-', "strand value incorrect {}".format(strand)
            strand = 1 if strand=='+' else -1
            start = row['start']
            stop = row['stop']
            target_feature = SeqFeature(FeatureLocation(start, stop, strand=strand))
            target = target_feature.location.extract(chr)
            inter_df.loc[index, 'target'] = str(target.seq)
        self.inter_df = inter_df.copy()

    def prepare_for_Unambiguous_Identification(self):
        self.inter_df.rename(
            columns={'mir_ID' : 'miRNA ID',
                     'miRNA_seq' : 'miRNA sequence',
                     'target' : 'target sequence'}	, inplace=True)
        self.inter_df['number of reads'] = 1
        self.inter_df.to_csv(self.pairing_beyond_to_unambiguous_identification)


    def final_file_formatting (self, in_df):

        JsonLog.add_to_json('paper', "Pairing beyond the Seed Supports MicroRNA Targeting Specificity")
        JsonLog.add_to_json('Organism', self.organism)

        in_df['Source'] = "Pairing_Beyond_Seed"
        in_df['Organism'] = self.organism

        # Remove miRNA with stars
        rows_without_stars = in_df['microRNA_name'].apply(lambda x: x.find('star') == -1)
        JsonLog.add_to_json("valid miRNA", sum(rows_without_stars))
        JsonLog.add_to_json("invalid miRNA (star)", (in_df.shape[0] - sum(rows_without_stars)))
        in_df = in_df[rows_without_stars]

        # reset the index
        in_df.reset_index(drop=True, inplace=True)

        # save to file
        in_df.to_csv(self.final_output)

    def run(self):
            #####################################################
            # Add the target
            #####################################################
            self.read_genome()
            self.extract_target()

            #####################################################
            # Add the miRNA
            #####################################################
            self.inter_df['mir_ID'] = self.inter_df['name'].apply(lambda x: x[:x.find("_")])
            self.inter_df = MBU.insert_mirna(self.inter_df, "Human/Raw/mature.fa", "mir_ID")
            self.inter_df.to_csv(self.elegans_with_target)


            #####################################################
            # Run Unambiguous_Identification.
            #####################################################
            self.prepare_for_Unambiguous_Identification()
            p_dir = Path("Celegans/Raw")
            ce = Unambiguous_Identification("elegans", self.pairing_beyond_to_unambiguous_identification, p_dir, True)
            ce.run()
            df = ce.file_formatting()
            self.final_file_formatting(df)


def main ():

    Pasquinelli_data = "Raw/1-s2.0-S1097276516305214-mmc3.xlsx"
    log_dir = "Logs/Datafiles_Prepare/"

    #####################################################
    # Celegans
    #####################################################
    JsonLog.set_filename(Path(log_dir) / "Pairing_Beyond_Seed_Celegans.json")
    ce = Pairing_Beyond_Seed (Pasquinelli_data, "Celegans")
    ce.run()


if __name__ == "__main__":
    main()


