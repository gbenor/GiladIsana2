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

class Unambiguous_Identification(object):

    def __init__(self, organism, input_file, data_dir, run_as_service=False):
        self.run_as_service = run_as_service
        if not run_as_service:
            print("Unambiguous Identification of miRNA:Target Site Interactions by Different Types of Ligation Reactions")
            print ("https://www.sciencedirect.com/science/article/pii/S1097276514003566#app3")
            print("#############################################")
            print("Organism: {}".format(organism))
            print("#############################################")

            JsonLog.add_to_json('paper', "Unambiguous Identification of miRNA:Target Site Interactions by Different Types of Ligation Reactions")
            JsonLog.add_to_json('Organism', organism)

        self.organism = organism
        self.input_file = input_file
        self.data_dir = data_dir

        self.Unambiguous_Identification_Files()
        self.read_paper_data()


    def Unambiguous_Identification_Files (self):
        self.blast_db_fasta_biomart = str(self.data_dir / Path (self.organism +"_biomart.fasta"))
        self.blast_db_fasta_biomart_valid = str(self.data_dir / Path (self.organism +"_biomart_valid.fasta"))
        self.blast_db_biomart = str(self.data_dir / "blast_files"  / "blastdb")
        self.blast_tmp = str(self.data_dir / "blast_files" /"tmp" / "tmp_result.xml")
        self.blast_result = str(self.data_dir / Path (self.organism +"_blast.csv"))
        self.blast_with_mRNA = str(self.data_dir / Path (self.organism +"_mRNA.csv"))
        self.final_output = str("Datafiles_Prepare/CSV"/ Path (self.organism + "_Unambiguous_Identification_Data.csv"))


    def read_paper_data(self):
        self.sheets = {"elegans" : "Sheet2",
                       "human" : "Sheet3",
                       "viral" : "Sheet4",
                       "mouse" : "Sheet5"}

        current_sheet = self.sheets[self.organism]

        # Read the interaction file
        if not self.run_as_service:
            inter_df = pd.read_excel(self.input_file, sheet_name=current_sheet, header=None)
            assert self.organism in inter_df.iloc[0, 0].lower(), "Read the wrong sheet. no {} in the first cell".format(self.organism)
            self.inter_df = pd.read_excel(self.input_file, sheet_name=current_sheet, skiprows=3)
        else:
            self.inter_df = pd.read_csv(self.input_file)

        JsonLog.add_to_json('Num of samples', self.inter_df.shape[0])



    def create_blast_db(self):
        fasta_sequences = SeqIO.parse(open(self.blast_db_fasta_biomart), 'fasta')
        valid_seq = []
        c = 0
        for s in fasta_sequences:
            c+=1
            if s.seq != 'Sequenceunavailable':
                valid_seq.append(s)
        SeqIO.write(valid_seq, self.blast_db_fasta_biomart_valid, 'fasta')
        print ("Finish filtering Biomart fasta")
        JsonLog.add_to_json("Total Biomart seq", c)
        JsonLog.add_to_json("Valid Biomart seq", len (valid_seq))

        BlastUtils.create_blast_db(self.blast_db_fasta_biomart_valid, self.blast_db_biomart)


    def blast_mRNA(self):
        blast_db_biomart = self.blast_db_biomart
        output_file = self.blast_result
        blast_tmp = self.blast_tmp
        inter_df = self.inter_df

        for index, row in self.inter_df.iterrows():
            seq_to_blast = row['target sequence']
            BlastUtils.run_blastn(seq_to_blast, blast_db_biomart, blast_tmp)
            full_match, title, sbjct_start, sbjct_end, identities = BlastUtils.parse_blast(blast_tmp, seq_to_blast)
            inter_df.loc[index, 'biomart_full_match'] = full_match
            inter_df.loc[index, 'biomart_title'] = title
            inter_df.loc[index, 'biomart_sbjct_start'] = sbjct_start
            inter_df.loc[index, 'biomart_sbjct_end'] = sbjct_end
            inter_df.loc[index, 'biomart_identities'] = identities
        inter_df.to_csv(output_file)

    def add_full_mRNA (self):
        blast_result = self.blast_result
        blast_db_fasta = self.blast_db_fasta_biomart
        output_file = self.blast_with_mRNA

        record_dict = SeqIO.index(blast_db_fasta, "fasta")
        inter_df = pd.read_csv(blast_result)
        for index, row in inter_df.iterrows():
            try:
                key = row['biomart_title'].split()[0]
                full_mrna = str(record_dict[key].seq)
                inter_df.loc[index, 'full_mrna'] = full_mrna
                inter_df.loc[index, 'blast_status'] = "OK" if full_mrna.find(row['target sequence']) != -1 else "Not OK"
            except AttributeError:
                inter_df.loc[index, 'blast_status'] = "Not OK"
        inter_df.to_csv(output_file)



    def file_formatting (self):
        in_df = pd.read_csv(self.blast_with_mRNA)
        in_df['Source'] = "Unambiguous_Identification"
        in_df['Organism'] = self.organism
        in_df['Creation_time'] = JsonLog.get_creation_time()

        # Take only the valid rows: Status=OK
        valid_rows = in_df['blast_status']=="OK"
        JsonLog.add_to_json("valid blast results", sum(valid_rows))
        JsonLog.add_to_json("invalid blast results",(in_df.shape[0] - sum(valid_rows)))
        in_df = in_df[valid_rows]

        #Remove miRNA with XXX
        rows_without_XXX = in_df['miRNA sequence'].apply(lambda x: x.find('X')==-1)
        JsonLog.add_to_json("valid miRNA" ,sum(rows_without_XXX))
        JsonLog.add_to_json("invalid miRNA (xxx)",(in_df.shape[0] - sum(rows_without_XXX)))
        in_df = in_df[rows_without_XXX]

        # Choose the necessary columns
        in_df_filter = in_df.filter(
            ['Source', 'Organism', 'miRNA ID', 'miRNA sequence', 'target sequence', 'number of reads',
             'biomart_title', 'biomart_sbjct_start', 'biomart_sbjct_end', 'full_mrna'], axis=1)

        in_df_filter.rename(columns={'miRNA ID': 'microRNA_name', 'biomart_title': 'mRNA_name', 'biomart_sbjct_start': 'mRNA_start',
                              'biomart_sbjct_end': 'mRNA_end'}, inplace=True)


        # reset the index
        in_df_filter.reset_index(drop=True, inplace=True)
        # self.log.append("miRNA statistics")
        # self.log.append(dict(in_df_filter['microRNA_name'].value_counts()))

        # save to file
        if not self.run_as_service:
            in_df_filter.to_csv(self.final_output)
        else:
            return in_df_filter

    def run (self):
        #####################################################
        # build the blast DB
        #####################################################
        self.create_blast_db()

        #####################################################
        # Run the blast against biomart
        #####################################################
        self.blast_mRNA()

        #####################################################
        # Add the full mRNA seq and report the status
        #####################################################
        self.add_full_mRNA()



def main ():
    interaction_file = str(Path("Raw/1-s2.0-S1097276514003566-mmc3.xls"))
    log_dir = "Datafiles_Prepare/Logs/"

    #####################################################
    # Celegans
    #####################################################
    p_dir = Path ("Celegans/Raw")
    JsonLog.set_filename(Path(log_dir)/"Unambiguous_Identification_Celegans.json")
    JsonLog.add_to_json('file name', interaction_file)
    ce = Unambiguous_Identification ("elegans", interaction_file, p_dir)
    ce.run()
    ce.file_formatting()


    #####################################################
    # Mouse
    #####################################################
    p_dir = Path ("Mouse/Raw")
    JsonLog.set_filename(Path(log_dir) /"Unambiguous_Identification_mouse.json")
    JsonLog.add_to_json('file name', interaction_file)
    mo = Unambiguous_Identification ("mouse", interaction_file, p_dir)
    mo.run()
    mo.file_formatting()

    #####################################################
    # Human
    #####################################################
    p_dir = Path ("Human/Raw")
    JsonLog.set_filename(Path(log_dir) / "Unambiguous_Identification_human.json")
    JsonLog.add_to_json('file name', interaction_file)
    ho = Unambiguous_Identification ("human", interaction_file, p_dir)
    ho.run()
    ho.file_formatting()




if __name__ == "__main__":
    main()


