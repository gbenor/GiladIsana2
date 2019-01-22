from Bio.Blast import NCBIWWW
import JsonLog
import mirBaseUtils as MBU

import Bio
from Bio.Blast import NCBIXML
import pandas as pd
import numpy as np
import datetime
from selenium import webdriver
import time
from Bio import SeqIO
from pathlib import Path



class mirmark(object):


    def __init__(self, input_file, data_dir):
        organism = "human"
        print ("mirMark: a site-level and UTR-level classifier for miRNA target prediction")
        print ("https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0500-5")
        print("#############################################")
        print("Organism: {}".format(organism))
        print("#############################################")

        self.organism = organism
        self.input_file = input_file
        self.data_dir = data_dir

        self.mirmark_Files()
        self.read_paper_data()

    def mirmark_Files (self):
        self.mirmark_with_mrna = "Human/Raw/mirmark_with_mrna_new_method.csv"
        self.mirmark_with_mrna_mirna = "Human/Raw/mirmark_with_mrna_mirna.csv"
        self.final_output = str("Datafiles_Prepare/CSV"/  Path (self.organism + "_Mirmark_Data.csv"))


    def read_paper_data(self):
        # Read the interaction file
        self.in_df = pd.read_csv(self.input_file)
        JsonLog.add_to_json('Num of samples', self.in_df.shape[0])


    def direct_web_query (driver, nm):
        link = "https://www.ncbi.nlm.nih.gov/nuccore/{}?report=fasta&log$=seqview&format=text".format(nm)
        driver.get(link)
        for i in range(5):
            time.sleep(1)
            content = driver.find_element_by_id("maincontent").text
            a = content[content.find("NM_"):]
            desc = a.split('\n')[0]
            seq = ''.join(a.split('\n')[1:])
            if len (seq)> 30:
                break
        return desc, seq


    def insert_full_mrna (self):
        cache={}
        driver = webdriver.Chrome()

        self.in_df['full_mrna_seq'] = np.nan
        self.in_df['mrna_name'] = np.nan
        for index, row in self.in_df.iterrows():
            print ("index: {} \t time {}".format(index, datetime.datetime.now()))
            nm =  row.mRNA_ID
            if nm in cache:
                print ("nm in cache. great success.")
                desc, seq = cache [nm]
            else:
                # web_result = web_query(nm)
                # desc, seq = parse_blast(web_result)
                desc, seq = self.direct_web_query(driver, nm)
                cache[nm] = (desc, seq)
            self.in_df.loc[index, 'full_mrna_seq'] = seq
            self.in_df.loc[index, 'mrna_name'] = desc
            print (desc)
            print (seq[0:30])
            # save each iteration
            self.in_df.to_csv(self.mirmark_with_mrna)
        self.in_df.to_csv(self.mirmark_with_mrna)

    def extract_target (self):
        def seq_slice(row):
            start = row['Start_position']
            end = row['End_position']
            full_seq = row['full_mrna_seq']
            return full_seq[start:end+1]


        in_df = pd.read_csv(self.mirmark_with_mrna)
        in_df['target'] = in_df.apply(seq_slice, axis=1)
        self.inter_df = in_df.copy()


    def file_formatting (self):
        in_df = pd.read_csv(self.mirmark_with_mrna_mirna)
        in_df['Source'] = "Mirmark"
        in_df['Organism'] = self.organism
        in_df['Creation_time'] = JsonLog.get_creation_time()
        in_df['number of reads'] = 1

        JsonLog.add_to_json("valid web-query results", in_df.shape[0])


        #Remove miRNA without match
        rows_with_match = in_df['miRNA_seq'].apply(lambda x: x.find('match')==-1)
        JsonLog.add_to_json("valid miRNA" ,sum(rows_with_match))
        JsonLog.add_to_json("invalid miRNA (no match)",(in_df.shape[0] - sum(rows_with_match)))
        in_df = in_df[rows_with_match]

        #Choose the necessary columns
        in_df_filter = in_df.filter(
            ['Source', 'Organism', 'miR_ID', 'miRNA_seq', 'mRNA_ID', 'target', 'number of reads',
             'Start_position', 'End_position', 'full_mrna_seq'], axis=1)

        in_df_filter.rename(columns={'miR_ID': 'microRNA_name', 'mRNA_ID': 'mRNA_name', 'Start_position': 'mRNA_start',
                              'End_position': 'mRNA_end','full_mrna_seq': 'full_mrna'}, inplace=True)


        # reset the index
        in_df_filter.reset_index(drop=True, inplace=True)

        # save to file
        in_df_filter.to_csv(self.final_output)



    def run (self):
        #####################################################
        # build the blast DB
        #####################################################
        #self.insert_full_mrna()
        self.extract_target()

        #####################################################
        # Add the miRNA
        #####################################################
        print (self.inter_df.columns)
        #self.inter_df['mir_ID'] = self.inter_df['name'].apply(lambda x: x[:x.find("_")])
        self.inter_df = MBU.insert_mirna(self.inter_df, "Human/Raw/mature.fa", "miR_ID")
        self.inter_df.to_csv(self.mirmark_with_mrna_mirna)

        #####################################################
        # File formatting
        #####################################################
        self.file_formatting()


def main ():
    interaction_file = str(Path("Human/Raw/13059_2014_500_MOESM1_ESM.csv"))
    log_dir = "Datafiles_Prepare/Logs/"

    #####################################################
    # Human
    #####################################################
    p_dir = Path ("human/Raw")
    JsonLog.set_filename(Path(log_dir)/"mirmark_human.json")
    JsonLog.add_to_json('file name', interaction_file)
    hu = mirmark (interaction_file, p_dir)
    hu.run()



if __name__ == "__main__":
    main()

