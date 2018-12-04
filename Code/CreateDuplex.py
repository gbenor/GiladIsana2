import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import pickle
from Bio import SeqIO
from collections import Counter
import RNA
from ViennaRNADuplex import *
from SeedFeatures import *
from MatchingFeatures import *
from MirandaDuplex import *
from DatafilesPreparation import *


CONFIG = {
    'minimum_pairs_for_interaction': 11

}

human_clash_data_utr3 = pd.read_csv("Data/Human/Parsed/human_clash_data_utr3.csv")
miranda_csv_file = "Data/Human/Parsed/human_clash_data_utr3_miranda.csv"
vienna_csv_file = "Data/Human/Parsed/human_clash_data_utr3_vienna.csv"


############################################################################
# Create duplex and filter out duplexes with small number of interactions
############################################################################

human_clash_data_utr3['vienna_num_of_pairs'] = np.nan
human_clash_data_utr3['miranda_num_of_pairs'] = np.nan

for index, row in human_clash_data_utr3.iterrows():
    print row.ensg
    # if row.ensg != "ENSG00000070444":
    #     continue
    try:
        mirnanda_dp = MirandaDuplex(row.miRNA_seq, row.mRNA_seq_extended, "Data/Human/Parsed")
        human_clash_data_utr3.loc[index,'miranda_num_of_pairs'] = mirnanda_dp.num_of_pairs
    except NoMirandaHits:
        human_clash_data_utr3.loc[index,'miranda_num_of_pairs'] = -1

    viennna_dp = ViennaRNADuplex (row.miRNA_seq, row.mRNA_seq_extended)
    human_clash_data_utr3.loc[index, 'vienna_num_of_pairs'] = viennna_dp.num_of_pairs

miranda_df = human_clash_data_utr3[human_clash_data_utr3.miranda_num_of_pairs>=CONFIG['minimum_pairs_for_interaction']]
vienna_df = human_clash_data_utr3[human_clash_data_utr3.vienna_num_of_pairs>=CONFIG['minimum_pairs_for_interaction']]

miranda_df.to_csv(miranda_csv_file)
vienna_df.to_csv(vienna_csv_file)
