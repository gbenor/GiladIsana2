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
from pathlib import Path


CONFIG = {
    'minimum_pairs_for_interaction': 11
}

files =[{"duplex_method": "vienna",
         "file_name" :"Data/Human/Parsed/human_clash_data_utr3_with_biomart_seq_vienna.csv"},
        {"duplex_method": "miranda",
                         "file_name":"Data/Human/Parsed/human_clash_data_utr3_with_biomart_seq_miranda.csv"}]

res =[]
for f in files:
    file_path = Path(f["file_name"])
    human_clash_data = pd.read_csv(file_path)
    human_clash_data ['GI_seed_type'] = np.nan

    for index, row in human_clash_data.iterrows():

        if f["duplex_method"] == "vienna" :
            dp = ViennaRNADuplex (row.miRNA_seq, row.mRNA_seq_extended)
        if f["duplex_method"] == "miranda":
            dp = MirandaDuplex(row.miRNA_seq, row.mRNA_seq_extended, "Data/Human/Parsed")


        c_seed = dp.IRP.extract_seed()

        seed_feature = SeedFeatures (c_seed)
        seed_feature.extract_seed_features()
        human_clash_data.loc[index, 'GI_seed_type'] = seed_feature.canonic

        #print seed_feature
        if seed_feature.canonic == "None" and len (seed_feature.seed_type) > 0 :
            raise Exception ("seed type error")


    # Save only the valid seeds
    #--------------------------
    seed_filter_df = human_clash_data[human_clash_data.GI_seed_type!="None"]
    seed_filter_df.to_csv(file_path.parent / Path(file_path.stem + "_valid_seeds" + file_path.suffix))
    res.append("Total rows with valid seeds: {}".format(seed_filter_df.shape[0]))

print res
