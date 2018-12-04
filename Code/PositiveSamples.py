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
from MrnaFeatures import *
from EnergyAccess import *
from DatafilesPreparation import *
from pathlib import Path
from Utils import extract_site
from NegativeSamples import *

CONFIG = {
    'minimum_pairs_for_interaction': 11
}

files =[{"duplex_method": "vienna",
        "file_name" :"Data/Human/Parsed/human_clash_data_utr3_vienna_valid_seeds_biomart_no_CDS_biomart_only"}
        # {"duplex_method": "miranda",
        #  "file_name":"Data/Human/Parsed/human_clash_data_utr3_with_biomart_seq_miranda_valid_seeds.csv"}
]







class TooManyMatches(Exception):
    pass



def extract_features (filename, duplex_method):
    df = pd.DataFrame()
    i=0
    human_clash_data = pd.read_csv(filename)

    for index, row in human_clash_data.iterrows():
        # print row.mRNA_name
        # i+=1
        # if i > 5:
        #     break
        # if row.ensg != "ENSG00000165572":
        #     continue

        ############################################################################
        # Generate the duplex
        ############################################################################
        if duplex_method == "vienna" :
            dp = ViennaRNADuplex (row.miRNA_seq, row.mRNA_seq_extended)
        if duplex_method == "miranda":
            dp = MirandaDuplex(row.miRNA_seq, row.mRNA_seq_extended, "Data/Human/Parsed")

        ############################################################################
        # Calc the site of interaction
        ############################################################################
        hint = None
        if duplex_method == "vienna" :
            hint = dp.mrna_coor
        site, site_start_loc, site_end_loc =  extract_site (dp.IRP.site, row.mRNA_seq_extended, row.full_mrna_seq_match_start, row.full_mrna_seq, row.mRNA_name, hint)
        print "site start: " + str(site_start_loc)
        dp.IRP.set_site(site[::-1])

        print (dp)


        #####################
        # Seed XX Features
        #####################
        c_seed = dp.IRP.extract_seed()
        seed_feature = SeedFeatures (c_seed)
        seed_feature.extract_seed_features()


        #####################
        # Matching & Pairing Features
        #####################
        matching_features = MatchingFeatures(dp.IRP)
        matching_features.extract_matching_features()

        #####################
        # mRNA Features
        #####################
        mrna_features = MrnaFeatures(row.miRNA_seq, site, site_start_loc, site_end_loc, row.full_mrna_seq, flank_number=70)
        mrna_features.extract_mrna_features()

        #####################
        # Free energy & Accessibility
        #####################
        enrgy_access = EnergyAccess (dp, site, site_start_loc, site_end_loc, row.full_mrna_seq)

        #####################
        # Data frame construction
        #####################
        features = [seed_feature.get_features(),
                    matching_features.get_features(),
                    mrna_features.get_features(),
                    enrgy_access.get_features()
               ]
        features_df = reduce (lambda x,y: pd.concat([x, y], axis=1, sort=False), features)

        information = pd.DataFrame(row).transpose()
        information['site_start'] = site_start_loc
        information_df = information [['seq_ID', 'microRNA_name','miRNA_seq', 'mRNA_name', 'site_start', 'mRNA_seq_extended', 'full_mrna_seq']]

        new_line_df = reduce (lambda x,y: pd.concat([x.reset_index(drop=True), y.reset_index(drop=True)], axis=1, sort=False), [information_df, features_df])

        # cons = conservation(mrna, mr_site_loc)

        # mmp = pd.DataFrame([self.mmp_dic])
        # mpc = pd.DataFrame([self.mpc_dic])

        df = pd.concat([df, new_line_df], sort=False)

    df.reset_index(drop=True, inplace=True)
    return df




#     # Save only the valid seeds
#     #--------------------------
#     seed_filter_df = human_clash_data[human_clash_data.GI_seed_type!="None"]
#     seed_filter_df.to_csv(file_path.parent / Path(file_path.stem + "_valid_seeds" + file_path.suffix))
#     res.append("Total rows with valid seeds: {}".format(seed_filter_df.shape[0]))
#
# print res



def main():
    data_dir = Path("Data/Human/Parsed/")

    #################################
    # Vienna
    #################################
    vienna_pos_input = data_dir / "human_clash_data_utr3_vienna_valid_seeds_biomart_no_CDS_biomart_only.csv"
    vienna_neg_input = data_dir / "human_vienna_neg_tmp.csv"
    vienna_pos_file = data_dir / "Samples/vienna_pos.csv"
    vienna_neg_file = data_dir / "Samples/vienna_neg.csv"
    vienna_positive_df = extract_features (vienna_pos_input, "vienna")
    vienna_positive_df.to_csv(vienna_pos_file)
    ns = NegativeSamples("Data/Human/Raw/mature.fa")
    ns.generate_negative_samples (vienna_pos_file, vienna_neg_input, duplex_method="vienna", num_of_tries=10000)
    vienna_negative_df = extract_features (vienna_neg_input, "vienna")
    vienna_negative_df.to_csv(vienna_neg_file)





if __name__ == "__main__":
    main()

