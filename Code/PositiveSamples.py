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
from pathlib import Path
from Utils import *
from NegativeSamples import *
import JsonLog



CONFIG = {
    'minimum_pairs_for_interaction': 11,
    'duplex_method' :  ["vienna", "miranda"]
 #   'duplex_method': ["vienna", "miranda", "vienna_with_constraint"]

    #'duplex_method' :  ["miranda"]
}



def extract_features (row, duplex_method, minimum_pairs_for_interaction):
    print row.mRNA_name
    # if not ("ENSMUSG00000028655|ENSMUST00000030408" in row.mRNA_name):
    #     return (False, np.nan, np.nan)




    ############################################################################
    # Generate the duplex
    ############################################################################
    if duplex_method == "vienna" :
        dp = ViennaRNADuplex (row['miRNA sequence'], row['target sequence'])
    if duplex_method == "vienna_with_constraint" :
        dp = ViennaRNADuplex (row['miRNA sequence'], row['target sequence'], constraint=".||||") #nt2-5 must be paired
    if duplex_method == "miranda":
        try:
            dp = MirandaDuplex(row['miRNA sequence'], row['target sequence'], "Data/Human/Parsed")
        except NoMirandaHits:
            return (False, np.nan, np.nan)

    if dp.num_of_pairs < minimum_pairs_for_interaction:
        return (False, np.nan, np.nan)

    ############################################################################
    # Calc the site of interaction
    ############################################################################
    site, site_start_loc, site_end_loc = \
        extract_site (dp.IRP.site, row['target sequence'], row['mRNA_start'],
                      row['full_mrna'], row['mRNA_name'], dp.mrna_coor)
    print "site start: " + str(site_start_loc)
    dp.IRP.set_site(site[::-1])

    print (dp)


    #####################
    # Seed XX Features
    #####################
    c_seed = dp.IRP.extract_seed()
    seed_feature = SeedFeatures (c_seed)
    seed_feature.extract_seed_features()
    valid_seed = (len(seed_feature.seed_type) > 0)

    #####################
    # Matching & Pairing Features
    #####################
    matching_features = MatchingFeatures(dp.IRP)
    matching_features.extract_matching_features()

    #####################
    # mRNA Features
    #####################
    mrna_features = MrnaFeatures(row['miRNA sequence'], site, site_start_loc, site_end_loc,
                                 row['full_mrna'], flank_number=70)
    mrna_features.extract_mrna_features()

    #####################
    # Free energy & Accessibility
    #####################
    enrgy_access = EnergyAccess (dp, site, site_start_loc, site_end_loc, row['full_mrna'])

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
    try:
        information = reduce (lambda x,y: pd.concat([x.reset_index(drop=True), y.reset_index(drop=True)], axis=1, sort=False), [information, dp.get_features()])
    except AttributeError:
        # 'MirandaDuplex' object has no attribute 'get_features'
        pass

    new_line_df = reduce (lambda x,y: pd.concat([x.reset_index(drop=True), y.reset_index(drop=True)], axis=1, sort=False), [information, features_df])
    # cons = conservation(mrna, mr_site_loc)

    # mmp = pd.DataFrame([self.mmp_dic])
    # mpc = pd.DataFrame([self.mpc_dic])

    return (True, valid_seed, new_line_df)

#     df = pd.concat([df, new_line_df], sort=False)
# #
# # df.reset_index(drop=True, inplace=True)
# return df





def main():
    #  sudo /sbin/service sshd start

    input_dir = Path("Data/Datafiles_Prepare/CSV")
    output_dir = Path("Data/Features/CSV")
    log_dir = Path("Data/Features/Logs")


    duplex_method = CONFIG['duplex_method']
    for f in input_dir.iterdir():
        print ("*"*80)
        print (f)
        in_df = pd.read_csv(f)
        pos_valid_seed_df = pd.DataFrame()
        pos_all_seed_df = pd.DataFrame()

        for dm in duplex_method:
            invalid_duplex = 0
            i=0
            for index, row in in_df.iterrows():
                i+=1

                if i>3:
                    continue

                valid_duplex, valid_seed, feature_row = extract_features(row, dm, CONFIG['minimum_pairs_for_interaction'])
                if not valid_duplex:
                    invalid_duplex+=1
                    continue
                pos_all_seed_df = pd.concat([pos_all_seed_df, feature_row], sort=False)
                if valid_seed:
                    pos_valid_seed_df = pd.concat([pos_valid_seed_df, feature_row], sort=False)

            pos_all_seed_df.reset_index(drop=True, inplace=True)
            pos_valid_seed_df.reset_index(drop=True, inplace=True)
            pos_all_seed_file = filename_date_append(output_dir / "{}_{}_pos_all_seeds.csv".format(f.stem, dm))
            pos_valid_seed_file = filename_date_append(output_dir / "{}_{}_pos_valid_seeds.csv".format(f.stem, dm))
            pos_all_seed_df.to_csv(pos_all_seed_file)
            pos_valid_seed_df.to_csv(pos_valid_seed_file)

            JsonLog.set_filename(Path(log_dir) / "{}.json".format(pos_all_seed_file.stem))
            JsonLog.add_to_json("source", f.stem)
            JsonLog.add_to_json("duplex_method", dm)
            JsonLog.add_to_json("label", "Positive")
            JsonLog.add_to_json("input size", in_df.shape[0])
            JsonLog.add_to_json("invalid duplex", invalid_duplex)
            JsonLog.add_to_json("output size", pos_all_seed_df.shape[0])

            JsonLog.set_filename(Path(log_dir) / "{}.json".format(pos_valid_seed_file.stem))
            JsonLog.add_to_json("input size", in_df.shape[0])
            JsonLog.add_to_json("invalid duplex", invalid_duplex)
            JsonLog.add_to_json("output size", pos_valid_seed_df.shape[0])

            #tbd : json update invalid_duplex
            print ("invalid duplex {}".format(invalid_duplex))
    exit (7)

    negative_seq = True
    vienna = False
    vienna_with_constraint = False
    miranda = False
    vienna_all = True
    vienna_with_constraint_all = False
    miranda_all = True

    #################################
    # Negative sequences
    #################################
    # Generate miRNA and mRNA to be the base for negative samples

    # 1. with no constraints: seed, num of pairs
    # 2. constraints: num of pairs
    # 3. constraints: num of pairs, seed (i think that support vienna is suffient)

    neg_input1 = input_dir / "human_neg_tmp1.csv"
    neg_input2 = input_dir / "human_neg_tmp2.csv"
    neg_input3 = input_dir / "human_neg_tmp3.csv"
    if negative_seq:
        pos_input = input_dir / "human_clash_data_utr3_vienna_all_seeds_biomart_no_CDS_biomart_only.csv"
        ns = NegativeSamples("Data/Human/Raw/mature.fa")
        ns.generate_negative_seq (pos_input, neg_input1,
                                  num_of_pairs=0,
                                  check_seed=False,
                                  seq_per_row=3)

        ns.generate_negative_seq(pos_input, neg_input2,
                                 num_of_pairs=11,
                                 check_seed=False,
                                 seq_per_row=3)

        ns.generate_negative_seq(pos_input, neg_input3,
                                 num_of_pairs=11,
                                 check_seed=True,
                                 seq_per_row=3)

    exit (7)
    #################################
    # Vienna
    #################################
    if vienna:
        vienna_pos_input = input_dir / "human_clash_data_utr3_vienna_valid_seeds_biomart_no_CDS_biomart_only.csv"
        vienna_pos_file = filename_date_append(input_dir / "Samples/vienna_pos.csv")
        vienna_neg_file = filename_date_append(input_dir / "Samples/vienna_neg.csv")
        vienna_positive_df = extract_features (vienna_pos_input, "vienna")
        vienna_positive_df.to_csv(vienna_pos_file)
        for i in range(3):
            neg_seq = filename_suffix_append(input_dir / "human_neg_tmp.csv", str(i))
            neg_sample = filename_suffix_append(vienna_neg_file, str(i))

            vienna_negative_df = extract_features (neg_seq, "vienna")
            vienna_negative_df.to_csv(neg_sample)

    #################################
    # Vienna with constraint
    #################################
    if vienna_with_constraint:
        vienna_pos_input = input_dir / "human_clash_data_utr3_vienna_valid_seeds_biomart_no_CDS_biomart_only.csv"
        vienna_neg_input = input_dir / "human_vienna_neg_tmp.csv"
        vienna_pos_file = filename_date_append(input_dir / "Samples/vienna_constraint_pos.csv")
        vienna_neg_file = filename_date_append(input_dir / "Samples/vienna_constraint_neg.csv")
        vienna_positive_df = extract_features(vienna_pos_input, "vienna_with_constraint")
        vienna_positive_df.to_csv(vienna_pos_file)
        ns = NegativeSamples("Data/Human/Raw/mature.fa")
        ns.generate_negative_seq(vienna_pos_file, vienna_neg_input, duplex_method="vienna_with_constraint", num_of_tries=10000)
        vienna_negative_df = extract_features(vienna_neg_input, "vienna")
        vienna_negative_df.to_csv(vienna_neg_file)

    #################################
    # miranda
    #################################
    if miranda:
        miranda_pos_input = input_dir / "human_clash_data_utr3_miranda_valid_seeds_biomart_no_CDS_biomart_only.csv"
        miranda_neg_input = input_dir / "human_miranda_neg_tmp.csv"
        miranda_pos_file = filename_date_append(input_dir / "Samples/miranda_pos.csv")
        miranda_neg_file = filename_date_append(input_dir / "Samples/miranda_neg.csv")
        miranda_positive_df = extract_features (miranda_pos_input, "miranda")
        miranda_positive_df.to_csv(miranda_pos_file)
        ns = NegativeSamples("Data/Human/Raw/mature.fa")
        ns.generate_negative_seq (miranda_pos_file, miranda_neg_input, duplex_method="miranda", num_of_tries=10000)
        miranda_negative_df = extract_features (miranda_neg_input, "miranda")
        miranda_negative_df.to_csv(miranda_neg_file)


    #################################
    # Vienna_all
    #################################
    if vienna_all:
        vienna_pos_input = input_dir / "human_clash_data_utr3_vienna_all_seeds_biomart_no_CDS_biomart_only.csv"
        vienna_neg_input = input_dir / "human_vienna_neg_tmp.csv"
        vienna_pos_file = filename_date_append(input_dir / "Samples/vienna_all_pos.csv")
        vienna_neg_file = filename_date_append(input_dir / "Samples/vienna_all_neg.csv")
        vienna_positive_df = extract_features (vienna_pos_input, "vienna")
        vienna_positive_df.to_csv(vienna_pos_file)
        ns = NegativeSamples("Data/Human/Raw/mature.fa")
        ns.generate_negative_seq (vienna_pos_file, vienna_neg_input, duplex_method="vienna", num_of_tries=10000)
        vienna_negative_df = extract_features (vienna_neg_input, "vienna")
        vienna_negative_df.to_csv(vienna_neg_file)

    #################################
    # Vienna with constraint
    #################################
    if vienna_with_constraint_all:
        vienna_pos_input = input_dir / "human_clash_data_utr3_vienna_all_seeds_biomart_no_CDS_biomart_only.csv"
        vienna_neg_input = input_dir / "human_vienna_neg_tmp.csv"
        vienna_pos_file = filename_date_append(input_dir / "Samples/vienna_constraint_all_pos.csv")
        vienna_neg_file = filename_date_append(input_dir / "Samples/vienna_constraint_all_neg.csv")
        vienna_positive_df = extract_features(vienna_pos_input, "vienna_with_constraint")
        vienna_positive_df.to_csv(vienna_pos_file)
        ns = NegativeSamples("Data/Human/Raw/mature.fa")
        ns.generate_negative_seq(vienna_pos_file, vienna_neg_input, duplex_method="vienna_with_constraint", num_of_tries=10000)
        vienna_negative_df = extract_features(vienna_neg_input, "vienna")
        vienna_negative_df.to_csv(vienna_neg_file)

    #################################
    # miranda
    #################################
    if miranda_all:
        miranda_pos_input = input_dir / "human_clash_data_utr3_miranda_all_seeds_biomart_no_CDS_biomart_only.csv"
        miranda_neg_input = input_dir / "human_miranda_neg_tmp.csv"
        miranda_pos_file = filename_date_append(input_dir / "Samples/miranda_all_pos.csv")
        miranda_neg_file = filename_date_append(input_dir / "Samples/miranda_all_neg.csv")
        miranda_positive_df = extract_features (miranda_pos_input, "miranda")
        miranda_positive_df.to_csv(miranda_pos_file)
        ns = NegativeSamples("Data/Human/Raw/mature.fa")
        ns.generate_negative_seq (miranda_pos_file, miranda_neg_input, duplex_method="miranda", num_of_tries=10000)
        miranda_negative_df = extract_features (miranda_neg_input, "miranda")
        miranda_negative_df.to_csv(miranda_neg_file)






if __name__ == "__main__":
    main()

