import pandas as pd
import numpy as np
from SeedFeaturesCompact import *
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
import itertools
from multiprocessing import Pool
import multiprocessing
import pprint





CONFIG = {
    'minimum_pairs_for_interaction': 11,
    'duplex_method' :  ["vienna"],
    'max_process' : 1
 #   'duplex_method': ["vienna", "miranda", "vienna_with_constraint"]

    #'duplex_method' :  ["miranda"]
}



def extract_features (row, duplex_method, minimum_pairs_for_interaction, seed_debug_file=None):

    pp = pprint.PrettyPrinter(indent=4, stream=seed_debug_file)

    print (row.mRNA_name)
    # if not ("WBGene00011573|T07C12.9" in row.mRNA_name):
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

    try:
        seed_feature_compact = SeedFeaturesCompact(c_seed)
        seed_feature_compact.extract_seed_features()
        valid_seed_compact = seed_feature_compact.valid_seed()
    except SeedException:
        valid_seed_compact = False
        #No seed. This is invalid duplex
        return (False, np.nan, np.nan)


    if seed_debug_file is not None:
        if valid_seed != valid_seed_compact:
            print ("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
            print ("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
            seed_debug_file.write("===============================\n")
            seed_debug_file.write("{}\n".format(row.mRNA_name))

            seed_debug_file.write ("valid_seed: {}\n".format(valid_seed))
            seed_debug_file.write ("valid_seed_compact: {}\n".format(valid_seed_compact))
            seed_debug_file.write (str(c_seed))
         #   pp.pprint(seed_feature_compact.smt_dic)

    valid_seed = valid_seed_compact
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
    features = [seed_feature_compact.get_features(),
                seed_feature.get_features(),
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


def worker (w):
    f, dm = w
    output_dir = Path("Data/Features/CSV")
    log_dir = Path("Data/Features/Logs")

    print ("Starting worker #{}".format(multiprocessing.current_process()))
    seed_debug_file_name = output_dir / "{}_seed_debug.txt".format(f.stem)
    seed_debug_file = seed_debug_file_name.open("wb")


    in_df = pd.read_csv(f)
    pos_valid_seed_df = pd.DataFrame()
    pos_all_seed_df = pd.DataFrame()


    invalid_duplex = 0
    # i=0
    for index, row in in_df.iterrows():
        # i+=1
        # if i> 3:
        #     continue
        valid_duplex, valid_seed, feature_row = extract_features(row, dm, CONFIG['minimum_pairs_for_interaction'], seed_debug_file)
        if not valid_duplex:
            invalid_duplex += 1
            continue
        pos_all_seed_df = pd.concat([pos_all_seed_df, feature_row], sort=False)
        if valid_seed:
            pos_valid_seed_df = pd.concat([pos_valid_seed_df, feature_row], sort=False)

    pos_all_seed_df.reset_index(drop=True, inplace=True)
    pos_valid_seed_df.reset_index(drop=True, inplace=True)
    pos_all_seed_file = filename_date_append(output_dir / "{}_{}_all_seeds.xlsx".format(f.stem, dm))
    pos_valid_seed_file = filename_date_append(output_dir / "{}_{}_valid_seeds.xlsx".format(f.stem, dm))

    for df in [pos_valid_seed_df, pos_all_seed_df]:
        for c in df.columns:
            if c.find ("Unnamed")!=-1:
                df.drop([c], axis=1, inplace=True)



    pos_all_seed_df.to_excel(pos_all_seed_file, sheet_name = 'POS')
    pos_valid_seed_df.to_excel(pos_valid_seed_file, sheet_name = 'POS')

    JsonLog.set_filename(Path(log_dir) / "{}.json".format(pos_all_seed_file.stem))
    JsonLog.add_to_json("source", f.stem)
    JsonLog.add_to_json("duplex_method", dm)
    JsonLog.add_to_json("seeds_filter", "all")
    JsonLog.add_to_json("label", "Positive")
    JsonLog.add_to_json("input size", in_df.shape[0])
    JsonLog.add_to_json("invalid duplex", invalid_duplex)
    JsonLog.add_to_json("output size", pos_all_seed_df.shape[0])

    JsonLog.set_filename(Path(log_dir) / "{}.json".format(pos_valid_seed_file.stem))
    JsonLog.add_to_json("source", f.stem)
    JsonLog.add_to_json("duplex_method", dm)
    JsonLog.add_to_json("seeds_filter", "valid")
    JsonLog.add_to_json("label", "Positive")
    JsonLog.add_to_json("input size", in_df.shape[0])
    JsonLog.add_to_json("invalid duplex", invalid_duplex)
    JsonLog.add_to_json("output size", pos_valid_seed_df.shape[0])
    seed_debug_file.close()
    print ("invalid duplex {}".format(invalid_duplex))
    print ("Finish worker #{}".format(multiprocessing.current_process()))




def main():
    #  sudo /sbin/service sshd start

    input_dir = Path("Data/Datafiles_Prepare/CSV")
    # output_dir = Path("Data/Features/CSV")
    # log_dir = Path("Data/Features/Logs")


    duplex_method = CONFIG['duplex_method']
    f = list(input_dir.iterdir())
    all_work = list(itertools.product(f, duplex_method))

    for w in all_work:
        worker(w)
    #
    # p = Pool(CONFIG['max_process'])
    # p.map(worker, all_work)
    # p.close()
    # p.join()
    #




if __name__ == "__main__":
    main()

