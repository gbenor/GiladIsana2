import pandas as pd
import numpy as np
from Bio.Blast.Applications import NcbiblastnCommandline
from ViennaRNADuplex import *
from Bio import SeqIO
import random
from SeedFeatures import *

from InteractionRichPresentation import *
from MirandaDuplex import *

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import itertools as it
from itertools import islice
from pathlib import Path
from Utils import *
import PositiveSamples as ps
import JsonLog
import itertools
from multiprocessing import Pool
import multiprocessing




CONFIG = {
    'minimum_pairs_for_interaction': 11,
    'max_process' : 1,
    'max_tries_to_get_valid_seed' : 1000

}


class NegativeSamples(object):

    def __init__(self, mirbase_file, organism):
        self.organism = organism
        self.mirna_prefix_table = {"human" : "hsa",
                                   "mouse" : "mmu",
                                   "elegans" : "cel"}
        self.organism_prefix = self.mirna_prefix_table[organism]
        self.read_mirbase_file(mirbase_file)


    def read_mirbase_file (self, mirbase_file):
        fasta_sequences = SeqIO.parse(open(mirbase_file), 'fasta')
        miRBase_seq_list = []
        miRBase_dic = {}
        for fasta in fasta_sequences:
            mi_name, ma_name, mi_seq = fasta.id, fasta.description.split()[1], str(fasta.seq)
            if mi_name.startswith(self.organism_prefix):
                miRBase_dic[ma_name] = [mi_name, mi_seq]
                miRBase_seq_list.append(mi_seq)
        print ('number of entrences in miRBase {} {}: '.format(self.organism, len(miRBase_dic)))
        self.miRBase_seq_list = miRBase_seq_list
        self.miRBase_dic = miRBase_dic

    #
    # # # 2. the positive pairs in CLASH experiments
    # f = open('data2/S2/S2_final_pairs_dic.pkl')
    # S2_final_pairs_dic = pk.load(f)
    # f.close()

    # # 4. a function used to generate mock mirna which their seed never equal to any mirna in the miRBase.
    def generate_mirna_mock(self, mirna, th=5):
        def seed_equal(a, b):
            return sum(a[i] == b[i] for i in range(len(a)))

        def equal_to_mirbase (mir, th=5):
            for seq in self.miRBase_seq_list:
                seq = list(seq)
                e27 = seed_equal(mir[1:7], seq[1:7])
                if e27 > th:
                    return True
                e38 = seed_equal(mir[2:8], seq[2:8])
                if e38 > th:
                    return True
            return False


        mirna_list_o = list (mirna.replace('T', 'U').upper())
        mirna_list_r = list (mirna.replace('T', 'U').upper())
        equal_to_itself = True

        num_shuffle = 0
        while equal_to_itself or eq_mirbase:
            random.shuffle(mirna_list_r)
            num_shuffle += 1
            if num_shuffle % 10000 == 0:
                print (num_shuffle)
            if num_shuffle > 100000:
                break

            # check if it equals to itself
            e27 = seed_equal(mirna_list_r[1:7], mirna_list_o[1:7])
            e38 = seed_equal(mirna_list_r[2:8], mirna_list_o[2:8])
            equal_to_itself = e27 > th or e38 > th
            # check against mirbase
            eq_mirbase = equal_to_mirbase(mirna_list_r)


        mirna_m = ''.join(mirna_list_r)
        return mirna_m.replace('U','T')


    def valid_negative_seq(self, mir, mrna, num_of_pairs=11, check_seed=True):
        dp = ViennaRNADuplex(mir, mrna)
        if check_seed:
            c_seed = dp.IRP.extract_seed()
            seed_feature = SeedFeatures (c_seed)
            seed_valid = seed_feature.canonic != "None"
        else:
            seed_valid = True

        pairs_valid = dp.num_of_pairs > num_of_pairs

        return seed_valid and pairs_valid

    def generate_negative_seq (self, pos_file, output, num_of_pairs=11, check_seed=True, seq_per_row=3, num_of_tries=10000):
        print ("start negative seq:\nnum_of pairs {}\nseed {}".format(num_of_pairs,check_seed))
        pos = pd.read_csv(pos_file)
        neg = pd.DataFrame()



        row_num = 0

        for index, row in pos.iterrows():

            print (index)
            for k in range(seq_per_row):
                valid = False
                for i in range(num_of_tries):
                    mock_mirna = self.generate_mirna_mock(row.miRNA_seq).replace('U', 'T')
                    full_mrna = row.full_mrna_seq
                    try:
                        valid = self.valid_negative_seq (mock_mirna, full_mrna, num_of_pairs, check_seed)
                        if valid:
                            break
                    except IndexError:
                        valid = False

                if not valid:
                    print ("Couldnt manage to create negative sample to row {}".format(row.seq_ID))
                    continue

                # #add the seq to the Dataframe
                # new_row = pd.DataFrame()
                # new_row.loc[0, 'seq_ID'] = row.seq_ID
                # new_row.loc[0, 'microRNA_name'] = "mock " + row.microRNA_name
                # new_row.loc[0, 'miRNA_seq'] = mock_mirna
                # new_row.loc[0, 'mRNA_name'] = row.mRNA_name
                # new_row.loc[0, 'mRNA_seq_extended'] = row.full_mrna_seq
                # new_row.loc[0, 'full_mrna_seq'] = row.full_mrna_seq
                # new_row.loc[0, 'full_mrna_seq_match_start'] = 0

                # add the seq to the Dataframe
                neg.loc[row_num, 'seq_ID'] = row.seq_ID
                neg.loc[row_num, 'microRNA_name'] = "mock " + row.microRNA_name
                neg.loc[row_num, 'miRNA_seq'] = mock_mirna
                neg.loc[row_num, 'mRNA_name'] = row.mRNA_name
                neg.loc[row_num, 'mRNA_seq_extended'] = row.full_mrna_seq
                neg.loc[row_num, 'full_mrna_seq'] = row.full_mrna_seq
                neg.loc[row_num, 'full_mrna_seq_match_start'] = 0
                row_num+=1

         #   neg = pd.concat([neg, new_row], sort=False)

        neg.reset_index(drop=True, inplace=True)
        neg.to_csv(output)



def worker (w):
    print ("Starting worker #{}".format(multiprocessing.current_process()))
    f, organism, dm = w

    output_dir = Path("Data/Features/CSV/Neg")
    log_dir = Path("Data/Features/Logs/Neg")

    ns = NegativeSamples("Data/all_mature_miRNA.fa", organism)
    in_df = pd.read_csv(f)
    in_df =  in_df.filter(['Source', 'Organism', 'microRNA_name', 'miRNA sequence',
          'target sequence', 'number of reads', 'mRNA_name', 'mRNA_start',
          'mRNA_end', 'full_mrna'], axis=1)

    neg_valid_seed_df = pd.DataFrame()
    neg_invalid_seed_df = pd.DataFrame()

    #
    # i=0
    for index, row in in_df.iterrows():
        # i+=1
        # if i> 3:
        #     continue
        ############################################################
        # Try to get valid duplex and valid/invalid seed
        ############################################################
        got_valid_seed = False
        got_invalid_seed = False
        for j in range (CONFIG['max_tries_to_get_valid_seed']):
            mock_mirna = ns.generate_mirna_mock(row['miRNA sequence'])
            row['miRNA sequence'] = mock_mirna
            valid_duplex, valid_seed, feature_row = ps.extract_features(row, dm, CONFIG['minimum_pairs_for_interaction'])
            if not valid_duplex:
                continue
            if valid_seed and not got_valid_seed:
                neg_valid_seed_df = pd.concat([neg_valid_seed_df, feature_row], sort=False)
                got_valid_seed = True
            if not valid_seed and not got_invalid_seed:
                neg_invalid_seed_df = pd.concat([neg_invalid_seed_df, feature_row], sort=False)
                got_invalid_seed = True
            if got_invalid_seed and got_valid_seed:
                break

    neg_invalid_seed_df.reset_index(drop=True, inplace=True)
    neg_valid_seed_df.reset_index(drop=True, inplace=True)
    neg_invalid_seed_file = filename_date_append(output_dir / "{}_{}_neg_invalid_seeds.csv".format(f.stem.split("_Data_")[0], dm))
    neg_valid_seed_file = filename_date_append(output_dir / "{}_{}_neg_valid_seeds.csv".format(f.stem.split("_Data_")[0], dm))
    neg_invalid_seed_df.to_csv(neg_invalid_seed_file)
    neg_valid_seed_df.to_csv(neg_valid_seed_file)

    JsonLog.set_filename(Path(log_dir) / "{}.json".format(neg_invalid_seed_file.stem))
    JsonLog.add_to_json("source", f.stem)
    JsonLog.add_to_json("duplex_method", dm)
    JsonLog.add_to_json("label", "Negative")
    JsonLog.add_to_json("seed", "invalid")
    JsonLog.add_to_json("input size", in_df.shape[0])
    JsonLog.add_to_json("output size", neg_invalid_seed_df.shape[0])

    JsonLog.set_filename(Path(log_dir) / "{}.json".format(neg_valid_seed_file.stem))
    JsonLog.add_to_json("source", f.stem)
    JsonLog.add_to_json("duplex_method", dm)
    JsonLog.add_to_json("label", "Negative")
    JsonLog.add_to_json("seed", "valid")
    JsonLog.add_to_json("input size", in_df.shape[0])
    JsonLog.add_to_json("output size", neg_valid_seed_df.shape[0])

    print ("Finish worker #{}".format(multiprocessing.current_process()))







def main():
    input_dir = Path("Data/Features/CSV")
    all_files = list(input_dir.iterdir())
    valid_seeds_files = [x for x in all_files if x.match("*_pos_valid_seeds_*")]
    files_with_params = [(x, x.name.split("_")[0], x.name.split("_Data_")[1].split("_")[0]) for x in valid_seeds_files]

    for w in files_with_params:
        worker(w)
    exit(3)
    p = Pool(CONFIG['max_process'])
    p.map(worker, files_with_params)
    p.close()
    p.join()

    exit(7)



if __name__ == "__main__":
    main()
