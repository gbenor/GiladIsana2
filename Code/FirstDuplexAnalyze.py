import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import pickle
from Bio import SeqIO
from collections import Counter
import RNA
from RichDuplex import *
from SeedFeatures import *
from MatchingFeatures import *
from Miranda import *
from DatafilesPreparation import *


human_clash_data_utr3 = pd.read_csv("Data/Human/Parsed/human_clash_data_utr3_with_biomart_seq_step2.csv")
c=0
for index, row in human_clash_data_utr3.iterrows():
    if index > 50:
        break

    print "ENSG: {} ENST: {} MIR: {}".format(row.ensg, row.enst, row.microRNA_name)
   # mirnanda_dp = Miranda(row.miRNA_seq, row.mRNA_seq_extended, "Data/Human/Parsed")
   # print "MIRANDA DUPLEX \n{}".format (mirnanda_dp)

    if row.ensg=="ENSG00000034063":
        print "debug"

    try:
        dp = RichDuplex (row.miRNA_seq,
                         row.mRNA_seq_extended,
                         row.full_mrna_seq,
                         row.mRNA_start,
                         row.mRNA_end_extended)
    except Exception:
        print "cccccccccc {}".format(c)
        c+=1

    print (dp)

    c_seed = dp.IRP.extract_seed()

    print "Seed:"
    print c_seed
    print "CLASH file seed type: {}".format(row.seed_type)
    seed_feature = SeedFeatures (c_seed)
    seed_feature.extract_seed_features()
    print seed_feature

    t = tuple(dp.IRP.mir_pairing_iterator())
    print "len = {} tuple = {}".format(len(t), t)

    matching_features = MatchingFeatures(dp.IRP)
    matching_features.extract_matching_features()
    print matching_features


    print ("**********************************************************************")
print "total cccccccccc {}".format(c)

print ("gilad")