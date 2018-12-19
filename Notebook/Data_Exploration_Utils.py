import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, accuracy_score
from sklearn.model_selection import GridSearchCV

def train_test_prepare(pos, neg, r_state=42):
    X = pd.concat([pos, neg])
    X.reset_index(drop=True, inplace=True)
    X.drop('seq_ID', axis=1, inplace=True)
    X.drop('microRNA_name', axis=1, inplace=True)
    X.drop('miRNA_seq', axis=1, inplace=True)
    X.drop('mRNA_name', axis=1, inplace=True)
    X.drop('site_start', axis=1, inplace=True)
    X.drop('mRNA_seq_extended', axis=1, inplace=True)
    X.drop('full_mrna_seq', axis=1, inplace=True)

    y_pos = pd.DataFrame(np.ones((pos.shape[0], 1)))
    y_neg = pd.DataFrame(np.zeros((pos.shape[0], 1)))
    Y = pd.concat([y_pos, y_neg])
    Y.reset_index(drop=True, inplace=True)

    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.1, random_state=r_state)
    return X_train, X_test, y_train, y_test



from sklearn.cluster import KMeans
def clustering (pos):
    # s = "d=b.loc[:,["
    # for i in range(1, 21):
    #     key = 'miRNA_match_position' + str(i)
    #     s = s + '\'' + key + '\''
    #     if i != 20:
    #         s += ','
    # s += "]]"
    # print (s)

    b = pos.loc[:, ['miRNA_match_position1', 'miRNA_match_position2', 'miRNA_match_position3', 'miRNA_match_position4',
                  'miRNA_match_position5', 'miRNA_match_position6', 'miRNA_match_position7', 'miRNA_match_position8',
                  'miRNA_match_position9', 'miRNA_match_position10', 'miRNA_match_position11', 'miRNA_match_position12',
                  'miRNA_match_position13', 'miRNA_match_position14', 'miRNA_match_position15',
                  'miRNA_match_position16', 'miRNA_match_position17', 'miRNA_match_position18',
                  'miRNA_match_position19', 'miRNA_match_position20']]

    b = b.applymap(lambda x: -10 if x != 5 else 10)

    kmeans = KMeans(n_clusters=5, random_state=0)
    c = kmeans.fit_predict(b)

    b['group'] = c*2
    d = b.sort_values('group')
    return d
    #return d.drop(['group'], axis=1)

import glob
import os
def get_latest_sample_file (method):
    p_pos = "*{}_pos*.csv".format(method)
    p_neg = "*{}_neg*.csv".format(method)
    list_of_files = glob.glob(p_pos)
    latest_pos = max(list_of_files, key=os.path.getctime)
    list_of_files = glob.glob(p_neg)
    latest_neg = max(list_of_files, key=os.path.getctime)

    print (latest_pos)
    print(latest_neg)
    return latest_pos, latest_neg