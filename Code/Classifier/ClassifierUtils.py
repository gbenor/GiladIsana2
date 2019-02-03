from pathlib import Path

import JsonLog
import itertools
from multiprocessing import Pool
import multiprocessing

import pandas as pd
import numpy as np

from sklearn.model_selection import train_test_split
from sklearn.metrics import f1_score, accuracy_score
from sklearn.model_selection import GridSearchCV
from scipy.stats import norm
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.naive_bayes import MultinomialNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import f1_score, accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.ensemble import GradientBoostingClassifier, RandomForestClassifier, AdaBoostClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import confusion_matrix






def corespondeing_neg_file (pos_file, l):
    term_to_find = str(pos_file.stem).split("_2019")[0]
    term_to_find = term_to_find.replace("pos", "neg")
    term_to_find = term_to_find.replace("_Data","")
    term_to_find = "*{}*".format(term_to_find)
    for x in l:
        if x.match(term_to_find):
            return x
    raise Exception ("Don't find the corespondeing_neg_file {}".format(pos_file))





def train_test_prepare(pos_files, list_of_neg_files, test_size=0.2, r_state=42):
    f0_pos = pos_files[0]
    f0_neg = corespondeing_neg_file(f0_pos, list_of_neg_files)
    f1_pos = pos_files[1]
    f1_neg = corespondeing_neg_file(f1_pos, list_of_neg_files)

    col_to_drop = ['Unnamed: 0', 'Source', 'Organism', 'microRNA_name', 'miRNA sequence',
                   'target sequence', 'number of reads', 'mRNA_name', 'mRNA_start',
                   'mRNA_end', 'full_mrna','site_start']

    if f0_pos == f1_pos:
        pos = pd.read_csv(f0_pos)
        neg = pd.read_csv(f0_neg)
        X = pd.concat([pos, neg])
        X.reset_index(drop=True, inplace=True)
        X.drop(col_to_drop, axis=1, inplace=True)
        for c in X.columns:
            if c.find ("Unnamed")!=-1:
                print (" ***************** delete")

                X.drop([c], axis=1, inplace=True)

        y_pos = pd.DataFrame(np.ones((pos.shape[0], 1)))
        y_neg = pd.DataFrame(np.zeros((neg.shape[0], 1)))
        Y = pd.concat([y_pos, y_neg])
        Y.reset_index(drop=True, inplace=True)

        X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=test_size ,random_state=r_state)
        return  X_train, X_test, y_train.values.ravel(), y_test.values.ravel()
    else:
        pos_train = pd.read_csv(f0_pos)
        neg_train = pd.read_csv(f0_neg)
        pos_test = pd.read_csv(f1_pos)
        neg_test = pd.read_csv(f1_neg)

        X_train = pd.concat([pos_train, neg_train])
        X_train.reset_index(drop=True, inplace=True)
        X_train.drop(col_to_drop, axis=1, inplace=True)

        y_pos = pd.DataFrame(np.ones((pos_train.shape[0], 1)))
        y_neg = pd.DataFrame(np.zeros((neg_train.shape[0], 1)))
        y_train = pd.concat([y_pos, y_neg])
        y_train.reset_index(drop=True, inplace=True)

        X_test = pd.concat([pos_test, neg_test])
        X_test.reset_index(drop=True, inplace=True)
        X_test.drop(col_to_drop, axis=1, inplace=True)

        y_pos = pd.DataFrame(np.ones((pos_test.shape[0], 1)))
        y_neg = pd.DataFrame(np.zeros((neg_test.shape[0], 1)))
        y_test = pd.concat([y_pos, y_neg])
        y_test.reset_index(drop=True, inplace=True)

        for c in X_test.columns:
            if c.find ("Unnamed")!=-1:
                print (" ***************** delete")
                X_test.drop([c], axis=1, inplace=True)

        for c in X_train.columns:
            if c.find ("Unnamed")!=-1:
                print (" ***************** delete")

                X_train.drop([c], axis=1, inplace=True)

        return  X_train, X_test, y_train.values.ravel(), y_test.values.ravel()

def model_run(training_config, train_X, test_x, train_y, test_y, scoring = "accuracy"):
    best_est = []
    exp_to_run = training_config.keys()
    results = {'name': [], 'f1': [], 'accuracy': []}
    for name in exp_to_run:
        conf = training_config[name]
        clf = conf['clf']
        parameters = conf['parameters']
        n_jobs = conf.get('n_jobs', 1)

        print('=' * 20)
        print('Starting training:', name)
        grid_obj = GridSearchCV(clf, parameters, scoring=scoring, cv=4, n_jobs=n_jobs, verbose=3)

        print('Number of Features:', train_X.columns.shape[0])
        grid_obj = grid_obj.fit(train_X, train_y)
        best_clf = grid_obj.best_estimator_
        best_est.append(best_clf)

        print('Best classifier:', repr(best_clf))
        model = best_clf.fit(train_X, train_y)
        pred_y = model.predict(test_x)

        f1 = f1_score(test_y, pred_y)
        acc = accuracy_score(test_y, pred_y)
        results['name'].append(name)
        results['f1'].append(f1)
        results['accuracy'].append(acc)
        JsonLog.add_to_json("{}_model_params".format(name), model.get_params())

        try:
            feature_importances = pd.DataFrame(best_clf.feature_importances_, index=train_X.columns,
                                               columns=['importance'])
            a = feature_importances.sort_values('importance', ascending=False)
            print (a.head(10))
            JsonLog.add_to_json("{}_feature_importances".format(name), a.to_dict())


        except Exception:
            print ("This model has no attribute 'feature_importances_")

        print (results)
    JsonLog.add_to_json("results", results)

    return best_est


def main():
    #  sudo /sbin/service sshd start
    log_dir = Path("Data/Results")

    input_pos = Path("Data/Features/CSV")
    input_neg = Path("Data/Features/CSV/Neg")
    all_neg = list(input_neg.iterdir())

    all_pos = list(input_pos.iterdir())
    all_pos_valid = [x for x in all_pos if x.match("*_pos_valid_seeds_*")]
    all_pos_vienna = [x for x in all_pos_valid if x.match("*vienna*")]
    all_pos_miranda = [x for x in all_pos_valid if x.match("*miranda*")]

    vienna_all_options = list(itertools.product(all_pos_vienna, all_pos_vienna))

    training_config = {
        'rf': {
            'clf': RandomForestClassifier(),
            'parameters': {
                'n_estimators': [10, 50, 200, 500],
                'criterion': ['gini', 'entropy'],
                'max_depth': [2, 4, 10, 20],
                'min_samples_leaf': [2, 3],
            },
            'n_jobs': 4,
            'one_hot': False
        }
        # 'rf': {
        #     'clf': RandomForestClassifier(),
        #     'parameters': {
        #         'n_estimators': [ 500],
        #         'criterion': [ 'entropy'],
        #         'max_depth': [20],
        #         'min_samples_leaf': [ 3],
        #     },
        #     'n_jobs': 4,
        #     'one_hot': False
        # }

    }
    dm="vienna"
    for files_pair in vienna_all_options:
        source =[]
        for f in files_pair:
            tmp = str(f.stem).split("_2019")[0]
            source.append(tmp.split("_Data")[0])

        json_name ="{}____{}.json".format(source[0], source[1])
        print (json_name)

        JsonLog.set_filename(Path(log_dir) / json_name)
        JsonLog.add_to_json("source",source[0])
        JsonLog.add_to_json("test",source[1])
        JsonLog.add_to_json("duplex_method", dm)


        X_train, X_test, y_train, y_test = train_test_prepare(files_pair, all_neg, test_size=0.2, r_state=42)
        JsonLog.add_to_json("train size", X_train.shape[0])
        JsonLog.add_to_json("test size", X_test.shape[0])
        best_est = model_run(training_config, X_train, X_test, y_train, y_test, scoring="accuracy")

    # valid_seeds_files = [x for x in all_files if x.match("*_pos_valid_seeds_*")]
    # files_with_params = [(x, x.name.split("_")[0], x.name.split("_Data_")[1].split("_")[0]) for x in valid_seeds_files]
    #
    # # output_dir = Path("Data/Features/CSV")
    # # log_dir = Path("Data/Features/Logs")
    #
    #
    # duplex_method = CONFIG['duplex_method']
    # f = list(input_dir.iterdir())
    # all_work = list(itertools.product(f, duplex_method))
    #
    # p = Pool(CONFIG['max_process'])
    # p.map(worker, all_work)
    # p.close()
    # p.join()
    #


if __name__ == "__main__":
    main()

