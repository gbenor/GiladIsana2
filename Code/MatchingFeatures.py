from collections import Counter
from InteractionRichPresentation import *
import pandas as pd

class MatchingFeatures(object):

    def __init__(self, irp):
        self.irp = irp
        self.irp.replace_T_U()

    def extract_matching_features(self):
        self.miRNA_match_position()
        self.miRNA_pairing_count()


    def get_features(self):
        f = [self.mmp_dic, self.mpc_dic]
        df = map(lambda x: pd.DataFrame([x]), f)
        r =  reduce (lambda x,y: pd.concat([x, y], axis=1, sort=False), df)
        return r



    def miRNA_match_position(self):  # 20
        AU = ['AU', 'UA']
        GC = ['GC', 'CG']
        GU = ['GU', 'UG']

        mmp_dic = {}
        pair_list = list(self.irp.mir_pairing_iterator())
        for i in range (1,21):
            key = 'miRNA_match_position' + str(i)
            try:
                pair = pair_list[i-1]
            except IndexError:
                pair = '  '

            if pair in AU:
                mmp_dic[key] = 2
            elif pair in GC:
                mmp_dic[key] = 1
            elif pair in GU:
                mmp_dic[key] = 3
            elif ' ' in pair:
                mmp_dic[key] = 5
            else:
                mmp_dic[key] = 4

        self.mmp_dic = mmp_dic


    def miRNA_pairing_count(self):  # 6*3=18
        AU = ['AU', 'UA']
        GC = ['GC', 'CG']
        GU = ['GU', 'UG']
        MM = ['AA', 'AG', 'AC', 'UU', 'UC', 'GA', 'GG', 'CA', 'CU', 'CC']

        mpc_dic = {'Seed_GC': 0,
                   'Seed_AU': 0,
                   'Seed_GU': 0,
                   'Seed_mismatch': 0,
                   'Seed_bulge': 0,
                   'Seed_bulge_nt': 0,
                   'Total_GC': 0,
                   'Total_AU': 0,
                   'Total_GU': 0,
                   'Total_mismatch': 0,
                   'Total_bulge': 0,
                   'Total_bulge_nt': 0,
                   'X3p_GC': 0,
                   'X3p_AU': 0,
                   'X3p_GU': 0,
                   'X3p_mismatch': 0,
                   'X3p_bulge': 0,
                   'X3p_bulge_nt': 0}

        i = 0
        for pair in self.irp.mir_pairing_iterator():
            i += 1
            if pair in AU:
                mpc_dic['Total_AU'] += 1
                if 0 < i < 9:
                    mpc_dic['Seed_AU'] += 1
                if i >=9 :
                    mpc_dic['X3p_AU'] += 1
            elif pair in GC:
                mpc_dic['Total_GC'] += 1
                if 0 < i < 9:
                    mpc_dic['Seed_GC'] += 1
                if i >= 9:
                    mpc_dic['X3p_GC'] += 1
            elif pair in GU:
                mpc_dic['Total_GU'] += 1
                if 0 < i < 9:
                    mpc_dic['Seed_GU'] += 1
                if i >= 9:
                    mpc_dic['X3p_GU'] += 1
            elif pair in MM:
                mpc_dic['Total_mismatch'] += 1
                if 0 < i < 9:
                    mpc_dic['Seed_mismatch'] += 1
                if i >= 9:
                    mpc_dic['X3p_mismatch'] += 1
            elif ' ' in pair:
                mpc_dic['Total_bulge_nt'] += 1
                if 0 < i < 9:
                    mpc_dic['Seed_bulge_nt'] += 1
                if i >= 9:
                    mpc_dic['X3p_bulge_nt'] += 1


        mir_total_bulges = self.irp.count_bulges()[0]
        mir_seed_bulges = self.irp.extract_seed().count_bulges()[0]
        mir_X3p_bulges = mir_total_bulges - mir_seed_bulges
        mpc_dic['Total_bulge'] = mir_total_bulges
        mpc_dic['Seed_bulge'] = mir_seed_bulges
        mpc_dic['X3p_bulge'] = mir_X3p_bulges

        #original code:
        # mirna = 'A' + mirna
        # for i in range(len(mirna) + 1)[1:]:
        #     if mirna[-i] == '-' and mirna[-i - 1] != '-':
        #         mpc_dic['Total_bulge'] += 1
        #         if -9 < i < -1:
        #             mpc_dic['Seed_bulge'] += 1
        #         if i <= -9:
        #             mpc_dic['X3p_bulge'] += 1
        self.mpc_dic = mpc_dic


    def tostring(self):
        mmp = pd.DataFrame([self.mmp_dic])
        mpc = pd.DataFrame([self.mpc_dic])
        pd.set_option('display.max_columns', None)

        classstr = ""

        classstr = classstr + str(mmp) + "\n"
        classstr = classstr + str(mpc) + "\n"

        return classstr


    def __str__(self):
        return self.tostring()
