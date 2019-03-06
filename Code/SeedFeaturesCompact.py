from collections import Counter
from InteractionRichPresentation import *
import pandas as pd
import pprint


class SeedException (Exception):
    pass



class SeedFeaturesCompact(object):

    def __init__(self, seed):
        self.seed = seed
        self.seed.replace_T_U()


    def extract_seed_features(self):
        self.seed_match_type()


    def get_features(self):
        f = [self.smt_dic]
        df = map(lambda x: pd.DataFrame([x]), f)
        r =  reduce (lambda x,y: pd.concat([x, y], axis=1, sort=False), df)
        return r


    def seed_complementary(self, seq1, seq2):
        count_c = 0
        count_w = 0
        count_mismatch = 0
        c = ['AU', 'UA', 'GC', 'CG']
        w = ['GU', 'UG']

        for i in range(len(seq1)):
            ss = seq1[i] + seq2[i]
            if ss in c:
                count_c += 1
            elif ss in w:
                count_w += 1
            else:
                count_mismatch+=1
        result = {'count_c': count_c, 'count_w': count_w, 'count_mismatch' : count_mismatch}
        return result

    def get_mirna_mrna (self):
        w = ['GU', 'UG']
        SEED_SITE_SIZE = 8
        mr = ''
        mi = ''
        for i in range(SEED_SITE_SIZE):
            if self.seed.mrna_inter[i] != ' ':
                mr += self.seed.mrna_inter[i]
                mi += self.seed.mir_inter[i]
                continue
            if self.seed.mir_bulge[i] != ' ' and self.seed.mrna_bulge[i] != ' ':
                if self.seed.mir_bulge[i] + self.seed.mrna_bulge[i] in w:
                    mr += self.seed.mrna_bulge[i]
                    mi += self.seed.mir_bulge[i]
                else:
                    mi += '-'
                    mr += '-'
                continue
            if self.seed.mir_bulge[i] != ' ':
                mi += self.seed.mir_bulge[i]
                mr += '-'
                continue
            mr += self.seed.mrna_bulge[i]
            mi += '-'
        assert len(mi)==SEED_SITE_SIZE, "mirna size is wrong {}".format(mi)
        assert len(mr)==SEED_SITE_SIZE, "mrna size is wrong {}".format(mr)
        return mi, mr

    def seed_match_type(self):  # 26
        c4 = ['AU', 'UA', 'GC', 'CG']
        smt_dic = {'Seed_match_compact_interactions': 0,
                   'Seed_match_compact_GU': 0,
                   'Seed_match_compact_A': 0,
                   'Seed_match_compact_start': 0,
                   'Seed_match_compact_target_bulge': 0,
                   'Seed_match_compact_mirna_bulge': 0}

        mirna, mrna = self.get_mirna_mrna()
        print ("mrna: " + mrna)
        print ("mirna:" + mirna)

        smt_dic['Seed_match_compact_interactions'] = self.seed_complementary(mirna[0:8], mrna[0:8])['count_c']
        smt_dic['Seed_match_compact_GU'] = self.seed_complementary(mirna[0:8], mrna[0:8])['count_w']

        if smt_dic['Seed_match_compact_interactions'] + smt_dic['Seed_match_compact_GU'] < 6:
            raise SeedException("not valid seed. Seed must have at least 6 combination of interactions and GUs")

        A_mirna = mirna[0] == 'A' and mirna[0] + mrna[0] not in c4
        A_mrna  = mrna[0] == 'A' and mirna[0] + mrna[0] not in c4
        smt_dic['Seed_match_compact_A'] = A_mirna or A_mrna

        bulge = mirna.find("-") or mrna.find("-")
        if not bulge:
            if self.seed_complementary(mirna[0:6], mrna[0:6])['count_mismatch'] == 0:
                smt_dic['Seed_match_compact_start'] = 0
            elif self.seed_complementary(mirna[1:7], mrna[1:7])['count_mismatch'] == 0:
                smt_dic['Seed_match_compact_start'] = 1
            elif self.seed_complementary(mirna[2:8], mrna[2:8])['count_mismatch'] == 0:
                smt_dic['Seed_match_compact_start'] = 2
        else:
            if smt_dic['Seed_match_compact_interactions']<6:
                raise SeedException("not valid seed. Seed with bulge must have at least 6 strong interactions")
            smt_dic['Seed_match_compact_target_bulge'] = min (1, max (0, mirna.find("-")))
            smt_dic['Seed_match_compact_mirna_bulge'] = min (1, max (0, mrna.find("-")))


            smt_dic['Seed_match_compact_start'] = 0 if mirna[0]+mrna[0] in c4 else 1

        #
        # if mrna[0] == 'A':
        #     mir_bulge = mirna[1:8].find("-")
        #     mrna_bulge = mrna[1:8].find("-")
        #     smt_dic['Seed_match_compact_target_bulge'] = min (1, max (0, mrna_bulge))
        #     smt_dic['Seed_match_compact_mirna_bulge'] = min (1, max (0, mir_bulge))

        ###############################################################
        # Update the dict
        ###############################################################
        self.smt_dic = smt_dic



def test_seed (seed, seed_type):
    pp = pprint.PrettyPrinter(indent=4)
    print ("**************************************************")
    print ("Test: " + seed_type)
    print(seed)

    s = SeedFeaturesCompact(seed)
    s.extract_seed_features()


    pp.pprint (s.smt_dic)

    #assert seed_type in s.seed_type, "test error"


def main ():

    s = InteractionRichPresentation("A       ", " GGGGGGG", " CCCCCCC", "C       ")
    test_seed(s, "Seed_match_8merA1")

    s = InteractionRichPresentation("A      C", " GGGGGG ", " CCCCCC ", "C      C")
    test_seed(s, "Seed_match_7merA1")

    s = InteractionRichPresentation("C       ", " GGGGGGG", " CCCCCCC", "C       ")
    test_seed(s, "Seed_match_7mer2")

    s = InteractionRichPresentation("C      C", " GGGGGG ", " CCCCCC ", "C      C")
    test_seed(s, "Seed_match_6mer2")

    s = InteractionRichPresentation("A       ", " GGGGGGG", " CCCCUCC", "C       ")
    test_seed(s, "Seed_match_6mer2GU1")

    s = InteractionRichPresentation("A       ", " GGGGUGG", " CCCCGCC", "C       ")
    test_seed(s, "Seed_match_6mer2GU1")

    s = InteractionRichPresentation("A    C  ", " GGGG GG", " CCCC CC", "C    C  ")
    test_seed(s, "Seed_match_6mer_LP")

    s = InteractionRichPresentation("A       ", " GGGG GG", " CCCC CC", "C    C  ")
    test_seed(s, "Seed_match_6mer_BM")

    s = InteractionRichPresentation("A C       ",
                                " G GGGGGGG",
                                " C CCCCCC",
                                "C       ")
    test_seed(s, "Seed_match_6mer_BT")

    s = InteractionRichPresentation("A C      ",
                                    " G GGGGGG",
                                    " C CCCCCC",
                                    "C        ")
    test_seed(s, "Seed_match_6mer_BT")


if __name__ == "__main__":
    main()

