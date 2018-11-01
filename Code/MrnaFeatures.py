from collections import Counter
from InteractionRichPresentation import *
import pandas as pd

class MrnaFeatures(object):

    def __init__(self, irp):
        self.irp = irp
        self.irp.replace_T_U()

    def extract_mrna_features(self):
        self.miRNA_match_position()
        self.miRNA_pairing_count()

    # # 3. location of target site (1)
    def distance_to_end(mrna, mr_site_loc):  # 1
        dte_dic = {'Dist_to_end': round(float(len(mrna) - int(mr_site_loc[1])) / len(mrna), 4)}
        return dte_dic

    # # 4. target composition (20+20+20)
    def target_composition(mr_site):  # 20
        mrna = mr_site.upper().replace('T', 'U')
        count_A = 0
        count_U = 0
        count_G = 0
        count_C = 0
        count_AA = 0
        count_AU = 0
        count_AG = 0
        count_AC = 0
        count_UA = 0
        count_UU = 0
        count_UG = 0
        count_UC = 0
        count_GA = 0
        count_GU = 0
        count_GG = 0
        count_GC = 0
        count_CA = 0
        count_CU = 0
        count_CG = 0
        count_CC = 0
        for i in range(len(mrna)):
            if mrna[i] == 'A':
                count_A += 1
            elif mrna[i] == 'U':
                count_U += 1
            elif mrna[i] == 'G':
                count_G += 1
            elif mrna[i] == 'C':
                count_C += 1
        for i in range(len(mrna) - 1):
            if mrna[i:i + 2] == 'AA':
                count_AA += 1
            elif mrna[i:i + 2] == 'AU':
                count_AU += 1
            elif mrna[i:i + 2] == 'AG':
                count_AG += 1
            elif mrna[i:i + 2] == 'AC':
                count_AC += 1
            elif mrna[i:i + 2] == 'UA':
                count_UA += 1
            elif mrna[i:i + 2] == 'UU':
                count_UU += 1
            elif mrna[i:i + 2] == 'UG':
                count_UG += 1
            elif mrna[i:i + 2] == 'UC':
                count_UC += 1
            elif mrna[i:i + 2] == 'GA':
                count_GA += 1
            elif mrna[i:i + 2] == 'GU':
                count_GU += 1
            elif mrna[i:i + 2] == 'GG':
                count_GG += 1
            elif mrna[i:i + 2] == 'GC':
                count_GC += 1
            elif mrna[i:i + 2] == 'CA':
                count_CA += 1
            elif mrna[i:i + 2] == 'CU':
                count_CU += 1
            elif mrna[i:i + 2] == 'CG':
                count_CG += 1
            elif mrna[i:i + 2] == 'CC':
                count_CC += 1
        all_monomer_count = count_A + count_U + count_G + count_C
        all_dimer_count = count_AA + count_AU + count_AG + count_AC + \
                          count_UA + count_UU + count_UG + count_UC + \
                          count_GA + count_GU + count_GG + count_GC + \
                          count_CA + count_CU + count_CG + count_CC
        tc_dic = {'Target_A_comp': round(float(count_A) / all_monomer_count, 4),
                  'Target_U_comp': round(float(count_U) / all_monomer_count, 4),
                  'Target_G_comp': round(float(count_G) / all_monomer_count, 4),
                  'Target_C_comp': round(float(count_C) / all_monomer_count, 4),
                  'Target_AA_comp': round(float(count_AA) / all_dimer_count, 4),
                  'Target_AU_comp': round(float(count_AU) / all_dimer_count, 4),
                  'Target_AG_comp': round(float(count_AG) / all_dimer_count, 4),
                  'Target_AC_comp': round(float(count_AC) / all_dimer_count, 4),
                  'Target_UA_comp': round(float(count_UA) / all_dimer_count, 4),
                  'Target_UU_comp': round(float(count_UU) / all_dimer_count, 4),
                  'Target_UG_comp': round(float(count_UG) / all_dimer_count, 4),
                  'Target_UC_comp': round(float(count_UC) / all_dimer_count, 4),
                  'Target_GA_comp': round(float(count_GA) / all_dimer_count, 4),
                  'Target_GU_comp': round(float(count_GU) / all_dimer_count, 4),
                  'Target_GG_comp': round(float(count_GG) / all_dimer_count, 4),
                  'Target_GC_comp': round(float(count_GC) / all_dimer_count, 4),
                  'Target_CA_comp': round(float(count_CA) / all_dimer_count, 4),
                  'Target_CU_comp': round(float(count_CU) / all_dimer_count, 4),
                  'Target_CG_comp': round(float(count_CG) / all_dimer_count, 4),
                  'Target_CC_comp': round(float(count_CC) / all_dimer_count, 4)}
        return tc_dic

    def flanking_up_composition(mrna, mr_site_loc, flank_number=70):  # 20
        mrna_full = mrna.upper().replace('T', 'U')
        mrna_up = mrna_full[max(0, mr_site_loc[0] - 70):mr_site_loc[0]]
        mrna_down = mrna_full[mr_site_loc[1] + 1:mr_site_loc[1] + 71]
        # print len(mrna_up), len(mrna_down)

        # # Up
        count_A = 0
        count_U = 0
        count_G = 0
        count_C = 0
        count_AA = 0
        count_AU = 0
        count_AG = 0
        count_AC = 0
        count_UA = 0
        count_UU = 0
        count_UG = 0
        count_UC = 0
        count_GA = 0
        count_GU = 0
        count_GG = 0
        count_GC = 0
        count_CA = 0
        count_CU = 0
        count_CG = 0
        count_CC = 0
        for i in range(len(mrna_up)):
            if mrna_up[i] == 'A':
                count_A += 1
            elif mrna_up[i] == 'U':
                count_U += 1
            elif mrna_up[i] == 'G':
                count_G += 1
            elif mrna_up[i] == 'C':
                count_C += 1
        for i in range(len(mrna_up) - 1):
            if mrna_up[i:i + 2] == 'AA':
                count_AA += 1
            elif mrna_up[i:i + 2] == 'AU':
                count_AU += 1
            elif mrna_up[i:i + 2] == 'AG':
                count_AG += 1
            elif mrna_up[i:i + 2] == 'AC':
                count_AC += 1
            elif mrna_up[i:i + 2] == 'UA':
                count_UA += 1
            elif mrna_up[i:i + 2] == 'UU':
                count_UU += 1
            elif mrna_up[i:i + 2] == 'UG':
                count_UG += 1
            elif mrna_up[i:i + 2] == 'UC':
                count_UC += 1
            elif mrna_up[i:i + 2] == 'GA':
                count_GA += 1
            elif mrna_up[i:i + 2] == 'GU':
                count_GU += 1
            elif mrna_up[i:i + 2] == 'GG':
                count_GG += 1
            elif mrna_up[i:i + 2] == 'GC':
                count_GC += 1
            elif mrna_up[i:i + 2] == 'CA':
                count_CA += 1
            elif mrna_up[i:i + 2] == 'CU':
                count_CU += 1
            elif mrna_up[i:i + 2] == 'CG':
                count_CG += 1
            elif mrna_up[i:i + 2] == 'CC':
                count_CC += 1
        all_monomer_count = count_A + count_U + count_G + count_C
        all_dimer_count = count_AA + count_AU + count_AG + count_AC + \
                          count_UA + count_UU + count_UG + count_UC + \
                          count_GA + count_GU + count_GG + count_GC + \
                          count_CA + count_CU + count_CG + count_CC
        if all_monomer_count == 0:
            all_monomer_count += 70
        if all_dimer_count == 0:
            all_dimer_count += 70
        fuc_dic = {'Up_A_comp': round(float(count_A) / all_monomer_count, 4),
                   'Up_U_comp': round(float(count_U) / all_monomer_count, 4),
                   'Up_G_comp': round(float(count_G) / all_monomer_count, 4),
                   'Up_C_comp': round(float(count_C) / all_monomer_count, 4),
                   'Up_AA_comp': round(float(count_AA) / all_dimer_count, 4),
                   'Up_AU_comp': round(float(count_AU) / all_dimer_count, 4),
                   'Up_AG_comp': round(float(count_AG) / all_dimer_count, 4),
                   'Up_AC_comp': round(float(count_AC) / all_dimer_count, 4),
                   'Up_UA_comp': round(float(count_UA) / all_dimer_count, 4),
                   'Up_UU_comp': round(float(count_UU) / all_dimer_count, 4),
                   'Up_UG_comp': round(float(count_UG) / all_dimer_count, 4),
                   'Up_UC_comp': round(float(count_UC) / all_dimer_count, 4),
                   'Up_GA_comp': round(float(count_GA) / all_dimer_count, 4),
                   'Up_GU_comp': round(float(count_GU) / all_dimer_count, 4),
                   'Up_GG_comp': round(float(count_GG) / all_dimer_count, 4),
                   'Up_GC_comp': round(float(count_GC) / all_dimer_count, 4),
                   'Up_CA_comp': round(float(count_CA) / all_dimer_count, 4),
                   'Up_CU_comp': round(float(count_CU) / all_dimer_count, 4),
                   'Up_CG_comp': round(float(count_CG) / all_dimer_count, 4),
                   'Up_CC_comp': round(float(count_CC) / all_dimer_count, 4)}
        return fuc_dic

    def flanking_down_composition(mrna, mr_site_loc, flank_number=70):  # 20
        mrna_full = mrna.upper().replace('T', 'U')
        mrna_up = mrna_full[max(0, mr_site_loc[0] - 70):mr_site_loc[0]]
        mrna_down = mrna_full[mr_site_loc[1] + 1:mr_site_loc[1] + 71]
        # print len(mrna_up), len(mrna_down)

        # # Down
        count_A = 0
        count_U = 0
        count_G = 0
        count_C = 0
        count_AA = 0
        count_AU = 0
        count_AG = 0
        count_AC = 0
        count_UA = 0
        count_UU = 0
        count_UG = 0
        count_UC = 0
        count_GA = 0
        count_GU = 0
        count_GG = 0
        count_GC = 0
        count_CA = 0
        count_CU = 0
        count_CG = 0
        count_CC = 0
        for i in range(len(mrna_down)):
            if mrna_down[i] == 'A':
                count_A += 1
            elif mrna_down[i] == 'U':
                count_U += 1
            elif mrna_down[i] == 'G':
                count_G += 1
            elif mrna_down[i] == 'C':
                count_C += 1
        for i in range(len(mrna_down) - 1):
            if mrna_down[i:i + 2] == 'AA':
                count_AA += 1
            elif mrna_down[i:i + 2] == 'AU':
                count_AU += 1
            elif mrna_down[i:i + 2] == 'AG':
                count_AG += 1
            elif mrna_down[i:i + 2] == 'AC':
                count_AC += 1
            elif mrna_down[i:i + 2] == 'UA':
                count_UA += 1
            elif mrna_down[i:i + 2] == 'UU':
                count_UU += 1
            elif mrna_down[i:i + 2] == 'UG':
                count_UG += 1
            elif mrna_down[i:i + 2] == 'UC':
                count_UC += 1
            elif mrna_down[i:i + 2] == 'GA':
                count_GA += 1
            elif mrna_down[i:i + 2] == 'GU':
                count_GU += 1
            elif mrna_down[i:i + 2] == 'GG':
                count_GG += 1
            elif mrna_down[i:i + 2] == 'GC':
                count_GC += 1
            elif mrna_down[i:i + 2] == 'CA':
                count_CA += 1
            elif mrna_down[i:i + 2] == 'CU':
                count_CU += 1
            elif mrna_down[i:i + 2] == 'CG':
                count_CG += 1
            elif mrna_down[i:i + 2] == 'CC':
                count_CC += 1
        all_monomer_count = count_A + count_U + count_G + count_C
        all_dimer_count = count_AA + count_AU + count_AG + count_AC + \
                          count_UA + count_UU + count_UG + count_UC + \
                          count_GA + count_GU + count_GG + count_GC + \
                          count_CA + count_CU + count_CG + count_CC
        if all_monomer_count == 0:
            all_monomer_count += 70
        if all_dimer_count == 0:
            all_dimer_count += 70
        fdc_dic = {'Down_A_comp': round(float(count_A) / all_monomer_count, 4),
                   'Down_U_comp': round(float(count_U) / all_monomer_count, 4),
                   'Down_G_comp': round(float(count_G) / all_monomer_count, 4),
                   'Down_C_comp': round(float(count_C) / all_monomer_count, 4),
                   'Down_AA_comp': round(float(count_AA) / all_dimer_count, 4),
                   'Down_AU_comp': round(float(count_AU) / all_dimer_count, 4),
                   'Down_AG_comp': round(float(count_AG) / all_dimer_count, 4),
                   'Down_AC_comp': round(float(count_AC) / all_dimer_count, 4),
                   'Down_UA_comp': round(float(count_UA) / all_dimer_count, 4),
                   'Down_UU_comp': round(float(count_UU) / all_dimer_count, 4),
                   'Down_UG_comp': round(float(count_UG) / all_dimer_count, 4),
                   'Down_UC_comp': round(float(count_UC) / all_dimer_count, 4),
                   'Down_GA_comp': round(float(count_GA) / all_dimer_count, 4),
                   'Down_GU_comp': round(float(count_GU) / all_dimer_count, 4),
                   'Down_GG_comp': round(float(count_GG) / all_dimer_count, 4),
                   'Down_GC_comp': round(float(count_GC) / all_dimer_count, 4),
                   'Down_CA_comp': round(float(count_CA) / all_dimer_count, 4),
                   'Down_CU_comp': round(float(count_CU) / all_dimer_count, 4),
                   'Down_CG_comp': round(float(count_CG) / all_dimer_count, 4),
                   'Down_CC_comp': round(float(count_CC) / all_dimer_count, 4)}
        return fdc_dic

    def tostring(self):
        # mmp = pd.DataFrame([self.mmp_dic])
        # mpc = pd.DataFrame([self.mpc_dic])
        # pd.set_option('display.max_columns', None)
        #
        # classstr = ""
        #
        # classstr = classstr + str(mmp) + "\n"
        # classstr = classstr + str(mpc) + "\n"
        #
        # return classstr


    def __str__(self):
        return self.tostring()
