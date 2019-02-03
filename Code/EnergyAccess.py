import RNA
import os
import pandas as pd
import re

class EnergyAccess(object):

    def __init__(self, dp, site, site_start_loc, site_end_loc, full_mrna):
        self.dp = dp
        self.site = site
        self.site_start_loc = site_start_loc
        self.site_end_loc = site_end_loc
        self.full_mrna = full_mrna

        self.extract_energy_()

    def extract_energy_(self):
        mrna = self.full_mrna
        mr_site_loc = [self.site_start_loc, self.site_end_loc]

        self.MFE = self.minimum_free_energy(self.dp, mrna, self.site, mr_site_loc)
        self.ACC = self.accessibility(mrna, mr_site_loc)


    def get_features(self):
        f = [self.MFE, self.ACC]
        df = map(lambda x: pd.DataFrame([x]), f)
        r =  reduce (lambda x,y: pd.concat([x, y], axis=1, sort=False), df)
        return r


    # # 5. energy  (5)
    def minimum_free_energy(self, dp, mrna, mrna_site, mr_site_loc):  # 5
        mrna_full = mrna.upper().replace('T', 'U')
        mrna_surrounding100 = mrna_full[max(0, mr_site_loc[0] - 50):mr_site_loc[1] + 50]
        seed = dp.IRP.extract_seed()
        mrna_seed = seed.site
        mrna_site = mrna_site[::-1].upper().replace('T', 'U')
        mrna_site_3p = mrna_site[len(mrna_seed):]

        MFE_Seed = RNA.fold(mrna_seed)
        MFE_surrounding100 = RNA.fold(mrna_surrounding100)
        MFE_surrounding100_norm = MFE_surrounding100[1] / len(mrna_surrounding100)
        #MEF_duplex = dp.duplex.energy #skip it to support vienna
        duplex = RNA.duplexfold(dp.mir, mrna_site[::-1])
        MEF_duplex = duplex.energy

        MFE_3p = RNA.fold(mrna_site_3p)

        constraint_low = "."*min(mr_site_loc[0],50)
        constraint_site = "x"*len(mrna_site)
        constraint_high = "."*min(len(mrna) - mr_site_loc[1], 50)
        constraint = constraint_low + constraint_site + constraint_high
        assert len(constraint)==len(mrna_surrounding100), "constraint and mrna_surrounding100 are not in the same length"

        mrna100file_in = 'mrna100_with_constraints{}.fa'.format(os.getpid())
        mrna100file_out = 'mrna100_with_constraints{}.result'.format(os.getpid())
        f = open(mrna100file_in, 'w')
        f.write(mrna_surrounding100 + "\n" + constraint + "\n")
        f.close()
        cmd = "RNAfold -C {infile} > {outfile}".format(infile=mrna100file_in, outfile=mrna100file_out)
        os.system(cmd)
        f = open(mrna100file_out, 'r')
        twolines = f.readlines()
        f.close()
        cmfe =float(re.findall("\d+\.\d+", twolines[1])[0])*(-1)
        # (struct, cmfe) = RNA.fold(mrna_surrounding100, constraint)
        os.remove(mrna100file_in)
        os.remove(mrna100file_out)

        MFE = {'MEF_Seed': round(MFE_Seed[1], 4),
               'MEF_3p': round(MFE_3p[1], 4),
               'MEF_local_target': round(MFE_surrounding100[1], 4),
               'MEF_local_target_normalized': round(MFE_surrounding100_norm, 4),
               'MEF_Duplex': round(MEF_duplex, 4),
               'MEF_cons_local_target': round(cmfe, 4),
               'MEF_cons_local_target_normalized': round(cmfe/len(mrna_surrounding100), 4)
               }
        return MFE


    # # 6. target site accessibility  (370)
    def accessibility(self, mrna, mr_site_loc):  # 37*10 = 370
        acc_file ="mrna_acc{}.fa".format(os.getpid())

        f = open(acc_file, 'w')
        f.write(mrna.replace('-', '') + '\n')
        f.close()

        os.system('RNAplfold -W 80 -L 40 -u 10 < {}'.format(acc_file))
        f = open('plfold_lunp','r')
        ACC_allstr = f.readlines()
        f.close()

        acc_score_matrix = []
        for line in ACC_allstr[2:]:
            l = line.strip().replace('NA', str(0.97)).split()
            acc_score_matrix.append(l)
        acc_score_matrix = [[0.0] * 11] * 37 + acc_score_matrix + [[0.0] * 11] * 37
        acc_score_matrix_segment = acc_score_matrix[mr_site_loc[1]+15:mr_site_loc[1]+52]
        # print mr_site_loc
        ACC = {}
        # print len(acc_score_matrix)
        # print acc_score_matrix_segment
        for i in range(1, 38):
            # print acc_score_matrix_segment[i-1]
            for j in range(1, 11):
                key = 'Acc_P%s_%sth' % (str(i), str(j))
                ACC[key] = float(acc_score_matrix_segment[i - 1][j])
        os.remove(acc_file)
        os.remove('plfold_lunp')
        os.remove('plfold_dp.ps')
        return ACC
