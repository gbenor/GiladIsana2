import re
from scipy.spatial import distance
import numpy as np


class TooManyMatches(Exception):
    pass


def extract_site (site_with_stars_and_hashtags, mRNA_seq_extended, mRNA_seq_extended_offset, full_mrna_seq, ensg, vienna_hint=None):
    site_stars = site_with_stars_and_hashtags.count("*")
    site_hashtags = site_with_stars_and_hashtags.count("#")
    site = site_with_stars_and_hashtags.replace("*", "")
    site = site.replace("#", "")
    site = site[::-1]

    one_match=True
    start_end_list = [(m.start(0), m.end(0)) for m in re.finditer(site, mRNA_seq_extended)]
    if len(start_end_list) == 0:
        raise Exception(
            "no match for the site within the full mRNA seq. \nENSG={}\nsite={}\nmrna={}".format(ensg, site,
                                                                                                 full_mrna_seq))
    if len(start_end_list) > 1:
        if not vienna_hint:
            raise TooManyMatches(
                "too many matches for the site within the full mRNA seq. \nENSG={}\nsite={}\nlist={}\nvienna_hint={}".format(ensg, site,
                                                                                                        start_end_list, vienna_hint))

        dist = [distance.euclidean(vienna_hint, x) for x in start_end_list]
        assert min(dist)<50, "I think i found incorrect site place\nENSG={}\nsite={}\nlist={}\nvienna_hint={}".format(ensg, site,
                                                                                                        start_end_list, vienna_hint)

        site_start_loc, site_end_loc = start_end_list[np.argmin(dist)]
        one_match = False
        # for p in start_end_list:
        #     s, e = p
        #     sv, ev = vienna_hint
        #     if s>= sv and e<=ev:
        #         site_start_loc, site_end_loc = p
        #         one_match = False
        # if one_match:
        #     raise TooManyMatches(
        #         "try to handle with vienna_hint. unfortuntly, too many matches for the site within the full mRNA seq. \nENSG={}\nsite={}\nlist={}\nvienna_hint={}".format(ensg, site,
        #                                                                                                  start_end_list, vienna_hint))


    if one_match:
        site_start_loc, site_end_loc = start_end_list[0]
    site_start_loc = max(0, site_start_loc + int(mRNA_seq_extended_offset) - site_stars)
    site_end_loc =  max(0, site_end_loc + int(mRNA_seq_extended_offset) + site_hashtags)

    return full_mrna_seq[site_start_loc:site_end_loc], site_start_loc, site_end_loc

