# cython: language_level=3

import cython

def group_insilico_amplicons(list in_silico_amplicon_list):
    cdef int i, lent
    cdef list amplicon_headers, bgc_set
    cdef dict bgc_dict

    amplicon_headers = [x[0] for x in in_silico_amplicon_list]
    bgc_set = list(set(['_'.join(x.split('_')[:-1]) for x in amplicon_headers]))
    bgc_dict = {}
    i = 1
    lent = len(bgc_set)
    for bgc in bgc_set:
        bgc_dict[bgc] = [x for x in amplicon_headers if x.startswith(bgc)]
        print(f"\r{i} / {lent}", end='')
        i += 1
    return bgc_dict