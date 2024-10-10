import numpy as np
from scipy import signal
from numba import njit
import numba as nb

def to_arrays(seqs):
    res = np.zeros((len(seqs), len(seqs[0])), dtype="uint8")
    for i in range(res.shape[0]):
        res[i] = np.array(tuple(seqs[i]), dtype="|S1").view("uint8")
    return res

@nb.njit
def similarity(s1, s2, _seq_len):
    s = 0
    for i in range(s1.shape[0]):
        if (s1[i] == s2[i]) and ((s1[i] != 45) or (s2[i] != 45)):
            s += 1
        else:
            continue
    return s/_seq_len

@nb.njit(parallel=True)
def _pairwise_distances(res, seq_arrays1, seq_arrays2, seq_len1, seq_len2,seqs1, seqs2):
    for i in nb.prange(res.shape[0]):
        for j in range(res.shape[1]):
            _seq_len = min(seq_len1[i], seq_len2[j])
            res[i,j] = similarity(seq_arrays1[i], seq_arrays2[j], _seq_len)

def pairwise_distances(seqs1, seqs2, seq_len1, seq_len2):
    seq_arrays1 = to_arrays(seqs1)
    seq_arrays2 = to_arrays(seqs2)
    res = np.zeros((seq_arrays1.shape[0], seq_arrays2.shape[0]))
    _pairwise_distances(res, seq_arrays1, seq_arrays2, seq_len1, seq_len2, seqs1, seqs2)
    return res

def pairwise_distances_self(seqs1, seqs2):
    seq_arrays1 = to_arrays(seqs1)
    seq_arrays2 = to_arrays(seqs2)
    res = np.zeros((seq_arrays1.shape[0], seq_arrays2.shape[0]))
    _pairwise_distances(res, seq_arrays1, seq_arrays2)
    for i in nb.prange(res.shape[0]):
        res[i, i] = 0
    return res

def binary_occurrence(dom1, dom2):
    d1 = dom1.astype(bool)*1
    d2 = dom2.astype(bool)*1
    denom = max(sum(d1), sum(d2))
    return signal.correlate(d1, d2, mode="valid")[0]/denom

def most_frequent(List):
    return max(set(List), key = List.count)
