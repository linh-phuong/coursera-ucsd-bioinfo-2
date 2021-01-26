from collections import defaultdict
from week4.peptide_with_mismatch import _find_leaderboard_with_masseslist


def dict_convolution(exp_spectrum):
    """find spectral convolution

    Args:
        exp_spectrum (list): spectrum obtained from experiments

    Returns:
        dict: spectral convolution
    """
    exp_spectrum = sorted(exp_spectrum)
    c_spectrum = defaultdict(lambda: 0)
    for s0 in exp_spectrum[1:]:
        for s1 in exp_spectrum:
            cs = s0 - s1
            if cs > 0:
                c_spectrum[cs] += 1
            else:
                break
    return c_spectrum


def dict_convolution_filter(exp_spectrum):
    exp_spectrum = sorted(exp_spectrum)
    c_spectrum = defaultdict(lambda: 0)
    for s0 in exp_spectrum[1:]:
        for s1 in exp_spectrum:
            cs = s0 - s1
            if 200 >= cs >= 57:
                c_spectrum[cs] += 1
            elif cs < 0:
                break
    return c_spectrum


def list_convolution(exp_spectrum):
    """return spectral convolution as a list

    """
    dc = dict_convolution(exp_spectrum)
    ls = []
    for i in dc:
        ls += [i] * dc[i]
    return ls


def most_common(convo_dict_filter, M):
    assert M >= 1
    ls = []
    sorted_d = sorted(convo_dict_filter.items(), key=lambda x: x[1], reverse=True)
    val_sorted = sorted(set(convo_dict_filter.values()), reverse=True)
    t = val_sorted[M - 1] if len(val_sorted) >= M else val_sorted[-1]
    for i, j in sorted_d:
        if j >= t:
            ls.append(i)
    return ls


def find_peptides_with_convolution(M, N, spectrum):
    if 0 not in spectrum:
        spectrum.append(0)
    convolution = dict_convolution_filter(spectrum)
    common_c = most_common(convolution, M)
    lb = _find_leaderboard_with_masseslist(spectrum, N, common_c)
    return lb
