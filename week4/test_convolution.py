from week4.convolution import (
    find_peptides_with_convolution,
    list_convolution,
    most_common,
)


def test_convolution():
    assert sorted(list_convolution([0, 137, 186, 323])) == [49, 137, 137, 186, 186, 323]
    assert list_convolution([]) == []
    assert list_convolution([0, 2]) == [2]
    assert list_convolution([0]) == []


def _parse_convolution():
    with open("week4/data/spectral_convolution.txt") as fd:
        fd.readline()
        L = fd.readline().strip().split()
        spectrum = [int(i.strip()) for i in L]
        fd.readline()
        L = fd.readline().strip().split()
        convolution = [int(i.strip()) for i in L]
    return spectrum, convolution


def test_convolution_large():
    spectrum, exp_c = _parse_convolution()
    c = list_convolution(spectrum)
    assert sorted(c) == sorted(exp_c)


def test_most_common():
    d = {100: 3, 101: 2, 102: 2, 103: 1}
    assert most_common(d, 2) == [100, 101, 102]
    assert most_common(d, 100) == [100, 101, 102, 103]
    d = {100: 3, 101: 3, 102: 2, 103: 1}
    assert most_common(d, 2) == [100, 101, 102]
    assert most_common(d, 1) == [100, 101]
    d = {100: 3, 101: 2, 102: 1}
    assert most_common(d, 2) == [100, 101]


def test_findpeptide_with_convolution():
    # fmt: off
    sp = [57, 57, 71, 99, 129, 137, 170, 186, 194, 208, 228, 265, 285, 299, 307, 323, 356, 364, 394, 422, 493]
    # fmt: on
    M = 20
    N = 60
    rlb = find_peptides_with_convolution(M, N, sp)
    for r in rlb:
        if r == (99, 71, 137, 57, 72, 57):
            return
    raise Exception("No correct solution")


def _parse_convolution_large():
    with open("week4/data/convolution_cyclopeptide_sequencing.txt") as fd:
        fd.readline()
        M = int(fd.readline().strip())
        N = int(fd.readline().strip())
        L = fd.readline().strip().split()
        spectrum = [int(i.strip()) for i in L]
        fd.readline()
        L = fd.readline().strip().split("-")
        exp = [int(i.strip()) for i in L]
    return M, N, spectrum, exp


def test_find_peptide_with_convolution_large():
    M, N, spectrum, exp = _parse_convolution_large()
    p = find_peptides_with_convolution(M, N, spectrum)
    for pi in p:
        if pi[0] in exp:
            i = exp.index(pi[0])
            exp_mdf1 = exp[i:] + exp[:i]
            exp_mdf2 = exp[i:][::-1] + exp[:i][::-1]
            if list(pi) == exp_mdf1 or list(pi) == exp_mdf2:
                return
    raise Exception("No solution matched")
