from collections import Counter
from week3.cyclopeptide import (
    count_peptide_given_mass_rc,
    expand_peptides,
    expand_peptides_within_spectrum,
    find_all_subpeptides,
    find_cyclopeptide,
    find_cyclopeptide_spectrum,
    find_linear_subpeptides_with_len,
    find_possible_spectrum,
    find_subpeptides_with_len,
    _parse_peptide_mass,
    _is_a_in_b,
    generate_peptide_given_mass,
)


def test_findsubpeptides_withlen():
    p = "NQ"
    assert sorted(find_subpeptides_with_len(p, 1)) == ["N", "Q"]
    p = "QQ"
    assert sorted(find_subpeptides_with_len(p, 1)) == ["Q", "Q"]
    p = "ELEL"
    assert sorted(find_subpeptides_with_len(p, 2)) == sorted(["EL", "LE", "EL", "LE"])


def test_findlinearsubpeptides_withlen():
    assert find_linear_subpeptides_with_len("NQ", 1) == ["N", "Q"]
    assert find_linear_subpeptides_with_len("NQEL", 1) == ["N", "Q", "E", "L"]
    assert find_linear_subpeptides_with_len("NQEL", 2) == ["NQ", "QE", "EL"]
    assert find_linear_subpeptides_with_len("NQEL", 3) == ["NQE", "QEL"]
    assert find_linear_subpeptides_with_len("ELEL", 1) == ["E", "L", "E", "L"]
    assert find_linear_subpeptides_with_len("ELEL", 2) == ["EL", "LE", "EL"]


def test_find_allsubpeptides():
    p = "NQ"
    assert sorted(find_all_subpeptides(p)) == ["N", "Q"]
    assert sorted(find_all_subpeptides(p, True)) == ["N", "Q"]
    p = "QQ"
    assert sorted(find_all_subpeptides(p)) == ["Q", "Q"]
    assert sorted(find_all_subpeptides(p, True)) == ["Q", "Q"]
    p = "ELEL"
    assert sorted(find_all_subpeptides(p)) == sorted(
        ["E", "L", "E", "L", "EL", "LE", "EL", "LE", "ELE", "LEL", "LEL", "ELE"]
    )
    assert sorted(find_all_subpeptides("ELEL", True)) == sorted(
        ["E", "L", "E", "L", "EL", "LE", "EL", "ELE", "LEL"]
    )


def test_possible_spectrum():
    p = "NQ"
    assert find_possible_spectrum(p) == [0, 114, 128, 242]
    assert find_possible_spectrum(p, True) == [0, 114, 128, 242]
    p = "NQEL"
    p_exp = [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
    assert find_possible_spectrum(p) == p_exp
    p_exp = [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]
    assert find_possible_spectrum(p, linear=True) == p_exp


def _parse_spectrum_data():
    with open("week3/data/theoretical_spectrum.txt") as fd:
        fd.readline()
        p = fd.readline().strip()
        fd.readline()
        s = [int(i) for i in fd.readline().strip().split()]
    return p, s


def test_spectrum_large():
    peptide, exp_m = _parse_spectrum_data()
    m = find_possible_spectrum(peptide)
    assert m == exp_m


def test_countpeptides_mass_rc():
    d = [0]
    assert count_peptide_given_mass_rc(0, d) == 1
    d = [2, 5, 4, 7, 9]
    assert count_peptide_given_mass_rc(0, d) == 1
    assert count_peptide_given_mass_rc(3, d) == 0
    assert count_peptide_given_mass_rc(4, d) == 2
    assert count_peptide_given_mass_rc(5, d) == 1
    assert count_peptide_given_mass_rc(7, d) == 3
    assert count_peptide_given_mass_rc(9, d) == 8


def test_countpeptides_mass():
    # d = [0]
    # assert count_peptide_given_mass(0, d, []) == []
    # d = [2, 5, 4, 7, 9]
    # assert count_peptide_given_mass(0, d, []) == []
    # assert count_peptide_given_mass(3, d, []) == []
    # assert generate_peptide_given_mass(4, d) == [[2, 2], [4]]
    d = [2, 3]
    assert sorted(generate_peptide_given_mass(6, d)) == [[2, 2, 2], [3, 3]]


def test_is_ainb():
    assert _is_a_in_b([1], Counter([2, 3])) is False
    assert _is_a_in_b([1], Counter([1, 2, 3])) is True
    assert _is_a_in_b([], Counter([])) is True
    assert _is_a_in_b([1], Counter([1, 1, 2])) is True
    assert _is_a_in_b([1, 1], Counter([1, 2])) is False


def test_expand_peptide():
    assert list(expand_peptides([[2]], [2, 3])) == [(2, 2), (2, 3)]
    assert list(expand_peptides([[3, 1], [0]], [2])) == [(3, 1, 2), (0, 2)]
    assert list(expand_peptides([[1]], [2])) == [(1, 2)]


def test_cyclopep():
    s = [0, 113, 128, 186, 241, 299, 314, 427]
    exp = [
        (113, 128, 186),
        (113, 186, 128),
        (128, 113, 186),
        (128, 186, 113),
        (186, 113, 128),
        (186, 128, 113),
    ]
    assert set(find_cyclopeptide(s)) == set(exp)


def _parse_cycelope():
    with open("week3/data/cyclopeptide_sequencing.txt") as fd:
        fd.readline()
        spectrum = [int(i) for i in fd.readline().strip().split(" ")]
        fd.readline()
        P = [i.strip().split("-") for i in fd.readline().strip().split(" ")]
        peptides = []
        for p in P:
            peptides.append(tuple([int(i) for i in p]))
    return spectrum, peptides


def test_cyclop_large():
    s, p_exp = _parse_cycelope()
    p = find_cyclopeptide(s)
    assert sorted(p) == sorted(p_exp)


def test_find_cyclopeptide_spectrum():
    assert find_cyclopeptide_spectrum([1, 2, 3], 1) == [7, 6, 4]
    assert find_cyclopeptide_spectrum([], 1) == []
    assert find_cyclopeptide_spectrum([0], 3) == [3]


def test_expand_peptide_withspectrum():
    cnt_spectrum = Counter([0, 1, 2])
    assert list(expand_peptides_within_spectrum([[2]], [2], cnt_spectrum)) == []
    cnt_spectrum = Counter([0, 1, 2, 3])
    assert list(expand_peptides_within_spectrum([[2]], [1], cnt_spectrum)) == [(2, 1)]


# def test_find_linear_spectrum():
#     assert find_linear_spectrum("NQ") == [0, 114, 128, 242]
#     assert find_linear_spectrum("NQEL") ==

