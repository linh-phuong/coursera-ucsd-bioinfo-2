from week4.peptide_with_mismatch import (
    find_leading_peptide,
    generate_spectrum_from_masses,
    score_peptide,
    score_theoretical_spectrum,
    select_best_candidates,
)
from week3.cyclopeptide import _parse_peptide_mass
import pytest

MASS_DICT = _parse_peptide_mass()
MASSES = MASS_DICT.values()


@pytest.mark.parametrize(
    "peptides, masses, expected_count", [("NQ", [0, 1, 2], 1), ("NQL", MASSES, 3), ("QQ", [], 0)]
)
def test_score_peptide(peptides, masses, expected_count):
    assert score_peptide(peptides, masses) == expected_count


def test_score_theoretical_spectrum():
    assert score_theoretical_spectrum([], MASSES) == 0


def test_select_best_candidates():
    candidates = [[0], [1, 2], [2, 1]]
    spectrum = [0, 3]
    assert select_best_candidates(candidates, spectrum, 1) == [[1, 2], [2, 1]]
    assert select_best_candidates(candidates, spectrum, 2) == [[1, 2], [2, 1]]
    assert select_best_candidates(candidates, spectrum, 6) == [[0], [1, 2], [2, 1]]
    candidates = [[0], [1], [1, 2]]
    spectrum = [0, 1, 2]
    assert select_best_candidates(candidates, spectrum, 1) == [[1, 2]]
    assert select_best_candidates(candidates, spectrum, 2) == [[1], [1, 2]]


def mass_to_str(masses, mtsdict):
    p_str = "".join([mtsdict.get(i, "*") for i in masses])
    return p_str


def str_to_mass(peptide_str, stmdict):
    return [stmdict[i] for i in peptide_str]


def _parse_best_candidate():
    with open("week4/data/trim.txt") as fd:
        fd.readline()
        P = fd.readline().strip().split()
        peptides_str = [p.strip() for p in P]
        L = fd.readline().strip().split()
        spectrum = [int(i.strip()) for i in L]
        N = int(fd.readline().strip())
        fd.readline()
        L = fd.readline().strip().split()
        exp = [l.strip() for l in L]
    return peptides_str, spectrum, N, exp


def test_best_candidate_large():
    p_to_m = _parse_peptide_mass()
    pstr, spectrum, N, exp = _parse_best_candidate()
    pm = [str_to_mass(i, p_to_m) for i in pstr]
    rm = select_best_candidates(pm, spectrum, N)
    exp_m = [str_to_mass(i, p_to_m) for i in exp]
    assert sorted(rm) == sorted(exp_m)


def test_best_candidates_with_str():
    p_to_m = _parse_peptide_mass()
    m_to_p = {v: k for k, v in p_to_m.items()}
    p_str = ["LAST", "ALST", "TLLT", "TQAS"]
    p_m = [str_to_mass(i, p_to_m) for i in p_str]
    spectrum = [0, 71, 87, 101, 113, 158, 184, 188, 259, 271, 372]
    best_m = select_best_candidates(p_m, spectrum, 2)
    best_str = [mass_to_str(m, m_to_p) for m in best_m]
    assert best_str == ["LAST", "ALST"]


def test_generate_spectrum_from_masses():
    p = [114, 128]
    assert generate_spectrum_from_masses(p) == [0, 114, 128, 242]
    assert generate_spectrum_from_masses(p, True) == [0, 114, 128, 242]
    p = [114, 128, 129, 113]
    p_exp = [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]
    assert generate_spectrum_from_masses(p) == p_exp
    p_exp = [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]
    assert generate_spectrum_from_masses(p, linear=True) == p_exp


def test_find_leading():
    spectrum = [0, 71, 113, 129, 147, 200, 218, 260, 313, 331, 347, 389, 460]
    exp = (113, 147, 71, 129)
    r = find_leading_peptide(spectrum, 10)
    i = r.index(exp[0])
    r = r[i:] + r[:i]
    assert exp == r


def _parse_leaderboard_large():
    with open("week4/data/leaderboard_cyclopeptide_sequencing.txt") as fd:
        fd.readline()
        N = int(fd.readline().strip())
        L = fd.readline().strip().split()
        spectrum = [int(i.strip()) for i in L]
        fd.readline()
        out = fd.readline().strip().split("-")
        exp_lb = [int(i.strip()) for i in out]
    return spectrum, N, exp_lb


def test_find_leading_large():
    spectrum, N, _ = _parse_leaderboard_large()
    results = find_leading_peptide(spectrum, N)
    # fmt: off
    exp_lb = (71, 99, 147, 97, 129, 97, 114, 163, 137, 101, 128, 87, 103, 128, 113, 71, 115, 163, 113, 71, 186)
    # fmt: on
    assert results == exp_lb
