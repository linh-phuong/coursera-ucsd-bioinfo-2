from week3.protein_translate import (
    find_all_peptides_encoding,
    find_peptide_encoding,
    translate_rna,
    find_reverse_complement,
    translate_rna_with_position,
)
import pytest


def test_protein_translate_small():
    rna = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
    exp_pr = "MAMAPRTEINSTRING"
    pr = translate_rna(rna)
    assert exp_pr == pr


def test_transcript_rna_withposition():
    rna = "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA"
    exp_pr = "MAMAPRTEINSTRING*"
    pr, po = translate_rna_with_position(rna)
    assert pr == exp_pr
    assert po == [i for i in range(0, len(rna), 3)]


def _parse_translate_large():
    with open("week3/data/protein_translation.txt") as fd:
        fd.readline()
        rna = fd.readline().strip()
        fd.readline()
        pr = fd.readline().strip()
    return rna, pr


def test_translation_large():
    rna, exp_pr = _parse_translate_large()
    pr = translate_rna(rna)
    assert exp_pr == pr


def test_find_reverse():
    rna = "ACGT"
    assert find_reverse_complement(rna) == "TGCA"

    rna = "A"
    assert find_reverse_complement(rna) == "T"

    rna = "CCC"
    assert find_reverse_complement(rna) == "GGG"


def test_findpeptide():
    ncoding = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
    ncoding = ncoding.replace("T", "U")
    exp_substr = ["AUGGCC", "AUGGCC"]  # "GGCCAT",
    substr = find_peptide_encoding(ncoding, "MA")
    assert sorted(exp_substr) == sorted(substr)


def test_find_all_encodint():
    ncoding = "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
    exp_substr = ["ATGGCC", "ATGGCC", "GGCCAT"]
    substr = find_all_peptides_encoding(ncoding, "MA")
    assert sorted(exp_substr) == sorted(substr)


def _parse_peptides_large(file):
    with open(file) as fd:
        fd.readline()
        noncoding = fd.readline().strip()
        proteins = fd.readline().strip()
        fd.readline()
        substring = []
        while fd.readline():
            substring.append(fd.readline().strip())
    return noncoding, proteins, substring


@pytest.mark.xfail(
    reason="multiple results can be accepted while the test provided by the excercise only gives one"
)
def test_findpeptide_large():
    noncoding, proteins, exp_substring = _parse_peptides_large("week3/data/peptide_encoding.txt")
    substrings = find_all_peptides_encoding(noncoding[0:], proteins)
    assert sorted(substrings) == sorted(exp_substring)
