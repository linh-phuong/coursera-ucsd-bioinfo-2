from week2.contigs import ContigGeneration, MaxBranch, get_max_out_count, in_out_count
from pathlib import Path
import pytest
import textwrap


def _parse_reads(file):
    with open(file) as fd:
        reads = [i.strip() for i in fd.readlines()]
        return reads


def _parse_contigs(file):
    with open(file) as fd:
        contigs = [i.strip() for i in fd.readline().split()]
        return contigs


def _parse_contigs_large(file):
    with open(file) as fd:
        fd.readline()
        reads = []
        contigs = []
        while True:
            line = fd.readline()
            if line.strip() != "Output:":
                reads.append(line.strip())
            else:
                break
        while True:
            line = fd.readline()
            if not line:
                break
            else:
                contigs.append(line.strip())
    return reads, contigs


def test_parsedata(tmp_path):
    fn = tmp_path / "a.txt"
    fn.write_text(
        textwrap.dedent(
            """
    Input
    ABCD
    DD
    Output:
    ABBBB
    MDD
    DD
    """
        ).strip()
    )

    assert (["ABCD", "DD"], ["ABBBB", "MDD", "DD"]) == _parse_contigs_large(fn)


def test_contig_indv():
    inp = "week2/data/Contigs/inputs/test3.txt"
    out = str(inp).replace("/inputs/", "/outputs/")
    reads = _parse_reads(inp)
    exp_c = _parse_contigs(out)
    c = ContigGeneration(reads)
    assert sorted(c) == sorted(exp_c)


TEST_DIR = Path("week2/data/Contigs/inputs")


@pytest.mark.parametrize("inp", sorted(list(TEST_DIR.glob("*"))))
def test_ContigGenerate(inp):
    out = str(inp).replace("/inputs/", "/outputs/")
    reads = _parse_reads(inp)
    c = ContigGeneration(reads)
    exp_c = _parse_contigs(out)
    assert sorted(c) == sorted(exp_c)


def test_get_max_out_count():
    assert "A" == get_max_out_count({"A": (2, 1), "B": (1, 2)})
    assert "B" == get_max_out_count({"A": (0, 1), "B": (1, 2)})


def test_in_out_count():
    assert {"A": [1, 0], "B": [0, 1]} == in_out_count({"A": ["B"]})
    assert {"A": [1, 1], "B": [0, 2], "C": [2, 0]} == in_out_count({"A": ["B"], "C": ["B", "A"]})


def test_contig_large():
    reads, exp_cotg = _parse_contigs_large("week2/data/contig_generation.txt")
    cotg = ContigGeneration(reads)
    assert sorted(cotg) == sorted(exp_cotg)


def _parse_dB_ind(file):
    with open(file) as fd:
        lines = fd.readlines()
        G = {}
        for l in lines:
            kv = l.strip().split("->")
            G[kv[0].strip()] = kv[1].strip().split(",")
    return G


def _parse_non_branch(file):
    with open(file) as fd:
        lines = fd.readlines()
        contigs = [l.strip() for l in lines]
    return contigs


def test_parsemaxbranch(tmp_path):
    fn = tmp_path / "a.txt"
    fn.write_text(
        textwrap.dedent(
            """
    Input
    0 -> 1,2
    2 -> 4
    Output
    0 -> 2 -> 1
    4 -> 4
    """
        ).strip()
    )

    assert ({"0": ["1", "2"], "2": ["4"]}, ["0->2->1", "4->4"]) == _parse_maxbranch_large(fn)


def test_max_branch_indv():
    inp = "week2/data/MaximalNonBranchingPaths/inputs/sample.txt"
    out = str(inp).replace("/inputs/", "/outputs/")
    G = _parse_dB_ind(inp)
    exp_nb = _parse_non_branch(out)
    nb = MaxBranch(G)
    assert (sorted(exp_nb)) == sorted(nb)


TEST_DIR = Path("week2/data/MaximalNonBranchingPaths/inputs")


@pytest.mark.parametrize("inp", sorted(list(TEST_DIR.glob("*"))))
def test_MaxBranch(inp):
    out = str(inp).replace("/inputs/", "/outputs/")
    G = _parse_dB_ind(inp)
    G_clone = {k: v for k, v in G.items()}
    c = MaxBranch(G_clone)
    exp_c = _parse_non_branch(out)
    assert sorted(c) == sorted(exp_c)


def _parse_maxbranch_large(file):
    with open(file) as fd:
        fd.readline()
        G = {}
        nb = []
        while True:
            line = fd.readline()
            if line.strip() != "Output":
                d = line.split("->")
                G[d[0].strip()] = d[1].strip().split(",")
            else:
                break
        while True:
            line = fd.readline()
            if not line:
                break
            else:
                nb.append(line.strip().replace(" ", ""))
        return G, nb


def test_maxbranch_large():
    G, exp_nb = _parse_maxbranch_large("week2/data/MaximalNonBranchingPaths.txt")
    nb = MaxBranch(G)
    assert sorted(exp_nb) == sorted(nb)
