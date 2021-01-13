from week2.euler_readpair import StringReconstructFromPairPath
from pathlib import Path
import pytest
import textwrap

TEST_DIR = Path("week2/data/PairedStringReconstruction/inputs")


@pytest.mark.parametrize("inp", sorted(list(TEST_DIR.glob("*"))))
def test_StringReconstuctPair(inp):
    out = str(inp).replace("/inputs/", "/outputs/")
    k, d, pairs = _parse_reads_pairs(inp)
    gene = StringReconstructFromPairPath(pairs, k, d)
    exp_gene = _parse_gene(out)
    assert gene == exp_gene


def _parse_reads_pairs(file):
    with open(file) as fd:
        nb = fd.readline().strip().split()
        k, d = [int(x.strip()) for x in nb]
        pairs = []
        for l in fd.readlines():
            pairs.append(l.strip().split("|"))
    return k, d, pairs


def _parse_gene(file):
    with open(file) as fd:
        return fd.readline().strip()


def test_StringReconstuctPair_ind():
    inp = "week2/data/PairedStringReconstruction/inputs/test3.txt"
    out = str(inp).replace("/inputs/", "/outputs/")
    k, d, pairs = _parse_reads_pairs(inp)
    gene = StringReconstructFromPairPath(pairs, k, d)
    exp_gene = _parse_gene(out)
    assert gene == exp_gene


def test_StringReconstructPair_large():
    file = "week2/data/StringReconstructionFromReadPairs.txt"
    k, d, exp_g, readpairs = _parse_StringReconstructPair_large(file)
    g = StringReconstructFromPairPath(readpairs, k, d)
    assert exp_g == g


def test_parselargedata(tmp_path):
    fn = tmp_path / "a.txt"
    fn.write_text(
        textwrap.dedent(
            """
    Input
    123 21
    hello|world
    hi|foo
    Output
    bar
    """
        ).strip()
    )

    assert (
        123,
        21,
        "bar",
        [["hello", "world"], ["hi", "foo"]],
    ) == _parse_StringReconstructPair_large(fn)


def _parse_StringReconstructPair_large(file):
    with open(file) as fd:
        fd.readline()
        line = fd.readline().strip()
        k, d = [int(x.strip()) for x in line.split()]

        readpairs = []
        while True:
            line = fd.readline().strip()
            if line == "Output":
                break
            else:
                readpairs.append(line.strip().split("|"))

        exp_g = fd.readline().strip()
        print(f"output is {exp_g}")
        return k, d, exp_g, readpairs
