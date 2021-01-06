from week1.debruijn import (
    DeBruijnGraph,
    DeBruijnGraphFromReads,
    OverlapGraph,
    PathToGenome,
    StringComposition,
)
from pathlib import Path
import pytest


def _parse_gene(file):
    with open(file) as fd:
        k = int(fd.readline())
        s = fd.readline()
        return k, s


def _parse_compositions(file):
    out = []
    with open(file) as output_file:
        lines_output = output_file.readlines()
        for l in lines_output:
            out.append(l.strip())
    return out


TEST_DIR = Path("week1/data/kmerComposition/inputs")


@pytest.mark.parametrize(
    "inp", list(TEST_DIR.glob("*")),
)
def test_StringComposition(inp):
    out = str(inp).replace("/inputs/", "/outputs/")
    k, string = _parse_gene(inp)
    c = StringComposition(k, string)
    exp_c = _parse_compositions(out)
    assert len(c) == len(exp_c)
    for i in c:
        assert i in exp_c, f"{i} is not correct in {inp}"


def _parse_reads(file):
    with open(file) as fd:
        return [s.strip() for s in fd.readlines()]


def _parse_gene_from_reads(file):
    with open(file) as fd:
        return fd.readline().strip()


TEST_DIR = Path("week1/data/GenomePath/inputs")


@pytest.mark.parametrize(
    "inp", list(TEST_DIR.glob("*")),
)
def test_GenomePath(inp):
    out = str(inp).replace("/inputs/", "/outputs/")
    reads = _parse_reads(inp)
    g = PathToGenome(reads)
    exp_g = _parse_gene_from_reads(out)
    g == exp_g


def _parse_graph(path):
    with open(path) as fd:
        ret = {}
        for s in fd.readlines():
            k, vs = s.strip().split("->")
            ret[k.strip()] = vs.strip().split(",")
        return ret


TEST_DIR = Path("week1/data/OverlapGraph/inputs")


@pytest.mark.parametrize(
    "inp", list(TEST_DIR.glob("*")),
)
def test_OverlapGraph(inp):
    out = str(inp).replace("/inputs/", "/outputs/")
    reads = _parse_reads(inp)
    G = OverlapGraph(reads)
    exp_G = _parse_graph(out)
    for k in G:
        assert sorted(G[k]) == sorted(exp_G[k])


TEST_DIR = Path("week1/data/deBruijnGraphString/inputs")


@pytest.mark.parametrize(
    "inp", list(TEST_DIR.glob("*")),
)
def test_DeBruijnGraph(inp):
    out = str(inp).replace("/inputs/", "/outputs/")
    k, reads = _parse_gene(inp)
    G = DeBruijnGraph(k, reads)
    exp_G = _parse_graph(out)
    for k in G:
        assert sorted(G[k]) == sorted(exp_G[k])


TEST_DIR = Path("week1/data/deBruijnGraphPatterns/inputs")


@pytest.mark.parametrize(
    "inp", list(TEST_DIR.glob("*")),
)
def test_DeBruijnGraphReads(inp):
    out = str(inp).replace("/inputs/", "/outputs/")
    reads = _parse_reads(inp)
    G = DeBruijnGraphFromReads(reads)
    exp_G = _parse_graph(out)
    for k in G:
        assert sorted(G[k]) == sorted(exp_G[k])

