from week1.debruijn import DeBruijnGraphFromReads, StringComposition
from week2.euler import (
    EulerianPath,
    IsEulerianPath,
    IsEulerianPathNb,
    PathOrdered,
    RandomCycle,
    EulerianCycle,
    IsEulerianCycle,
    StringReconstruction,
    UnbalancedNodes,
    UniversalCircular,
)
from pathlib import Path
import pytest
from itertools import product


def test_RandomCycle():
    test = {0: [1], 1: [2], 2: [0]}
    stest1 = RandomCycle(test)
    assert stest1 == [[1, 2, 0, 1]]
    test = {0: [3, 1], 1: [2], 2: [0], 3: [0]}
    stest2 = RandomCycle(test)
    assert stest2 == [[3, 0, 3], [1, 2, 0, 1]]
    test = {1: [2], 2: [1, 2]}
    stest3 = RandomCycle(test)
    assert stest3 == [[2, 1, 2], [2, 2]]


@pytest.mark.parametrize(
    "path,start,expected",
    [
        ([0, 5, 1, 2, 0], 2, [2, 0, 5, 1, 2]),
        ([0, 5, 1, 2, 0], 0, [0, 5, 1, 2, 0]),
        ([1, 1, 1], 1, [1, 1, 1]),
        ([1, 3, 2, 3, 1], 3, [3, 2, 3, 1, 3]),
        ([4], 4, [4]),
    ],
)
def test_PathOrdered(path, start, expected):
    assert expected == PathOrdered(start, path)


TEST_DIR = Path("week2/data/EulerianCycle/inputs/")


@pytest.mark.parametrize(
    "inp", list(TEST_DIR.glob("*")),
)
def test_run(inp):
    G = _parse_graph(inp)
    cycle = EulerianCycle(G)
    assert IsEulerianCycle(cycle, G)


def _parse_graph(path):
    with open(path) as fd:
        ret = {}
        for s in fd.readlines():
            k, vs = s.strip().split("->")
            ret[k.strip()] = vs.strip().split(",")
        return ret


def test_UnbalancedNode():
    test = {0: [2], 2: [2, 3], 4: [2], 3: [4]}
    assert UnbalancedNodes(test) == (2, 0)
    test = {0: [2], 1: [3], 2: [1], 3: [0, 4], 6: [3, 7], 7: [8], 8: [9], 9: [6]}
    assert UnbalancedNodes(test) == (4, 6)
    test = {0: [1], 1: [2], 2: [3]}
    assert UnbalancedNodes(test) == (3, 0)


def test_EulerianPath():
    test = {0: [2], 1: [3], 2: [1], 3: [0, 4], 6: [3, 7], 7: [8], 8: [9], 9: [6]}
    assert EulerianPath(test) == [6, 7, 8, 9, 6, 3, 0, 2, 1, 3, 4]
    assert IsEulerianPathNb([6, 7, 8, 9, 6, 3, 0, 2, 1, 3, 4], test)

    test = {0: [1], 1: [2], 2: [3]}
    assert EulerianPath(test) == [0, 1, 2, 3]
    assert IsEulerianPathNb([0, 1, 2, 3], test)

    test = {0: [1], 1: [2, 5], 2: [3], 3: [4], 4: [1]}
    assert EulerianPath(test) == [0, 1, 2, 3, 4, 1, 5]
    assert IsEulerianPathNb([0, 1, 2, 3, 4, 1, 5], test)

    test = {2: [1], 1: [3, 4, 0], 3: [1, 4], 4: [3, 1]}
    assert EulerianPath(test) == [2, 1, 3, 4, 3, 1, 4, 1, 0]
    assert IsEulerianPathNb([2, 1, 3, 4, 3, 1, 4, 1, 0], test)

    test = {0: [1], 1: [14, 17], 14: [2, 3, 4], 2: [1], 3: [14], 4: [5], 5: [14]}
    assert EulerianPath(test) == [0, 1, 14, 3, 14, 4, 5, 14, 2, 1, 17]
    assert IsEulerianPathNb([0, 1, 14, 3, 14, 4, 5, 14, 2, 1, 17], test)

    test = {2: [3, 5], 3: [4], 4: [2], 5: [6], 6: [2], 1: [2, 0], 0: [1]}
    assert EulerianPath(test) == [1, 0, 1, 2, 5, 6, 2, 3, 4, 2]
    assert IsEulerianPathNb([1, 0, 1, 2, 5, 6, 2, 3, 4, 2], test)


TEST_DIR = Path("week2/data/StringReconstruction/inputs/")


@pytest.mark.parametrize("inp", list(TEST_DIR.glob("*")))
def test_StringReconstruction(inp):
    out = str(inp).replace("/inputs/", "/outputs/")
    k, reads = _parse_reads_and_k(inp)
    exp_gene = _parse_gene(out)
    gene = StringReconstruction(reads)
    G = DeBruijnGraphFromReads(reads)
    assert IsEulerianPath(gene, G)
    assert IsEulerianPath(exp_gene, G)


def test_StringReconstruction1():
    inp = "week2/data/StringReconstruction/inputs/test5.txt"
    out = str(inp).replace("/inputs/", "/outputs/")
    k, reads = _parse_reads_and_k(inp)
    exp_gene = _parse_gene(out)
    gene = StringReconstruction(reads)
    G = DeBruijnGraphFromReads(reads)
    assert IsEulerianPath(gene, G)
    assert IsEulerianPath(exp_gene, G)


def _parse_reads_and_k(file):
    with open(file) as fd:
        k = fd.readline().strip()
        reads = [s.strip() for s in fd.readlines()]
    return int(k), reads


def _parse_gene(file):
    with open(file) as fd:
        return fd.readline().strip()


def test_StringReconstruct_large():
    exp_gene = None
    with open("week2/data/StringReconstructionProblem.txt") as fd:
        fd.readline()
        fd.readline()
        reads = []
        for read in fd.readlines():
            if read.strip() != "Output:":
                reads.append(read.strip())
            else:
                break
            exp_gene = fd.readline()
    gene = StringReconstruction(reads)
    G = DeBruijnGraphFromReads(reads)
    assert IsEulerianPath(gene, G)
    assert IsEulerianPath(exp_gene, G)


def test_UniversalCircle():
    assert UniversalCircular(2) == "0011"
    k = 3
    c = UniversalCircular(k)
    comp = sorted(["".join(i) for i in product("01", repeat=k)])
    cs = StringComposition(k, c + c[0 : k - 1])  # noqa: E203
    assert sorted(cs) == sorted(comp)
    k = 4
    c = UniversalCircular(k)
    comp = sorted(["".join(i) for i in product("01", repeat=k)])
    cs = StringComposition(k, c + c[0 : k - 1])  # noqa: E203
    assert sorted(cs) == sorted(comp)


def test_UniversalCircle_large():
    with open("week2/data/universal_string.txt") as fd:
        fd.readline()
        k = int(fd.readline().strip())
        fd.readline()
        s = fd.readline().strip()
    ms = UniversalCircular(k)
    exp_c = s + s[0 : k - 1]
    c = ms + ms[0 : k - 1]
    exp_comp = StringComposition(k, exp_c)
    c_comp = StringComposition(k, c)
    correct_comp = sorted(["".join(i) for i in product("01", repeat=k)])
    assert sorted(exp_comp) == correct_comp
    assert sorted(c_comp) == correct_comp
