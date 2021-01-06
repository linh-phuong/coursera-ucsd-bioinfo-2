from week2.euler import PathOrdered, RandomCycle, EulerianCycle, IsEulerianCycle
from pathlib import Path
import pytest


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
