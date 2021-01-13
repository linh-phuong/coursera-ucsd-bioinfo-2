from week2.euler import EulerianPath


# def StringReconstructFromPairPath(pairs, k, d):
#     pairs_graph = PairCompositionGraph(pairs)
#     PairPath = EulerianPath(pairs_graph)
#     s0 = PairPath[0][0]
#     s1 = PairPath[0][1]
#     for p in PairPath[1:]:
#         s0 += p[0][-1]
#         s1 += p[1][-1]
#     s = len(s1) - (d + k)
#     nonoverlap = s1[s:]
#     return s0 + nonoverlap


def StringReconstructFromPairPath(pairs, k, d):
    overlaps0 = 0
    overlaps1 = 1
    while overlaps0 != overlaps1:
        pairs_graph = PairCompositionGraph(pairs)
        PairPath = EulerianPath(pairs_graph)
        s0 = PairPath[0][0]
        s1 = PairPath[0][1]
        for p in PairPath[1:]:
            s0 += p[0][-1]
            s1 += p[1][-1]
        overlaps0 = s0[k + d :]
        s = len(s1) - (d + k)
        overlaps1 = s1[:s]
    nonoverlap = s1[s:]
    return s0 + nonoverlap


def PairCompositionGraph(pairs):
    pair_graph = dict()
    for p in pairs:
        key = (p[0][:-1], p[1][:-1])
        value = (p[0][1:], p[1][1:])
        if key not in pair_graph:
            pair_graph[key] = [value]
        else:
            pair_graph[key].append(value)
    return pair_graph
