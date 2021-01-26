from week2.euler import EulerianPath


def StringReconstructFromPairPath(pairs, k, d):
    # This function still needs revision,
    # EulerianPath can produce many results when the debruijn graph is highly tangled
    # However only cycles that show overlaps between the two pairs can be accepted
    overlaps0 = 0
    overlaps1 = 1
    c = 0
    while overlaps0 != overlaps1:
        c += 1
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
        if c == 100:
            return None
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
