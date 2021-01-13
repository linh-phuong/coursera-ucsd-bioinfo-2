from week2.euler import DeBruijnGraphFromReads, _is_empty
from collections import Counter


def ContigGeneration(reads):
    dB = DeBruijnGraphFromReads(reads)
    in_out = in_out_count(dB)
    start_nodes = []
    for k in in_out.keys():
        if in_out[k][0] != in_out[k][1] or in_out[k] > [1, 1]:
            start_nodes += [k] * in_out[k][0]
    contigs = []
    kmer = len(start_nodes[0])
    while sorted(start_nodes):
        n_start = start_nodes.pop()
        n_to = dB[n_start].pop()
        n_start += n_to[-1]
        while in_out[n_to] == [1, 1] and dB[n_to]:
            n_to = dB[n_to].pop()
            n_start += n_to[-1]
        contigs.append(n_start)
    if not _is_empty(dB):
        loop_nodes = sorted([k for k in dB.keys() if dB[k]])
    else:
        return contigs
    # while loop_nodes:
    #     n_start = loop_nodes.pop()
    #     n_to = dB[n_start].pop()
    #     n_start += n_to[-1]
    #     loop_nodes.remove(n_to)
    #     n_to = dB[n_to].pop()
    #     n_start += n_to[-1]
    #     contigs.append(n_start)
    while loop_nodes:
        n_start = loop_nodes.pop()
        n_to = dB[n_start].pop()
        while n_to != n_start[0:kmer]:
            n_start += n_to[-1]
            loop_nodes.remove(n_to)
            n_to = dB[n_to].pop()
        n_start += n_to[-1]
        contigs.append(n_start)

    assert _is_empty(dB)
    return contigs


def MaxBranch(dB):
    in_out = in_out_count(dB)
    start_nodes = []
    for k in in_out.keys():
        if in_out[k][0] != in_out[k][1] or in_out[k] > [1, 1]:
            start_nodes += [k] * in_out[k][0]
    contigs = []
    while sorted(start_nodes):
        n_start = start_nodes.pop()
        n_to = dB[n_start].pop()
        n_start += "->"
        n_start += n_to
        while in_out[n_to] == [1, 1] and dB[n_to]:
            n_to = dB[n_to].pop()
            n_start += "->"
            n_start += n_to
        contigs.append(n_start)
    if not _is_empty(dB):
        loop_nodes = sorted([k for k in dB.keys() if dB[k]])
    else:
        return contigs

    while loop_nodes:
        new_start = loop_nodes.pop()
        n_to = dB[new_start].pop()
        l_contig = new_start
        while n_to != new_start:
            l_contig += "->"
            l_contig += n_to
            loop_nodes.remove(n_to)
            n_to = dB[n_to].pop()
        l_contig += "->"
        l_contig += n_to
        contigs.append(l_contig)
    assert _is_empty(dB)
    return contigs


def in_out_count(debruijn):
    nodes_to = [iv for v in debruijn.values() for iv in v]
    incount = Counter(nodes_to)
    allnodes = set(list(debruijn.keys()) + nodes_to)
    d = dict()
    for n in allnodes:
        if n in debruijn and n in incount:
            d[n] = [len(debruijn[n]), incount[n]]
        elif n in debruijn and n not in incount:
            d[n] = [len(debruijn[n]), 0]
        elif n not in debruijn and n in incount:
            d[n] = [0, incount[n]]
        else:
            raise Exception("the node does not exist")
    return d


def get_max_out_count(count):
    # return sorted(count.items(), key=lambda i: i[1][0])[0]
    c = 0
    maxkey = None
    for k, (outcount, _) in sorted(count.items(), key=lambda i: i[0]):
        if outcount > c:
            maxkey = k
            c = outcount
    return maxkey
