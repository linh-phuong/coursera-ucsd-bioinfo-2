from week1.debruijn import DeBruijnGraphFromReads, PathToGenome, StringComposition
from itertools import product


def EulerianCycle(graph):
    graph = {k: list(v) for k, v in graph.items()}
    subcycles = RandomCycle(graph)
    cycle = ConnectCycles(subcycles)
    return cycle


def RandomCycle(graph):
    edges_total = [iv for v in graph.values() for iv in v]
    subcycles = []

    while edges_total:
        newstart_node = edges_total[0]
        cycle = [newstart_node]
        next_node = None
        start_node = newstart_node
        while next_node != newstart_node:
            next_node = graph[start_node][0]
            cycle.append(next_node)
            edges_total.remove(start_node)
            graph[start_node].remove(next_node)
            start_node = next_node
        subcycles.append(cycle)
    return subcycles


def ConnectCycles(cycles):
    fc = cycles[0]
    rc = cycles[1:]
    while rc:
        for c in rc:
            e = list(set(fc).intersection(c))
            if e:
                i = fc.index(e[0])
                oc = PathOrdered(e[0], c)
                fc = fc[:i] + oc + fc[i + 1 :]  # noqa: E203
                rc.remove(c)
    return fc


def PathOrdered(first_element, path):
    assert type(path) is list, "path has to be a list"
    assert first_element in path, "first element not in the path"
    i = path.index(first_element)
    if path[i] == path[0]:
        return path
    new_path = path[i:] + path[1:i] + [path[i]]
    assert len(new_path) == len(path)
    return new_path


def _is_empty(G):
    return all(not v for v in G.values())


def IsEulerianCycle(path, graph):
    graph = {k: list(v) for k, v in graph.items()}
    for node in path:
        if graph[node]:
            graph[node].pop()
    return _is_empty(graph)


def UnbalancedNodes(graph):
    """Find imbalances nodes
    Return (None, None) if non imbalance nodes found
    Return (from, to) if found exactly two imbalance nodes
    Raise exception if found more than 2 imbalance
    """
    values_flat = sum(graph.values(), [])
    all_nodes = set(values_flat + list(graph.keys()))
    fr, to = None, None
    for n in all_nodes:
        out_nodes = len(graph[n]) if n in graph else 0
        in_nodes = values_flat.count(n)
        if out_nodes > in_nodes:
            if to is not None:
                raise Exception("Found more than one 'to' node")
            to = n
        elif out_nodes < in_nodes:
            if fr is not None:
                raise Exception("Found more than one 'fr' node")
            fr = n
    return fr, to


def EulerianPath(graph):
    graph = {k: list(v) for k, v in graph.items()}

    # Discover unbalanced nodes
    fr, to = UnbalancedNodes(graph)

    # Modify only if the graph is unbalanced
    if fr is not None and to is not None:
        # Connect unbalanced node
        if fr in graph:
            graph[fr].append(to)
        else:
            graph.update({fr: [to]})
    # Eulerian cycle
    cycle = EulerianCycle(graph)

    # Order cycle where the node with more outs in the beginning
    # and the node with more in locate in the end
    # only when the graph is unbalanced
    if fr is None or to is None:
        return cycle[:-1]
    for i in range(len(cycle) - 1):
        if cycle[i] == fr and cycle[i + 1] == to:
            si = i + 1
            return cycle[si:] + cycle[1:si]
    raise Exception("there are no eularian path")


def IsEulerianPathNb(path, graph):
    graph = {k: list(v) for k, v in graph.items()}
    for node in path:
        if graph[node]:
            next_node = graph[node].pop()
            if next_node not in graph.keys():
                graph.update({next_node: []})
    return _is_empty(graph)


def IsEulerianPath(path, graph):
    graph = {k: list(v) for k, v in graph.items()}
    k = len(list(graph.keys())[0])
    reads = StringComposition(k, path)
    for r in reads:
        if r not in graph:
            graph.update({r: []})
        if graph[r]:
            next_node = graph[r].pop()
            if next_node not in graph:
                graph.update({next_node: []})
    return _is_empty(graph), graph


def StringReconstruction(reads):
    dB = DeBruijnGraphFromReads(reads)
    path = EulerianPath(dB)
    Text = PathToGenome(path)
    return Text


def UniversalCircular(k):
    if k < 2:
        raise Exception("k has to be greater than 2 for universal circular from binary")
    strings = ["".join(i) for i in product("01", repeat=k)]
    G = DeBruijnGraphFromReads(strings)
    Ep = EulerianPath(G)
    univer_c = Ep[0]
    for s in Ep:
        univer_c += s[-1]
    return univer_c[: -k + 1]
