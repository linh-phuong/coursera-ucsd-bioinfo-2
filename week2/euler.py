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
