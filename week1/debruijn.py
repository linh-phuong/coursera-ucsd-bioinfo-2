def StringComposition(k, text):
    assert type(text) is str, "input is not a string"
    assert " " not in text, "there is space in the string"
    scanlen = len(text) - k + 1
    return sorted([text[i : i + k] for i in range(scanlen)])  # noqa: E203


def PathToGenome(path):
    assert type(path) is list, "input is not a list"
    text = path[0]
    for pattern in path[1:]:
        text += pattern[-1]
    return text


def OverlapGraph(patterns):
    assert type(patterns) is list, "input is not a list"
    prefix = []
    suffix = []
    for p in patterns:
        prefix.append(p[:-1])
        suffix.append(p[1:])
    connect_dict = dict()
    for i, suf_i in enumerate(suffix):
        for j, pref_j in enumerate(prefix):
            if suf_i == pref_j:
                text = patterns[i]
                if text in connect_dict.keys():
                    if patterns[j] not in connect_dict[text]:
                        connect_dict[text].append(patterns[j])
                else:
                    connect_dict[text] = [patterns[j]]
    return connect_dict


def DeBruijnGraph(k, text):
    assert type(text) is str, "input is not a string"
    patterns = StringComposition(k, text)
    prefix = []
    suffix = []
    for p in patterns:
        prefix.append(p[:-1])
        suffix.append(p[1:])
    path_dict = dict()
    for i, read in enumerate(patterns):
        node_fr = read[:-1]
        node_to = read[1:]
        if node_fr not in path_dict.keys():
            path_dict[node_fr] = [node_to]
        else:
            path_dict[node_fr].append(node_to)
    return path_dict


def DeBruijnGraphFromReads(patterns):
    assert type(patterns) is list, "input is not a list"
    prefix = []
    suffix = []
    for p in patterns:
        prefix.append(p[:-1])
        suffix.append(p[1:])
    path_dict = dict()
    for i, contig in enumerate(patterns):
        node = contig[:-1]
        edge = contig[1:]
        if node not in path_dict.keys():
            path_dict[node] = [edge]
        else:
            path_dict[node].append(edge)
    return path_dict
