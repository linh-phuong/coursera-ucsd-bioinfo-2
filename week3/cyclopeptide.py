from collections import Counter
from itertools import product
from tqdm import tqdm
from math import factorial


def _parse_peptide_mass():
    massdict = {}
    with open("week3/integer_mass_table.txt") as fd:
        data = fd.readlines()
        for line in data:
            kv = line.strip().split(" ")
            massdict[kv[0].strip()] = int(kv[1].strip())
    return massdict


def find_subpeptides_with_len(peptides, length):
    """find all subpeptides of a cyclo peptide

    Args:
        peptides (str): a cyclo peptide (in linear format)
        length (int): length of the subpeptides

    Returns:
        list: all possible subpeptides
    """
    assert len(peptides) >= 2
    assert len(peptides) > length, "number of subpeptides greater than the length of the peptide"
    newp = peptides + peptides[:length]
    subpeptides = []
    scanlen = len(newp) - length
    for i in range(scanlen):
        subpeptides.append(newp[i : i + length])
    return subpeptides


def find_linear_subpeptides_with_len(peptide, length):
    """find all subpeptides with a certain length of a cyclo peptide that has been cut one time

    Args:
        peptide (sequence): the peptide can be given by a string or list of single masses  
        length (int): lenght of the subpeptides 

    Returns:
        list: all possible subpeptides with a certain length
    """
    subpeptides = []
    scanlen = len(peptide) - length + 1
    for i in range(scanlen):
        subpeptides.append(peptide[i : i + length])
    return subpeptides


def find_all_subpeptides(peptide, linear=False):
    """find all possible subpeptide of a cyclopeptide
    Args:
    peptide (sequence): a cyclopeptide that can be in the format of a string of single amino acid 
                            or a list of single masse of the amino acid
    linear (bool, optional): False (default): all possible subpeptides of a cyclopeptide
                             True: all possible subpeptides of a cyclopeptide that is cut once

    Returns:
        list: all possible subpeptides 
    """
    psizes = len(peptide)
    subpeptides = []
    for s in range(1, psizes):
        if not linear:
            subpeptides += find_subpeptides_with_len(peptide, s)
        else:
            subpeptides += find_linear_subpeptides_with_len(peptide, s)
    return subpeptides


def find_possible_spectrum(peptide, linear=False):
    """find the theoretical spectrum of a peptide

    Args:
        peptide (str): a cyclopeptide 
        linear (bool, optional): False(default): all possible subpeptides of a cyclopeptide
                                 True: all possible subpeptides of a cyclopeptide that is cut once


    Returns:
        list: theoretical spectrum of the peptide
    """
    massdict = _parse_peptide_mass()
    subpeptides = (
        find_all_subpeptides(peptide) if not linear else find_all_subpeptides(peptide, linear=True)
    )
    subpeptides.append(peptide)
    mass = [0]
    for p in subpeptides:
        mass.append(sum([massdict[i] for i in p]))
    return sorted(mass)


def count_peptide_given_mass_rc(mass, mass_list):
    """count the number of peptide strains with the given mass

    Args:
        mass ([int]): total mass of a peptide strain

    Returns:
        [int]: number of peptide strains with the mass

    """
    if mass == 0:
        return 1
    if mass < 0:
        return 0
    c = 0
    for m in mass_list:
        cm = mass - m
        # new_mass = [i for i in mass_list if i <= m]
        # c += count_peptide_given_mass(cm, new_mass)
        c += count_peptide_given_mass_rc(cm, mass_list)
    return c


def find_nb_combinations(a_list):
    cnt = Counter(a_list)
    counts = cnt.values()
    all_nb = factorial(sum(counts))
    rep_nb = 1
    for i in counts:
        rep_nb = rep_nb * factorial(i)
    return all_nb / rep_nb


def generate_peptide_given_mass(mass, mass_list):
    return _generate_peptide_given_mass_recursive(mass, sorted(mass_list, reverse=True), [])


def _generate_peptide_given_mass_recursive(mass, mass_list, selected):
    if mass == 0:
        return [selected]
    if mass < 0:
        return []
    solutions = []
    for m in mass_list:
        new_mass_list = [i for i in mass_list if i <= m]
        new_selected = [*selected, m]
        branch_solutions = _generate_peptide_given_mass_recursive(
            mass - m, new_mass_list, new_selected
        )
        solutions += branch_solutions
    return solutions


def count_peptides_given_mass_rc2(mass, mass_list):
    unique_peptides = generate_peptide_given_mass(mass, mass_list)
    c = 0
    for p in unique_peptides:
        c += find_nb_combinations(p)
    return c


def _is_a_in_b(a, b):
    cnt_a = Counter(a)
    cnt_b = Counter(b)
    return all(v <= cnt_b.get(k, 0) for k, v in cnt_a.items())


def expand_peptides(peptides, aa):
    c_pd = product(peptides, aa)
    for l0, l1 in c_pd:
        yield (*l0, l1)


def find_cyclopeptide_spectrum(peptide, aa):
    """ find the spectrum of a cyclopeptide spectrum when adding an amino acid

    Args:
        peptide (tuple): a cylopeptide that is cut into a string, each element is the mass of an amino acid
        aa (int): integer mass of additional aminoacid
    """
    isinstance(peptide, list)
    isinstance(aa, int)
    scanlen = len(peptide)
    return [sum(peptide[i:]) + aa for i in range(scanlen)]


def expand_peptides_within_spectrum(peptides, aa, cnt_spectrum):
    c_pd = product(peptides, aa)
    for l0, l1 in c_pd:
        peptide_spectrum = find_cyclopeptide_spectrum(l0, l1)
        if _is_a_in_b(peptide_spectrum, cnt_spectrum):
            yield (*l0, l1)


def find_cyclopeptide(spectrum):
    peptides = set()
    sp = _parse_peptide_mass()
    parent = spectrum[-1]
    spectrum_counter = Counter(spectrum)
    continuation = sorted([i for i in spectrum_counter if i in set(sp.values())])

    candidates = [()]

    while candidates:
        potential_candidates = expand_peptides_within_spectrum(
            candidates, continuation, spectrum_counter
        )
        desc = f"Iter {len(candidates[0]) + 1:3}"
        total = len(candidates) * len(continuation)
        candidates = []
        for c in tqdm(potential_candidates, total=total, desc=desc):
            mass = sum(c)
            if mass > parent or not _is_a_in_b(c, spectrum_counter):
                continue
            if mass < parent and _is_a_in_b([mass], spectrum_counter):
                candidates.append(c)
            elif mass == parent:
                peptides.add(c)
    return peptides


# def find_combinations_given_mass(mass: int, elements: Sequence[int]) -> Sequence[Tuple[int]]:
#     pass
