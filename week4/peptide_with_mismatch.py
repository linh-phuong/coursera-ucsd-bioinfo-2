from collections import Counter
from week3.cyclopeptide import (
    expand_peptides,
    find_all_subpeptides,
    find_possible_spectrum,
    _parse_peptide_mass,
)
from tqdm import tqdm


def score_peptide(peptide, spectrum, linear=False):
    """score a peptide a gainst a given spectrum

    Args:
        peptide (str): a hypothetical peptide
        spectrum (list): a spectrum from an experiment

    Return:
        int: score of the peptide
    """
    theoretical_sp = find_possible_spectrum(peptide, linear=linear)
    return score_theoretical_spectrum(theoretical_sp, spectrum)


def score_theoretical_spectrum(th_spectrum, spectrum):
    """score a theoretical spectrum against a given parent spectrum
    Args:
        th_spectrum (list): the theoretical spectrum
        spectrum (list): the spectrum from the experiment
    Return:
        [int]: score
    """
    cnt_sp = Counter(spectrum)
    cnt_thsp = Counter(th_spectrum)
    c = 0
    for i in cnt_thsp:
        c += min([cnt_thsp[i], cnt_sp.get(i, 0)])
    return c


def count_common_multiset(a, b):
    pass


def find_leading_peptide(spectrum, N):
    """find the peptide with the highest score against an experiment spectrum using most common aa

    Args:
        spectrum (list): spectrum generated from experiments
        N (int): threshold to trim best candidates

    Returns:
        list: the highest score peptide (in the form of single masses)
    """
    massdict = _parse_peptide_mass()
    sp = sorted(set(massdict.values()))
    ls = _find_leaderboard_with_masseslist(spectrum, N, sp)
    return ls[0] if ls else None


def find_leaderboard(spectrum, N):
    """find the leaderboard with common amino acid

    Args:
        spectrum (list): spectrum generated from experiments
        N (int): threshold to trim best candidates

    Returns:
        list: a leaderboard of peptides with highest score (peptides are in the form of single masses)
    """
    massdict = _parse_peptide_mass()
    sp = sorted(set(massdict.values()))
    return _find_leaderboard_with_masseslist(spectrum, N, sp)


def _find_leaderboard_with_masseslist(spectrum, N, continuation):
    """find a leaderboard of peptide with highest score against an experiment spectrum

    Args:
        spectrum (list): spectrum generated from experiments
        N (int): threshold to trim best candidates
        continuation (list): possible amino acids that form the peptides

    Returns:
        list: a leader board of peptides (in single masses form)
    """
    leader_peptide = []
    parent_mass = sorted(spectrum)[-1]
    best_score = 0
    candidates_lb = [()]
    while candidates_lb:
        potential_candidates = expand_peptides(candidates_lb, continuation)
        desc = f"Iter {len(candidates_lb[0]) + 1:3}"
        total = len(candidates_lb) * len(continuation)
        candidates_lb = []
        for c in tqdm(potential_candidates, total=total, desc=desc):
            mass = sum(c)
            if mass > parent_mass:
                continue
            elif mass == parent_mass:
                th_spectrum = generate_spectrum_from_masses(c, linear=False)
                score = score_theoretical_spectrum(th_spectrum, spectrum)
                if score > best_score:
                    leader_peptide = [c]
                    best_score = score
                elif score == best_score:
                    leader_peptide.append(c)
                candidates_lb.append(c)
            else:
                candidates_lb.append(c)
        candidates_lb = select_best_candidates(candidates_lb, spectrum, N)
        # candidates_lb = (
        #     select_best_candidates(candidates_lb, spectrum, N) if candidates_lb else candidates_lb
        # )
    return leader_peptide


def generate_spectrum_from_masses(masses, linear=False):
    """generate the theoretical spectrum of a peptide (in the form of single masses)

    Args:
        masses (list): a peptide (in the form of singles masses)
        linear (bool, optional): False (Default): spectrum of a cyclopeptide
                                 True: spectrum of a linear peptide

    Returns:
        list: theoretical spectrum
    """
    subpeptides = (
        find_all_subpeptides(masses) if not linear else find_all_subpeptides(masses, linear=True)
    )
    subpeptides.append(masses)
    spectrum_masses = [0]
    for p in subpeptides:
        spectrum_masses.append(sum(p))
    return sorted(spectrum_masses)


def select_best_candidates(candidates, spectrum, N):
    """select the best candidate peptides

    Args:
        candidates (list): potential peptides to be selected
        spectrum (list): spectrum obtained from experiment
        N (int): select the Nth highest score peptides (including ties)

    Returns:
        list: best candidates
    """
    if not candidates:
        return []
    candidates_with_score = []
    sc_list = []
    for c in candidates:
        sp = generate_spectrum_from_masses(c, linear=True)
        sc = score_theoretical_spectrum(sp, spectrum)
        candidates_with_score.append((c, sc))
        sc_list.append(sc)
    sc_list = sorted(sc_list, reverse=True)
    thd = sc_list[N - 1] if len(sc_list) >= N else sc_list[-1]
    return [c for c, s in candidates_with_score if s >= thd]
