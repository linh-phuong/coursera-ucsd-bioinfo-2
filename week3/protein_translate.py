def _parse_pr_dict():
    RNA_to_Pr = {}
    with open("week3/RNA_codon_table_1.txt") as fd:
        lines = fd.readlines()
        for l in lines:
            d = l.split(" ")
            RNA_to_Pr[d[0].strip()] = d[1].strip()
    return RNA_to_Pr


def translate_rna(rna):
    assert type(rna) is str, "input is not a string"
    assert " " not in rna, "there are spaces in the RNA"
    pr_string = ""
    codon_len = 3
    translate_dict = _parse_pr_dict()
    for i in range(0, len(rna), codon_len):
        codon = rna[i : i + codon_len]
        pr = translate_dict[codon] if len(codon) == 3 else ""
        pr_string += pr
    return pr_string


def find_reverse_complement(dna):
    assert type(dna) is str, "dna has to be string type"
    assert " " not in dna, "there is space in the dna"
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join(complement[n] for n in dna)


def find_all_peptides_encoding(noncoding_str, peptide):
    assert isinstance(noncoding_str, str)
    assert " " not in noncoding_str, "the noncoding strain should not have space"
    assert isinstance(peptide, str)
    assert " " not in peptide, "the peptide should not have space"
    assert "U" not in peptide, "the noncoding strain has to be a DNA"

    encoding_strains = []

    rna = noncoding_str.replace("T", "U")
    for i, fr in enumerate([rna, rna[1:], rna[2:]]):
        es = find_peptide_encoding(fr, peptide)
        for e in es:
            ori_e = e.replace("U", "T")

            # We only accept duplicated in the first frame
            if i == 0 or ori_e not in encoding_strains:
                assert ori_e in noncoding_str
                encoding_strains.append(ori_e)

    rv = find_reverse_complement(noncoding_str)[::-1].replace("T", "U")
    for fr in [rv, rv[1:], rv[2:]]:
        es = find_peptide_encoding(fr, peptide)
        for e in es:
            # transform to non coding strain
            ori_e = find_reverse_complement(e.replace("U", "T")[::-1])
            if ori_e not in encoding_strains:
                assert ori_e in noncoding_str, "the substrain is not in the noncoding strain"
                encoding_strains.append(ori_e)
    return encoding_strains


def find_peptide_encoding(transcripted_str, peptide):
    assert type(transcripted_str) is str
    assert " " not in transcripted_str, "the transcripted string should not have space"
    assert type(peptide) is str, "peptide has to be a string"
    assert " " not in peptide, "peptide should not have space"
    assert "T" not in transcripted_str, "the string has to be an RNA"

    codon_len = 3
    substrings = []
    pep_len = len(peptide)
    peptide_str, pr_indice = translate_rna_with_position(transcripted_str)
    pep_indice = find_substring_positions(peptide_str, peptide)
    if pep_indice:
        on_str_id = [pr_indice[i] for i in pep_indice]
    else:
        return []
    for i in on_str_id:
        ss = transcripted_str[i : i + codon_len * pep_len]
        substrings.append(ss)
    return substrings


def find_substring_positions(string, substring):
    return [i for i in range(len(string)) if string.startswith(substring, i)]


def translate_rna_with_position(rna):
    assert type(rna) is str, "input is not a string"
    assert " " not in rna, "there are spaces in the RNA"
    pr_string = ""
    codon_len = 3
    translate_dict = _parse_pr_dict()
    position = []
    for i in range(0, len(rna), codon_len):
        codon = rna[i : i + codon_len]
        pr = translate_dict[codon] if len(codon) == 3 else ""
        pr_string += pr if pr else "*"
        position.append(i)
    return pr_string, position
