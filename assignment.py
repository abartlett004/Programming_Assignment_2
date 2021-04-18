from itertools import product, combinations
from collections import Counter

translation_table = {'TTT': 'F', 'TCT': 'S', 'TAT': 'Y', 'TGT': 'C',
                     'TTC': 'F', 'TCC': 'S', 'TAC': 'Y', 'TGC': 'C',
                     'TTA': 'L', 'TCA': 'S', 'TAA': '*', 'TGA': '*',
                     'TTG': 'L', 'TCG': 'S', 'TAG': '*', 'TGG': 'W',
                     'CTT': 'L', 'CCT': 'P', 'CAT': 'H', 'CGT': 'R',
                     'CTC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',
                     'CTA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',
                     'CTG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',
                     'ATT': 'I', 'ACT': 'T', 'AAT': 'N', 'AGT': 'S',
                     'ATC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',
                     'ATA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',
                     'ATG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',
                     'GTT': 'V', 'GCT': 'A', 'GAT': 'D', 'GGT': 'G',
                     'GTC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',
                     'GTA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',
                     'GTG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'
                     }

# nomenclature for degenerate codons
expanded_code = {'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'],
                 'W': ['A', 'T'], 'S': ['C', 'G'], 'M': ['A', 'C'], 'K': ['G', 'T'], 'R': ['A', 'G'], 'Y': ['C', 'T'],
                 'B': ['C', 'G', 'T'], 'D': ['A', 'G', 'T'], 'H': ['A', 'C', 'T'], 'V': ['A', 'C', 'G'],
                 'N': ['A', 'C', 'G', 'T']
                 }

# helpful for validating input
valid_nucleotides = 'ACGTWSMKRYBDHVN'
valid_aa = 'GAVLIMFWPSTCYNQDEKRH*'

# calculate translation table including codons with nontraditional bases
# value of each codon will be a counter where the keys are the amino acids it can
# generate and the value will be the number of ways it can generate that AA
expanded_tt = {}

# examine each possible codon using expanded code
for codon in product(expanded_code.keys(), repeat = 3):
    base1 = expanded_code[codon[0]]
    base2 = expanded_code[codon[1]]
    base3 = expanded_code[codon[2]]
    
    # if current codon is degenerate, find all traditional codons it can encode for
    trad_codons = []
    for base1_opt in base1:
        for base2_opt in base2:
            for base3_opt in base3:
                trad_codons.append(base1_opt + base2_opt + base3_opt)

    # use possible traditional codons to determine what amino acids current codon can code for
    poss_aas = [translation_table[codon] for codon in trad_codons]
    expanded_tt[''.join(codon)] = Counter(poss_aas)

def get_codon_for_amino_acids(amino_acids):
    """
    :param amino_acids: set
        the amino acids we want to code for, i.e. {'A','I','V'}
    :rtype: set, float
        returns two values the set of most efficient codons for the input set list, e.g. {'RYA', 'RYH', 'RYC', 'RYW', 'RYM', 'RYY', 'RYT'} and the achieved efficiency e.g. 0.75
    """
    # validate input
    for aa in amino_acids:
        if aa not in valid_aa:
            raise('Input contains an invalid amino acid')
    
    # initalize tracker for maximum efficiency and codons that achieve it
    max_eff = 0
    max_eff_codons = set()
    # check each possible degenerate codon
    for codon, aa in expanded_tt.items():

        # if current degenerate codon can generate all input AAs and cannot generate a stop
        if all(aa in expanded_tt[codon] for aa in amino_acids) and '*' not in expanded_tt[codon]:
            
            # determine efficiency of using current codon to generate input AAs
            counter = expanded_tt[codon]
            occurrences = sum(counter[aa] for aa in set(amino_acids))
            eff = occurrences/sum(counter[aa] for aa in counter.keys())

            # if efficiency of current codon is greater than previously observed max efficiency,
            # reassign max efficiency, clear list of max efficient codons, and add current codon to list
            if eff > max_eff:
                max_eff_codons = set([codon])
                max_eff = eff
            # if efficiency of current codon matches previously observed max efficiency, add it to list
            elif eff == max_eff:
                max_eff_codons.add(codon)

    return(max_eff_codons, max_eff)


def truncate_list_of_amino_acids(amino_acids):
    """
    :param amino_acids: set
        the amino acids we want to code for, i.e. {'A','I','V'}
    :rtype: set
        the set of sets of amino acids that can be coded with 100% efficiency, i.e. {frozenset({'V', 'A'}), frozenset({'V', 'I'})}
    """
    # validate input
    for aa in amino_acids:
        if aa not in valid_aa:
            raise('Input contains an invalid amino acid')

    # determine whether complete input set can be encoded by codons with 100% efficiency
    eff = get_codon_for_amino_acids(amino_acids)[1]
    eff_aas = set()
    if eff == 1:
        eff_aas = frozenset(amino_acids)
    # if not, use helper function to examine subsets of input set of AAs that exclude 1 AA
    else: 
        eff_aas = truncate_helper(amino_acids, 1)

    return(eff_aas)

            
def truncate_helper(amino_acids, num_to_exclude):
    # if all possible combinations of AAs have been tested, return empty set
    eff_aas = set()
    if (num_to_exclude == len(amino_acids)):
        return eff_aas

    # find subsets of input AAs where given number of AAs are excluded
    subsets = combinations(amino_acids, len(amino_acids) - num_to_exclude)

    # determine whether any of these subsets can be encoded with 100% efficiency
    for subset in subsets:
        eff = get_codon_for_amino_acids(subset)[1]
        # if so, return the subsets that can 
        if eff == 1:
            eff_aas.add(frozenset(subset))
    # if not, examine combinations that exclude an additional amino acid
    if len(eff_aas) == 0:
        eff_aas = truncate_helper(amino_acids, num_to_exclude + 1)

    return eff_aas

if __name__ == "__main__":
    # # using sets instead of lists throughout the code since the order doesn't matter and all items should be unique
    assert get_codon_for_amino_acids({'A', 'I', 'V'}) == ({'RYA', 'RYH', 'RYC', 'RYW', 'RYM', 'RYY', 'RYT'}, 0.75)
    assert get_codon_for_amino_acids({'M', 'F'}) == ({'WTS', 'WTK', "WTB"}, 0.5)
    
    # "frozenset" here since this seems to be the only way to get a set of sets - see https://stackoverflow.com/questions/5931291/how-can-i-create-a-set-of-sets-in-python
    assert truncate_list_of_amino_acids({'A', 'V', 'I'}) == {frozenset({'V', 'A'}), frozenset({'V', 'I'})}
