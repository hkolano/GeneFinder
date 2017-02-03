# -*- coding: utf-8 -*-
"""

@author: Hannah Kolano
hannah.kolano@students.olin.edu

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):         # This one works
    """ Returns the complementary nucleotide
        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    nuc = list(nucleotide)
    count = 0
    complement = ''
    for element in nuc:
        if element == 'A':
            nuc[count] = 'T'
        elif element == 'T':
            nuc[count] = 'A'
        elif element == 'C':
            nuc[count] = 'G'
        elif element == 'G':
            nuc[count] = 'C'
        complement = complement + nuc[count]
        count = count + 1
    return complement


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    dna2 = get_complement(dna)
    dna3 = dna2[::-1]
    return str(dna3)


def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    num_codons = int(len(dna)/3)
    num = 0
    list_codons = []
    cut_dna = ''
    stop_index = -1
    while num < num_codons:
        num_start = int(num*3)
        num_end = int(num*3 + 3)
        list_codons.append(dna[num_start:num_end])
        num = num + 1
    for element in list_codons:
        if element == 'TAA' or element == 'TAG' or element == 'TGA':
            stop_index = list_codons.index(element)
            break

    if stop_index != -1:
        this_ORF = list_codons[0:stop_index]
        for element in this_ORF:
            cut_dna = cut_dna + element
        return cut_dna
    else:
        return dna


def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    condition = 1
    dna2 = dna
    list_of_ORFs = []
    while condition == 1:
        num_codons = int(len(dna2)/3)
        num = 0
        list_codons = []
        start_index = -1
        while num < num_codons:
            num_start = int(num*3)
            num_end = int(num*3 + 3)
            list_codons.append(dna2[num_start:num_end])
            num = num + 1
        for element in list_codons:
            if element == 'ATG':
                start_index = list_codons.index(element)
                break
        if start_index != -1:
            this_orf = rest_of_ORF(dna2[start_index*3:])
            list_of_ORFs.append(this_orf)
            dna2 = dna2[start_index*3+len(this_orf)+3:]
        else:
            condition = 2

    return list_of_ORFs


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    dnatwo = dna[1:]
    dnathree = dna[2:]
    all_ORFs = find_all_ORFs_oneframe(dna)
    all_ORFs.extend(find_all_ORFs_oneframe(dnatwo))
    all_ORFs.extend(find_all_ORFs_oneframe(dnathree))
    if [] in all_ORFs:
        all_ORFs.remove([])
    return all_ORFs


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    reverse_dna = get_reverse_complement(dna)
    orfs = find_all_ORFs(dna)
    backward_orfs = find_all_ORFs(reverse_dna)
    orfs.extend(backward_orfs)
    return orfs


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    MyList = find_all_ORFs_both_strands(dna)
    if MyList == []:
        maximum_orf = 'none'
    else:
        maximum_orf = max(MyList, key=len)
    return maximum_orf


def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    count = 0
    orfs = []
    while count < num_trials:
        current_strand = shuffle_string(dna)
        longest_in_shuffle = longest_ORF(current_strand)
        orfs.append(longest_in_shuffle)
        count = count + 1
    longest = max(orfs, key=len)
    return len(longest)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    num_codons = int(len(dna)/3)
    num = 0
    list_codons = []
    aacids = ''
    while num < num_codons:
        num_start = int(num*3)
        num_end = int(num*3 + 3)
        list_codons.append(dna[num_start:num_end])
        num = num + 1
    for element in list_codons:
        thing = aa_table[element]
        aacids = aacids + thing
    return aacids


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    all_orfs_both_strands = find_all_ORFs_both_strands(dna)
    longest_fake_orf = longest_ORF_noncoding(dna, 20)
    for element in all_orfs_both_strands:
        if len(element) > longest_fake_orf:
            a_a_string = coding_strand_to_AA(element)
        else:
            a_a_string = 'not longer than shuffle'
        print(a_a_string)


if __name__ == "__main__":
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna))
    # import doctest
    # sdoctest.testmod()
    # doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)
