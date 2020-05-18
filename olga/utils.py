#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Module for utility functions.

    Copyright (C) 2018 Zachary Sethna

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
    
For documentation on how to define custom amino acid alphabets, look at the
docstring of construct_codons_dict.

@author: zacharysethna
"""
from __future__ import print_function, division
import numpy as np

    
#%% Translation functions
def nt2aa(ntseq):
    """Translate a nucleotide sequence into an amino acid sequence.

    Parameters
    ----------
    ntseq : str
        Nucleotide sequence composed of A, C, G, or T (uppercase or lowercase)

    Returns
    -------
    aaseq : str
        Amino acid sequence
    
    Example
    --------
    >>> nt2aa('TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC')
    'CAWSVAPDRGGYTF'
        
    """
    nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}
    aa_dict ='KQE*TPASRRG*ILVLNHDYTPASSRGCILVFKQE*TPASRRGWMLVLNHDYTPASSRGCILVF'
    
    return ''.join([aa_dict[nt2num[ntseq[i]] + 4*nt2num[ntseq[i+1]] + 16*nt2num[ntseq[i+2]]] for i in range(0, len(ntseq), 3) if i+2 < len(ntseq)])

def nt2codon_rep(ntseq):
    """Represent nucleotide sequence by sequence of codon symbols.

    'Translates' the nucleotide sequence into a symbolic representation of
    'amino acids' where each codon gets its own unique character symbol. These 
    characters should be reserved only for representing the 64 individual 
    codons --- note that this means it is important that this function matches 
    the corresponding function in the preprocess script and that any custom 
    alphabet does not use these symbols. Defining symbols for each individual
    codon allows for Pgen computation of inframe nucleotide sequences. 

    Parameters
    ----------

    ntseq : str
        A Nucleotide sequence (normally a CDR3 nucleotide sequence) to be 
        'translated' into the codon - symbol representation. Can be either
        uppercase or lowercase, but only composed of A, C, G, or T.

    Returns
    -------
    codon_rep : str
        The codon - symbolic representation of ntseq. Note that if 
        len(ntseq) == 3L --> len(codon_rep) == L
    
    Example
    --------
    >>> nt2codon_rep('TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC')
    '\xbb\x96\xab\xb8\x8e\xb6\xa5\x92\xa8\xba\x9a\x93\x94\x9f'
        
    """
    
    #
    nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'a': 0, 'c': 1, 'g': 2, 't': 3}
    #Use single characters not in use to represent each individual codon --- this function is called in constructing the codon dictionary
    codon_rep ='\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf'
    
    return ''.join([codon_rep[nt2num[ntseq[i]] + 4*nt2num[ntseq[i+1]] + 16*nt2num[ntseq[i+2]]] for i in range(0, len(ntseq), 3) if i+2 < len(ntseq)])

#%% Functions that delete/insert reverse palindromes nucleotide sequences
def cutR_seq(seq, cutR, max_palindrome):
    """Cut genomic sequence from the right.

    Parameters
    ----------
    seq : str
        Nucleotide sequence to be cut from the right
    cutR : int
        cutR - max_palindrome = how many nucleotides to cut from the right.
        Negative cutR implies complementary palindromic insertions.
    max_palindrome : int
        Length of the maximum palindromic insertion.

    Returns
    -------
    seq : str
        Nucleotide sequence after being cut from the right
    
    Examples
    --------
    >>> cutR_seq('TGCGCCAGCAGTGAGTC', 0, 4)
    'TGCGCCAGCAGTGAGTCGACT'
    >>> cutR_seq('TGCGCCAGCAGTGAGTC', 8, 4)
    'TGCGCCAGCAGTG'
    
    """
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} #can include lower case if wanted
    if cutR < max_palindrome:
        seq = seq + ''.join([complement_dict[nt] for nt in seq[cutR - max_palindrome:]][::-1]) #reverse complement palindrome insertions
    else:
        seq = seq[:len(seq) - cutR + max_palindrome] #deletions
    
    return seq
    
def cutL_seq(seq, cutL, max_palindrome):
    """Cut genomic sequence from the left.

    Parameters
    ----------
    seq : str
        Nucleotide sequence to be cut from the right
    cutL : int
        cutL - max_palindrome = how many nucleotides to cut from the left.
        Negative cutL implies complementary palindromic insertions.
    max_palindrome : int
        Length of the maximum palindromic insertion.

    Returns
    -------
    seq : str
        Nucleotide sequence after being cut from the left
    
    Examples
    --------
    >>> cutL_seq('TGAACACTGAAGCTTTCTTT', 8, 4)
    'CACTGAAGCTTTCTTT'
    >>> cutL_seq('TGAACACTGAAGCTTTCTTT', 0, 4)
    'TTCATGAACACTGAAGCTTTCTTT'
    
    """
    
    complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} #can include lower case if wanted
    if cutL < max_palindrome:
        seq = ''.join([complement_dict[nt] for nt in seq[:max_palindrome - cutL]][::-1]) + seq #reverse complement palindrome insertions
    else:
        seq = seq[cutL-max_palindrome:] #deletions
    
    return seq

#%% Define codon dicts that define the 'amino acid' alphabet
def construct_codons_dict(alphabet_file = None):
    """Generate the sub_codons_right dictionary of codon suffixes.

    syntax of custom alphabet_files:
        
    char: list,of,amino,acids,or,codons,separated,by,commas
    
    
    Parameters
    ----------
    alphabet_file : str
        File name for a custom alphabet definition. If no file is provided, the 
        default alphabet is used, i.e. standard amino acids, undetermined amino 
        acids (B, J, X, and Z), and single codon symbols.

    Returns
    -------
    codons_dict : dict
        Dictionary, keyed by the allowed 'amino acid' symbols with the values 
        being lists of codons corresponding to the symbol.
    
    """
    
    #Some symbols can't be used in the CDR3 sequences in order to allow for
    #regular expression parsing and general manipulation.
    protected_symbols = [' ', '\t', '\n', '\x0b', '\x0c', '\r', ':', ',', ';', '[', ']', '{', '}', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9']
    
    
    #construct list of all 64 codons
    codons = [i + j + k for i in 'ACGT' for j in 'ACGT' for k in 'ACGT']
    
    codons_dict = {}
    #add standard amino acids symbols to the dict (i.e. 'ACDEFGHIKLMNPQRSTVWY*').
    #these symbols CANNOT be overwritten by custom alphabet files
    for codon in codons:
        codons_dict[nt2aa(codon)] = codons_dict.get(nt2aa(codon), []) + [codon]
                 
    #add single codon symbols to allow for inframe ntseq pgen computation
    #'\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f\xa0\xa1\xa2\xa3\xa4\xa5\xa6\xa7\xa8\xa9\xaa\xab\xac\xad\xae\xaf\xb0\xb1\xb2\xb3\xb4\xb5\xb6\xb7\xb8\xb9\xba\xbb\xbc\xbd\xbe\xbf'
    #these symbols CANNOT be overwritten by custom alphabet files
    for codon in codons:
        codons_dict[nt2codon_rep(codon)] = [codon]  
    
    #Check to see if custom alphabet file is supplied, else use default alphabet
    
    #Include standard ambigious amino acids.
    #these symbols CAN be overwritten by custom alphabet files
    expanded_alphabet = {}
    expanded_alphabet['B'] = ['D','N']
    expanded_alphabet['J'] = ['I', 'L']
    expanded_alphabet['X'] = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    expanded_alphabet['Z'] = ['E', 'Q']
    if alphabet_file is not None: #Use custom alphabet file definitions
        alphabet_f = open(alphabet_file, 'r')
        for line in alphabet_f:
            #assumed syntax is of a line is:
            #s: a1, a2, a3, a4, a5, ..., aN
            #where s is a single character symbol that isn't reserved, and all 
            #of the a's are either amino acid symbols or codons. Whitespaces
            #will be stripped as will brackets if the a's are presented as a
            #list.
            c_symbol = line.split(':', 1)[0].strip(''.join(protected_symbols))
            #Note there shouldn't be any additional colons -- this is a protected symbol.
            c_aa_codon_list_str = line.split(':', 1)[1]
            expanded_alphabet[c_symbol] = [x.strip(''.join(protected_symbols)) for x in c_aa_codon_list_str.split(',')]
        
        alphabet_f.close()
    
    
    for symbol in expanded_alphabet.keys():
        #Double check that the symbol isn't already used (important particularly for the single codon representation)
        if symbol in codons_dict.keys():
            print(symbol + " is already used as an 'amino acid' symbol for codons: ")
            print(codons_dict[symbol])
            continue
        elif not len(symbol) == 1: #Check that the custom symbol is a single character
            print("Can't use " + symbol + " as a custom 'amino acid' definitions as such symbols must be single characters.")
            continue
        elif symbol in protected_symbols: #This elif shouldn't trigger due to the stripping of protected symbols.
            print(symbol + " is a protected character")
        current_codon_collection = set()
        for x in expanded_alphabet[symbol]:
            if x in codons_dict.keys(): #Check if reference to an amino acid or other amino acid symbol
                current_codon_collection = current_codon_collection.union(codons_dict[x]) #If so, add those codons to the new collection
            elif x.upper() in codons: #Check if specifying a single codon
                current_codon_collection.add(x.upper()) #If so, add the codon to the new collection
            elif len(x) == 0: #fully stripped away
                continue
            else: #If not, don't recognize the addition and continue.
                print('Unfamiliar amino acid symbol or codon: ' + x)
                continue
        codons_dict[symbol] = list(current_codon_collection)
    
                
    return codons_dict

def generate_sub_codons_left(codons_dict):
    """Generate the sub_codons_left dictionary of codon prefixes.

    Parameters
    ----------
    codons_dict : dict
        Dictionary, keyed by the allowed 'amino acid' symbols with the values 
        being lists of codons corresponding to the symbol.

    Returns
    -------
    sub_codons_left : dict
        Dictionary of the 1 and 2 nucleotide prefixes (read from 5') for 
        each codon in an 'amino acid' grouping
        
    """
    sub_codons_left = {}
    for aa in codons_dict.keys():
        sub_codons_left[aa] = list(set([x[0] for x in codons_dict[aa]] + [x[:2] for x in codons_dict[aa]]))
        
    return sub_codons_left
    
def generate_sub_codons_right(codons_dict):
    """Generate the sub_codons_right dictionary of codon suffixes.

    Parameters
    ----------
    codons_dict : dict
        Dictionary, keyed by the allowed 'amino acid' symbols with the values 
        being lists of codons corresponding to the symbol.

    Returns
    -------
    sub_codons_right : dict
        Dictionary of the 1 and 2 nucleotide suffixes (read from 5') for 
        each codon in an 'amino acid' grouping.
        
    """
    sub_codons_right = {}
    for aa in codons_dict.keys():
        sub_codons_right[aa] = list(set([x[-1] for x in codons_dict[aa]] + [x[-2:] for x in codons_dict[aa]]))
        
    return sub_codons_right        

def determine_seq_type(seq, aa_alphabet):
    """Determine the type of a sequence.
    
    Parameters
    ----------

    seq : str
        Sequence to be typed.
    aa_alphabet : str
        String of all characters recoginized as 'amino acids'. (i.e. the keys
        of codons_dict: aa_alphabet = ''.join(codons_dict.keys())  )

    Returns
    -------
    seq_type : str
        The type of sequence (ntseq, aaseq, regex, None) seq is.
    
    Example
    --------
    >>> determine_seq_type('TGTGCCAGCAGTTCCGAAGGGGCGGGAGGGCCCTCCCTGAGAGGTCATGAGCAGTTCTTC', aa_alphabet)
    'ntseq'
    >>> determine_seq_type('CSARDX[TV]GNX{0,}', aa_alphabet)
    'regex
    
    """
    
    if all([x in 'ACGTacgt' for x in seq]):
        return 'ntseq'
    elif all([x in aa_alphabet for x in seq]):
        return 'aaseq'
    elif all([x in aa_alphabet + '[]{}0123456789,' for x in seq]):
        return 'regex'
    
#%%
#If using the steady-state distribution for first nucleotide probabilities we include a function to compute it
def calc_steady_state_dist(R):
    """Calculate the steady state dist of a 4 state markov transition matrix.

    Parameters
    ----------
    R : ndarray
        Markov transition matrix

    Returns
    -------
    p_ss : ndarray
        Steady state probability distribution
        
    """
    #Calc steady state distribution for a dinucleotide bias matrix
    
    w, v = np.linalg.eig(R)
    
    for i in range(4):
        if np.abs(w[i] - 1) < 1e-8:
            return np.real(v[:, i] / np.sum(v[:, i]))
    return -1

#%% Gene recognition
def gene_to_num_str(gene_name, gene_type):
    """Strips excess gene name info to number string.
    
    Parameters
    ----------
    gene_name : str
        Gene or allele name
    gene_type : char
        Genomic cassette type. (i.e. V, D, or J)

    Returns
    -------
    num_str : str
        Reduced gene or allele name with leading zeros and excess
        characters removed.
        
    """

    num_str = gene_type.lower().join([g.lstrip('0') for g in gene_name.lower().split(gene_type.lower())[1:]])
    num_str = '-'.join([g.lstrip('0') for g in num_str.split('-')])
    num_str = '*'.join([g.lstrip('0') for g in num_str.split('*')])
    return gene_type.lower() + num_str.replace('/', '') # get rid of '/' too

#%% Entropy Functions
def calc_S(P, base):
    """Calculate the entropy of a probability distribution in base units.
    
    Parameters
    ----------
    P : ndarray
        Probability distribution
    base : float
        Base used to set entropy units. Suggestions:
        base = 2 (entropy in bits)
        base = e (entropy in nats)
        base = 10 (entropy in dits)

    Returns
    -------
    S : float
        Entropy of P
        
    """
    return np.sum([-p*np.log(p) for p in P.flatten() if p > 0])/np.log(base)

def calc_S_single_gene(PG, PdelG_given_G, base = 2):
    """Calculate the entropy of an independent gene usage and deletion dists.
    
    Parameters
    ----------
    PG : ndarray
        Probability distribution of the usage of G-type genes (e.g PV)
    PdelG_given_G : ndarray
        Deletion profiles of each G gene conditioned on the gene 
        (e.g. PdelV_given_V or PdelDldelDr_given_D)
    base : float
        Base used to set entropy units. Suggestions:
        base = 2 (entropy in bits, default)
        base = e (entropy in nats)
        base = 10 (entropy in dits)

    Returns
    -------
    SG : float
        Entropy of PG (entropy of gene choice)
    SdelG : float
        Conditional entropy of PdelG_given_G (entropy of deletions)
        
    """
    
    SG = calc_S(PG, base)
    SdelG = 0
    for g, pg in enumerate(PG):
        if pg == 0:
            continue
        SdelG += pg * calc_S(PdelG_given_G[..., g], base)
        
    return SG, SdelG

def calc_S_joint_genes(PG1G2, PdelG1_given_G1, PdelG2_given_G2, base = 2):
    """Calculate the entropy of a joint gene usage and deletion dists.
    
    Parameters
    ----------
    PG1G2 : ndarray
        Joint probability distribution of the usage of G1-type genes and 
        G2-type genes (e.g. PDJ or PVJ)
    PdelG1_given_G1 : ndarray
        Deletion profiles of each G1 gene conditioned on the gene 
        (e.g. PdelV_given_V or PdelDldelDr_given_D)
    PdelG2_given_G2 : ndarray
        Deletion profiles of each G2 gene conditioned on the gene 
        (e.g. PdelV_given_V or PdelDldelDr_given_D)
    base : float
        Base used to set entropy units. Suggestions:
        base = 2 (entropy in bits, default)
        base = e (entropy in nats)
        base = 10 (entropy in dits)

    Returns
    -------
    SG1G2 : float
        Entropy of PG1G2 (entropy of gene choice)
    SdelG1 : float
        Conditional entropy of PdelG1_given_G1 (entropy of G1 deletions)
    SdelG2 : float
        Conditional entropy of PdelG2_given_G2 (entropy of G2 deletions)
        
    """

    PG1 = np.sum(PG1G2, axis = 1)
    PG2 = np.sum(PG1G2, axis = 0)
    
    SG1G2 = calc_S(PG1G2, base)
    SdelG1 = 0
    for g1, pg1 in enumerate(PG1):
        if pg1 == 0:
            continue
        SdelG1 += pg1 * calc_S(PdelG1_given_G1[..., g1], base)

    SdelG2 = 0
    for g2, pg2 in enumerate(PG2):
        if pg2 == 0:
            continue
        SdelG2 += pg2 * calc_S(PdelG2_given_G2[..., g2], base)
        
    return SG1G2, SdelG1, SdelG2

def calc_Sins(Pins, R, p_0 = None, base = 2):
    """Calculate the entropy of an insertion junction.

    Entropy of nucleotide choice in the Markov model is given by:
    H(p_0)(1-Pins[0]) + \sum_{l=2}Pins[l] \sum_{k=2}^l \sum_m R^{k-1}p_0(m) H(R[:, m])
    
    
    Parameters
    ----------
    Pins : ndarray
        Length distribution of the junction (e.g. PinsVD, PinsDJ, or PinsVJ)
    R : ndarray
        Markov transition matrix
    base : float
        Base used to set entropy units. Suggestions:
        base = 2 (entropy in bits, default)
        base = e (entropy in nats)
        base = 10 (entropy in dits)

    Returns
    -------
    Sins : float
        Entropy of Pins (entropy of length distribution)
    Smarkov : float
        Entropy of nucleotide choice in the Markov model.
        
    """
    Sins = calc_S(Pins, base)
    if p_0 is None: #Use steady state if no p_0 given
        p_0 = calc_steady_state_dist(R)
    cr_Pins = Pins[::-1].cumsum()[::-1]
    
    Smarkov = calc_S(p_0, base)*(1-Pins[0])
    p_k = p_0
    for k in range(2, len(Pins)):
        p_k = np.dot(R, p_k)
        Smarkov += cr_Pins[k] * np.sum([p_k[m]*calc_S(R[:, m], base) for m in range(4)])
    
    return Sins, Smarkov