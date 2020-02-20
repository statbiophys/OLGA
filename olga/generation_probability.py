#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Definition of classes GenerationProbabilityV(D)J to compute Pgen of a CDR3 seq.

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


This module defines three classes. The first, and parent of the other two,
GenerationProbability has wrapper methods for formatting and calling a pgen
computation for a nucleotide, amino acid, or regular expression CDR3 sequence.

GenerationProbabilityVDJ and GenerationProbabilityVJ each have a method, 
compute_CDR3_pgen, which implements the respective dynamic programming 
algorithm for either a VDJ or VJ generative recombination model.

In order to instantiate GenerationProbabilityV(D)J, instances of 
GenerativeModelV(D)J and GenomicDataV(D)J are needed. 
GenerationProbabilityV(D)J will inherit the processed parameters of the
provided generative model and genomic data through the 
PreprocessedParametersV(D)J classes.


Example
-------
>>> import olga.load_model as load_model
>>> import olga.generation_probability as pgen
>>> 
>>> params_file_name = './models/human_T_beta/model_params.txt'
>>> marginals_file_name = './models/human_T_beta/model_marginals.txt'
>>> V_anchor_pos_file ='./models/human_T_beta/V_gene_CDR3_anchors.csv'
>>> J_anchor_pos_file = './models/human_T_beta/J_gene_CDR3_anchors.csv'
>>> 
>>> genomic_data = load_model.GenomicDataVDJ()
>>> genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
>>> 
>>> generative_model = load_model.GenerativeModelVDJ()
>>> generative_model.load_and_process_igor_model(marginals_file_name)
>>> 
>>> pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)
>>> 
>>> pgen_model.compute_regex_CDR3_template_pgen('CASSAX{0,5}SARPEQFF')
6.846877804096558e-10
>>> 
>>> pgen_model.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF', 'TRBV30*01', 'TRBJ1-2*01')
1.203646865765782e-10
>>> 
>>> pgen_model.compute_nt_CDR3_pgen('TGTGCCAGTAGTATAACAACCCAGGGCTTGTACGAGCAGTACTTC')
3.9945642868171824e-14


@author: zacharysethna

"""
from __future__ import print_function, division
import numpy as np
import re
from olga.utils import nt2codon_rep, gene_to_num_str
from olga.preprocess_generative_model_and_data import PreprocessedParametersVDJ, PreprocessedParametersVJ

#Set input = raw_input for python 2
try:
    import __builtin__
    input = getattr(__builtin__, 'raw_input')
except (ImportError, AttributeError):
    pass

class GenerationProbability(object):
    """Class used to define Pgen functions and sequence formatting.

    This class is used to define three types of functions that are used by 
    both the VDJ pgen algorithm and the VJ pgen algorithm.
    
    The first type is functions that wrap around the core 'amino acid' 
    algorithms to allow for computing Pgen of regular expression, amino acid,
    and nucleotide CDR3 sequences (etc).
    
    The second type are functions that format some of the inputs (V/J mask,
    lists seqs for regular expressions) of the first type.
    
    The last group of functions are alignment/matching scripts that are used to 
    check how much of an 'amino acid' CDR3 is consistent with a given
    nucleotide sequence. These methods are used in both core algorithms when
    including the V and J contributions.

    Attributes
    ----------
    codons_dict : dict
        Dictionary, keyed by the allowed 'amino acid' symbols with the values 
        being lists of codons corresponding to the symbol.
    
    d_V_usage_mask : list of int
        Default V usage mask of indices of all productive V genes/alleles.                   
    V_mask_mapping : dict
        Dictionary mapping allowed keywords (strings corresponding to V gene 
        and V allele names) to the indices of the V alleles they refer to.

    d_J_usage_mask : list of int
        Default J usage mask of indices of all productive J genes/alleles.             
    J_mask_mapping : dict
        Dictionary mapping allowed keywords (strings corresponding to J gene 
        and J allele names) to the indices of the J alleles they refer to.

    Methods
    ----------
    compute_regex_CDR3_template_pgen(regex_seq, V_usage_mask_in = None, J_usage_mask_in = None, print_warnings = True, raise_overload_warning = True)
        Compute Pgen for all seqs consistent with regular expression regex_seq.
        
    compute_aa_CDR3_pgen(CDR3_seq, V_usage_mask_in = None, J_usage_mask_in = None, print_warnings = True)
        Compute Pgen for the amino acid sequence CDR3_seq
        
    compute_hamming_dist_1_pgen(CDR3_seq, V_usage_mask_in = None, J_usage_mask_in = None, print_warnings = True)
        Compute Pgen of all seqs hamming dist 1 (in amino acids) from CDR3_seq
        
    compute_nt_CDR3_pgen(CDR3_ntseq, V_usage_mask_in = None, J_usage_mask_in = None, print_warnings = True)
        Compute Pgen for the inframe nucleotide sequence CDR3_ntseq.
        
    compute_CDR3_pgen(CDR3_seq, V_usage_mask, J_usage_mask)
        Dummy function that is replaced in classes GenerationProbabilityV(D)J.
        The methods that replace it implement the different algorithms for 
        computing Pgen on a VDJ CDR3 sequence or a VJ CDR3 sequence.
        
    format_usage_masks(V_usage_mask_in, J_usage_mask_in, print_warnings = True)
        Format raw usage masks into lists of indices.
        
    list_seqs_from_regex(regex_seq, print_warnings = True, raise_overload_warning = True)
        List sequences that match regular expression template.
        
    max_nt_to_aa_alignment_left(CDR3_seq, ntseq)
        Find maximum match between CDR3_seq and ntseq from the left.
        
    max_nt_to_aa_alignment_right(CDR3_seq, ntseq)
        Find maximum match between CDR3_seq and ntseq from the right.
            
    """
    def __init__(self):
        """Initialize class GenerationProbability.
        
        Only define dummy attributes for this class. The children classes 
        GenerationProbabilityVDJ and GenerationProbabilityVJ will initialize
        the actual attributes.
        
        """
        
        self.codons_dict = None
        self.d_V_usage_mask = None
        self.V_mask_mapping = None
        
        self.d_J_usage_mask = None
        self.J_mask_mapping = None
        
    def compute_regex_CDR3_template_pgen(self, regex_seq, V_usage_mask_in = None, J_usage_mask_in = None, print_warnings = True, raise_overload_warning = True):
        """Compute Pgen for all seqs consistent with regular expression regex_seq.
    
        Computes Pgen for a (limited vocabulary) regular expression of CDR3 
        amino acid sequences, conditioned on the V genes/alleles indicated in 
        V_usage_mask_in and the J genes/alleles in J_usage_mask_in. Please note
        that this function will list out all the sequences that correspond to the 
        regular expression and then calculate the Pgen of each sequence in
        succession. THIS CAN BE SLOW. Consider defining a custom alphabet to 
        represent any undetermined amino acids as this will greatly speed up the 
        computations. For example, if the symbol ^ is defined as [AGR] in a custom 
        alphabet, then instead of running 
        compute_regex_CDR3_template_pgen('CASS[AGR]SARPEQFF', ppp),
        which will compute Pgen for 3 sequences, the single sequence 
        'CASS^SARPEQFF' can be considered. (Examples are TCRB sequences/model)
        
    
        Parameters
        ----------
        regex_seq : str
            The regular expression string that represents the CDR3 sequences to be 
            listed then their Pgens computed and summed.
        V_usage_mask_in : str or list
            An object to indicate which V alleles should be considered. The default
            input is None which returns the list of all productive V alleles.
        J_usage_mask_in : str or list
            An object to indicate which J alleles should be considered. The default
            input is None which returns the list of all productive J alleles.
        print_warnings : bool
            Determines whether warnings are printed or not. Default ON.
        raise_overload_warning : bool
            A flag to warn of more than 10000 seqs corresponding to the regex_seq
    
        Returns
        -------
        pgen : float
            The generation probability (Pgen) of the sequence
        
        Examples
        --------
        >>> generation_probability.compute_regex_CDR3_template_pgen('CASS[AGR]SARPEQFF')
        8.1090898050318022e-10
        >>> generation_probability.compute_regex_CDR3_template_pgen('CASSAX{0,5}SARPEQFF')
        6.8468778040965569e-10
            
        """
        
        V_usage_mask, J_usage_mask = self.format_usage_masks(V_usage_mask_in, J_usage_mask_in, print_warnings)
        
        CDR3_seqs = self.list_seqs_from_regex(regex_seq, print_warnings, raise_overload_warning)
        
        pgen = 0
        for  CDR3_seq in CDR3_seqs:
            if len(CDR3_seq) == 0:
                continue
            pgen += self.compute_CDR3_pgen(CDR3_seq, V_usage_mask, J_usage_mask)
    
        return pgen
    
    
    def compute_aa_CDR3_pgen(self, CDR3_seq, V_usage_mask_in = None, J_usage_mask_in = None, print_warnings = True):
        """Compute Pgen for the amino acid sequence CDR3_seq.
    
        Conditioned on the V genes/alleles indicated in V_usage_mask_in and the 
        J genes/alleles in J_usage_mask_in. (Examples are TCRB sequences/model)
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' -- the standard amino acids, 
            plus any custom symbols for an expanded codon alphabet (note the 
            standard ambiguous amino acids -- B, J, X, and Z -- are included by
            default).
        V_usage_mask_in : str or list
            An object to indicate which V alleles should be considered. The default
            input is None which returns the list of all productive V alleles.
        J_usage_mask_in : str or list
            An object to indicate which J alleles should be considered. The default
            input is None which returns the list of all productive J alleles.
        print_warnings : bool
            Determines whether warnings are printed or not. Default ON.
    
        Returns
        -------
        pgen : float
            The generation probability (Pgen) of the sequence
        
        Examples
        --------
        >>> generation_probability.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF')
        1.5756106696284584e-10
        >>> generation_probability.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF', 'TRBV30*01', 'TRBJ1-2*01')
        1.203646865765782e-10
        >>> generation_probability.compute_aa_CDR3_pgen('CAWXXXXXXXGYTF')
        7.8102586432014974e-05
            
        """
        if len(CDR3_seq) == 0:
            return  0
        for aa in CDR3_seq:
            if aa not in self.codons_dict.keys():
                #Check to make sure all symbols are accounted for
                if print_warnings:
                    print('Invalid amino acid CDR3 sequence --- unfamiliar symbol: ' + aa)
                return 0    
        
        V_usage_mask, J_usage_mask = self.format_usage_masks(V_usage_mask_in, J_usage_mask_in, print_warnings)
        
        return self.compute_CDR3_pgen(CDR3_seq, V_usage_mask, J_usage_mask)
    
    def compute_hamming_dist_1_pgen(self, CDR3_seq, V_usage_mask_in = None, J_usage_mask_in = None, print_warnings = True):
        """Compute Pgen of all seqs hamming dist 1 (in amino acids) from CDR3_seq.
    
        Please note that this function will list out all the 
        sequences that are hamming distance 1 from the base sequence and then 
        calculate the Pgen of each sequence in succession. THIS CAN BE SLOW 
        as it computes Pgen for L+1 sequences where L = len(CDR3_seq). (Examples 
        are TCRB sequences/model)
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of amino acids (ONLY the standard amino acids).
            Pgens for all sequences of hamming distance 1 (in amino acid sequence) 
            are summed.
        V_usage_mask_in : str or list
            An object to indicate which V alleles should be considered. The default
            input is None which returns the list of all productive V alleles.
        J_usage_mask_in : str or list
            An object to indicate which J alleles should be considered. The default
            input is None which returns the list of all productive J alleles.
        print_warnings : bool
            Determines whether warnings are printed or not. Default ON.
    
        Returns
        -------
        pgen : float
            The sum of generation probabilities (Pgens) of the sequences at most
            hamming distance 1 (in amino acids) from CDR3_seq.
            
        """
        
        
        #make sure that the symbol X is defined as the fully undetermined amino acid:
        #X ~ ACDEFGHIKLMNPQRSTVWY
        
        V_usage_mask, J_usage_mask = self.format_usage_masks(V_usage_mask_in, J_usage_mask_in, print_warnings)
        
        if len(CDR3_seq) == 0:
            return  0
        for aa in CDR3_seq:
            if aa not in 'ACDEFGHIKLMNPQRSTVWY':
                #Check to make sure all symbols are accounted for
                if print_warnings:
                    print('Invalid amino acid CDR3 sequence --- unfamiliar symbol: ' + aa)
                return 0
        tot_pgen = 0
        for i in range(len(CDR3_seq)):
            tot_pgen += self.compute_CDR3_pgen(CDR3_seq[:i] + 'X' + CDR3_seq[i+1:], V_usage_mask, J_usage_mask)
        tot_pgen += -(len(CDR3_seq) - 1)*self.compute_CDR3_pgen(CDR3_seq, V_usage_mask, J_usage_mask)
        return tot_pgen
    
    def compute_nt_CDR3_pgen(self, CDR3_ntseq, V_usage_mask_in = None, J_usage_mask_in = None, print_warnings = True):
        """Compute Pgen for the inframe nucleotide sequence CDR3_ntseq.
    
        Conditioned on the V genes/alleles indicated in V_usage_mask_in and the 
        J genes/alleles in J_usage_mask_in. (Examples are TCRB sequences/model)
    
        Parameters
        ----------
        CDR3_ntseq : str
            Inframe nucleotide sequence composed of ONLY A, C, G, or T (either 
            uppercase or lowercase).
        V_usage_mask_in : str or list
            An object to indicate which V alleles should be considered. The default
            input is None which returns the list of all productive V alleles.
        J_usage_mask_in : str or list
            An object to indicate which J alleles should be considered. The default
            input is None which returns the list of all productive J alleles.
        print_warnings : bool
            Determines whether warnings are printed or not. Default ON.
    
        Returns
        -------
        pgen : float64
            The generation probability (Pgen) of the sequence
        
        Examples
        --------
        >>> generation_probability.compute_nt_CDR3_pgen('TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC')
        3.2674893012379071e-12
        >>> generation_probability.compute_nt_CDR3_pgen('TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC', 'TRBV30*01', 'TRBJ1-2*01')
        2.3986503758867323e-12
            
        """
        
        if not len(CDR3_ntseq)%3 == 0:
            #Make sure sequence is inframe
            if print_warnings:
                print('Invalid nucleotide CDR3 sequence --- out of frame sequence')
            return 0
        elif len(CDR3_ntseq) == 0:
            return 0
        else:
            for nt in CDR3_ntseq:
                if nt not in 'ACGTacgt':
                    if print_warnings:
                        print('Invalid nucleotide CDR3 sequence --- unfamiliar nucleotide: ' + nt)
                    return 0
        
        V_usage_mask, J_usage_mask = self.format_usage_masks(V_usage_mask_in, J_usage_mask_in, print_warnings)
        
        return self.compute_CDR3_pgen(nt2codon_rep(CDR3_ntseq), V_usage_mask, J_usage_mask)
    
    def compute_CDR3_pgen(self, CDR3_seq, V_usage_mask, J_usage_mask):
        """Dummy function that is replaced in classes GenerationProbabilityV(D)J."""
        #Proxy for the actual function that will call either the VDJ algorithm
        #or the VJ algorithm
        pass
    
    #Formatting methods for the top level Pgen computation calls
    def format_usage_masks(self, V_usage_mask_in, J_usage_mask_in, print_warnings = True):
        """Format raw usage masks into lists of indices.
    
        Usage masks allows the Pgen computation to be conditioned on the V and J 
        gene/allele identities. The inputted masks are lists of strings, or a 
        single string, of the names of the genes or alleles to be conditioned on. 
        The default mask includes all productive V or J genes.
    
        Parameters
        ----------
        V_usage_mask_in : str or list
            An object to indicate which V alleles should be considered. The default
            input is None which returns the list of all productive V alleles.
        J_usage_mask_in : str or list
            An object to indicate which J alleles should be considered. The default
            input is None which returns the list of all productive J alleles.
        print_warnings : bool
            Determines whether warnings are printed or not. Default ON.
    
        Returns
        -------
        V_usage_mask : list of integers
            Indices of the V alleles to be considered in the Pgen computation
        J_usage_mask : list of integers
            Indices of the J alleles to be considered in the Pgen computation
        
        Examples
        --------
        >>> generation_probability.format_usage_masks('TRBV27*01','TRBJ1-1*01')
        ([34], [0])
        >>> generation_probability.format_usage_masks('TRBV27*01', '')
        ([34], [0, 1, 2, 3, 4, 7, 8, 9, 10, 11, 12, 13])
        >>> generation_probability.format_usage_masks(['TRBV27*01', 'TRBV13*01'], 'TRBJ1-1*01')
        ([34, 18], [0])
            
        """
        #Format the V usage mask
        if isinstance(V_usage_mask_in, str):
            V_usage_mask_in = [V_usage_mask_in]
        
        if V_usage_mask_in is None: #Default case, use all productive V genes with non-zero probability
            #V_usage_mask = [v for v, V in enumerate(ppp['cutV_genomic_CDR3_segs']) if len(V) > 0]
            V_usage_mask = self.d_V_usage_mask
        elif isinstance(V_usage_mask_in, list):
            e_V_usage_mask = set()
            for v in V_usage_mask_in:
                try:
                    e_V_usage_mask = e_V_usage_mask.union(self.V_mask_mapping[gene_to_num_str(v, 'V')])
                except:
                    if print_warnings:
                        print('Unfamiliar V gene/allele: ' + v)
                    pass
            if len(e_V_usage_mask) == 0:
                if print_warnings:
                    print('No recognized V genes/alleles. Using default V_usage_mask')
                V_usage_mask = self.d_V_usage_mask
            else:
                V_usage_mask = list(e_V_usage_mask)
        else:
            if print_warnings:
                print('Unfamiliar typed V usage mask: ' + str(V_usage_mask_in) + '. Using default V_usage_mask')
            V_usage_mask = self.d_V_usage_mask
                
                
        #Format the J usage mask
        if isinstance(J_usage_mask_in, str):
            J_usage_mask_in = [J_usage_mask_in]
            
        if J_usage_mask_in is None: #Default case, use all productive J genes with non-zero probability
            #J_usage_mask = [j for j, J in enumerate(ppp['cutJ_genomic_CDR3_segs']) if len(J) > 0]
            J_usage_mask = self.d_J_usage_mask
        elif isinstance(J_usage_mask_in, list):
            e_J_usage_mask = set()
            for j in J_usage_mask_in:
                try:
                    e_J_usage_mask = e_J_usage_mask.union(self.J_mask_mapping[gene_to_num_str(j, 'J')])
                except:
                    if print_warnings:
                        print('Unfamiliar J gene/allele: ' + j)
                    pass
            if len(e_J_usage_mask) == 0:
                if print_warnings:
                    print('No recognized J genes/alleles. Using default J_usage_mask')
                J_usage_mask = self.d_J_usage_mask
            else:
                J_usage_mask = list(e_J_usage_mask)
        else:
            if print_warnings:
                print('Unfamiliar typed J usage mask: ' + str(J_usage_mask_in) + '. Using default J_usage_mask')
            J_usage_mask = self.d_J_usage_mask
                
        return V_usage_mask, J_usage_mask
    
    def list_seqs_from_regex(self, regex_seq, print_warnings = True, raise_overload_warning = True):
        """List sequences that match regular expression template.
    
        This function parses a limited regular expression vocabulary, and 
        lists all the sequences consistent with the regular expression. Supported 
        regex syntax: [] and {}. Cannot have two {} in a row. Note we can't use 
        Kline star (*) as this is the symbol for a stop codon --- use {}.
    
        Parameters
        ----------
        regex_seq : str
            The regular expression string that represents the sequences to be 
            listed.
        print_warnings : bool
            Determines whether warnings are printed or not. Default ON.
        raise_overload_warning : bool
            A flag to warn of more than 10000 seqs corresponding to the regex_seq
    
        Returns
        -------
        CDR3_seqs : list
            A list of CDR3 sequences that correspond to the regex_seq
        
        Examples
        --------
        >>> generation_probability.list_seqs_from_regex('CASS[AGR]SARPEQFF')
        ['CASSGSARPEQFF', 'CASSRSARPEQFF', 'CASSASARPEQFF']
        >>> generation_probability.list_seqs_from_regex('CASSAX{0,5}SARPEQFF')
        ['CASSASARPEQFF',
         'CASSAXXXXSARPEQFF',
         'CASSAXXSARPEQFF',
         'CASSAXXXXXSARPEQFF',
         'CASSAXXXSARPEQFF',
         'CASSAXSARPEQFF']
            
        """
    
        aa_symbols = ''.join(self.codons_dict)
        
        default_max_reps = 40
    
        #Check to make sure that expression is of the right form/symbols
        
        
        #Identify bracket expressions
        bracket_ex = [x for x in re.findall('\[[' + aa_symbols + ']*?\]|\{\d+,{0,1}\d*\}', regex_seq)]
        
        
        split_seq  = re.split('\[[' + aa_symbols + ']*?\]|\{\d+,{0,1}\d*\}', regex_seq)
        #Check that all remaining characters are in the codon dict
        for aa in ''.join(split_seq):
            if aa not in aa_symbols:
                if print_warnings:
                    print('Unfamiliar symbol representing a codon:' + aa + ' --- check codon dictionary or the regex syntax')
                return []
        
        
        regex_list = [split_seq[i//2] if i%2 == 0 else bracket_ex[i//2] for i in range(len(bracket_ex) + len(split_seq)) if not (i%2 == 0 and len(split_seq[i//2]) ==  0)]
        
        max_num_seqs = 1
        for l, ex in enumerate(regex_list[::-1]):
            i = len(regex_list) - l - 1
            if ex[0] == '[': #bracket expression
                #check characters
                for aa in ex.strip('[]'):
                    if aa not in aa_symbols:
                        if print_warnings:
                            print('Unfamiliar symbol representing a codon:' + aa + ' --- check codon dictionary')
                        return []
                max_num_seqs *= len(ex) - 2
            elif ex[0] == '{': #curly bracket
                if i == 0:
                    if print_warnings:
                        print("Can't have {} expression at start of sequence")
                    return []
                elif isinstance(regex_list[i-1], list):
                    if print_warnings:
                        print("Two {} expressions in a row is not supported")
                    return []
                elif regex_list[i-1][0] == '[':
                    syms = regex_list[i-1].strip('[]')
                    regex_list[i-1] = ''                
                else:
                    syms = regex_list[i-1][-1]
                    regex_list[i-1] = regex_list[i-1][:-1]
                if ',' not in ex:
                    new_expression = [int(ex.strip('{}')), int(ex.strip('{}')), syms]
                    max_num_seqs *= len(syms)**new_expression[0]
                else:
                    try:
                        new_expression = [int(ex.strip('{}').split(',')[0]), int(ex.strip('{}').split(',')[1]), syms]
                    except ValueError: #No max limit --- use default
                        new_expression = [int(ex.strip('{}').split(',')[0]), default_max_reps, syms]
                    if new_expression[0] > new_expression[1]:
                        if print_warnings:
                            print('Check regex syntax --- should be {min,max}')
                        return []
                    max_num_seqs *= sum([len(syms)**n for n in range(new_expression[0], new_expression[1]+1)])/len(syms)
                #print new_expression
                regex_list[i] = new_expression
                
        if max_num_seqs > 10000 and raise_overload_warning:
            if print_warnings:
                answer = input('Warning large number of sequences (estimated ' + str(max_num_seqs) + ' seqs) match the regular expression. Possible memory and time issues. Continue? (y/n)')
                if not answer == 'y':
                    print('Canceling...')
                    return []
            else:
                return []
        #print regex_list
        CDR3_seqs = ['']
        for l, ex in enumerate(regex_list[::-1]):
            i = len(regex_list) - l - 1
            if isinstance(ex, list): #curly bracket case
                c_seqs = ['']
                f_seqs = []
                for j in range(ex[1] + 1):
                    if j in range(ex[0], ex[1]+1):
                        f_seqs += c_seqs
                    c_seqs = [aa + c_seq for aa in ex[2] for c_seq in c_seqs]
                CDR3_seqs = [f_seq + CDR3_seq for f_seq in f_seqs for CDR3_seq in CDR3_seqs]
            elif len(ex) == 0:
                pass
            elif ex[0] == '[': #square bracket case
                CDR3_seqs = [aa + CDR3_seq for aa in ex.strip('[]') for CDR3_seq in CDR3_seqs]
            else:
                CDR3_seqs = [ex + CDR3_seq for CDR3_seq in CDR3_seqs]
    
    
        return list(set(CDR3_seqs))
    
    # Alignment/Matching methods
    def max_nt_to_aa_alignment_left(self, CDR3_seq, ntseq):    
        """Find maximum match between CDR3_seq and ntseq from the left.
    
        This function returns the length of the maximum length nucleotide
        subsequence of ntseq contiguous from the left (or 5' end) that is 
        consistent with the 'amino acid' sequence CDR3_seq.
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' (single character symbols
            each corresponding to a collection of codons as given by codons_dict).
        ntseq : str
            Genomic (V locus) nucleotide sequence to match.
    
        Returns
        -------
        max_alignment : int
            Maximum length (in nucleotides) nucleotide sequence that matches the 
            CDR3 'amino acid' sequence.
        
        Example
        --------
        >>> generation_probability.max_nt_to_aa_alignment_left('CASSSEGAGGPSLRGHEQFF', 'TGTGCCAGCAGTTTATCGATA')
        13
            
        """
        
        max_alignment = 0
        if len(ntseq) == 0:
            return 0
        aa_aligned = True
        while aa_aligned:
            if ntseq[max_alignment:max_alignment+3] in self.codons_dict[CDR3_seq[max_alignment//3]]:
                max_alignment += 3
                if max_alignment//3 == len(CDR3_seq):
                    return max_alignment
            else:
                break
                aa_aligned = False
        last_codon = ntseq[max_alignment:max_alignment+3]
        codon_frag = ''
        for nt in last_codon:
            codon_frag += nt
            if codon_frag in self.sub_codons_left[CDR3_seq[max_alignment//3]]:
                max_alignment += 1
            else:
                break
        return max_alignment
        
    def max_nt_to_aa_alignment_right(self, CDR3_seq, ntseq):
        """Find maximum match between CDR3_seq and ntseq from the right.
    
        This function returns the length of the maximum length nucleotide
        subsequence of ntseq contiguous from the right (or 3' end) that is 
        consistent with the 'amino acid' sequence CDR3_seq
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' (single character symbols
            each corresponding to a collection of codons as given by codons_dict).
        ntseq : str
            Genomic (J locus) nucleotide sequence to match. 
    
        Returns
        -------
        max_alignment : int
            Maximum length (in nucleotides) nucleotide sequence that matches the 
            CDR3 'amino acid' sequence.
        
        Example
        --------
        >>> generation_probability.max_nt_to_aa_alignment_right('CASSSEGAGGPSLRGHEQFF', 'TTCATGAACACTGAAGCTTTCTTT')
        6
            
        """
        r_CDR3_seq = CDR3_seq[::-1] #reverse CDR3_seq
        r_ntseq = ntseq[::-1] #reverse ntseq
        max_alignment = 0
        if len(ntseq) == 0:
            return 0
        aa_aligned = True
        while aa_aligned:
            if r_ntseq[max_alignment:max_alignment+3][::-1] in self.codons_dict[r_CDR3_seq[max_alignment//3]]:
                max_alignment += 3
                if max_alignment//3 == len(CDR3_seq):
                    return max_alignment
            else:
                break
                aa_aligned = False
        r_last_codon = r_ntseq[max_alignment:max_alignment+3]
        codon_frag = ''
        for nt in r_last_codon:
            codon_frag = nt + codon_frag
            if codon_frag in self.sub_codons_right[r_CDR3_seq[max_alignment//3]]:
                max_alignment += 1
            else:
                break
        return max_alignment

#%%    
class GenerationProbabilityVDJ(GenerationProbability, PreprocessedParametersVDJ):
    """Class used to compute the Pgen of CDR3 sequences from a  VDJ model.
    
    All of the attributes of GenerationProbabilityVDJ are inherited from the
    class PreprocessedParametersVDJ. The methods of the class are used to 
    compute the Pgen of an 'amino acid' sequence from a VDJ generative model.

    Attributes
    ----------
    codons_dict : dict
        Dictionary, keyed by the allowed 'amino acid' symbols with the values 
        being lists of codons corresponding to the symbol.        
    sub_codons_left : dict
        Dictionary of the 1 and 2 nucleotide prefixes (read from 5') for 
        each codon in an 'amino acid' grouping        
    sub_codons_right : dict
        Dictionary of the 1 and 2 nucleotide suffixes (read from 5') for 
        each codon in an 'amino acid' grouping
        
    V_allele_names : list of strings
        List of V allele names in genomic_data        
    d_V_usage_mask : list of int
        Default V usage mask of indices of all productive V genes/alleles.        
    V_mask_mapping : dict
        Dictionary mapping allowed keywords (strings corresponding to V gene 
        and V allele names) to the indices of the V alleles they refer to.        
    cutV_genomic_CDR3_segs : list of strings
        List of the V germline nucleotide sequences, trimmed to begin at the
        CDR3 region (includes the conserved C residue) with the maximum number
        of reverse complementary palindromic insertions appended.        
    PVdelV_nt_pos_vec : list of ndarrays
        For each V allele, format P(V)*P(delV|V) into the correct form for a Pi 
        array or V_{x_1}. This is only done for the first and last position in 
        each codon.         
    PVdelV_2nd_nt_pos_per_aa_vec : list of dicts
        For each V allele, and each 'amino acid', format P(V)*P(delV|V) for 
        positions in the middle of a codon into the correct form for a Pi array 
        or V_{x_1} given the 'amino acid'.
        
    D_allele_names : list of strings
        List of D allele names in genomic_data        
    cutD_genomic_CDR3_segs : list of strings
        List of the D germline nucleotide sequences, with the maximum number
        of reverse complementary palindromic insertions appended to both ends.        
    PD_nt_pos_vec : list of ndarrays
        For each D allele, format P(delDl, delDr|D) into the correct form for a 
        Pi array as if each position were the first in a codon.    
    PD_2nd_nt_pos_per_aa_vec : list of dicts
        For each D allele, and each 'amino acid', format P(delDl, delDr|D) for 
        positions in the middle of a codon into the correct form for a Pi array 
        as if each position were the middle of a codon corresponding to the 
        'amino acid'.        
    min_delDl_given_DdelDr : list of lists
        minimum delDl for each delDr, D combination.        
    max_delDl_given_DdelDr : list of lists
        maximum delDl for each delDr, D combination.
    zeroD_given_D : list of floats
        The probability that a given D allele is fully deleted away.        
    PdelDldelDr_given_D : ndarray
        Joint probability distribution of the D deletions given the D allele,
        i.e. P(delDl, delDr |D)
        
    J_allele_names : list of strings
        List of J allele names in genomic_data        
    d_J_usage_mask : list of int
        Default J usage mask of indices of all productive J genes/alleles.        
    J_mask_mapping : dict
        Dictionary mapping allowed keywords (strings corresponding to V gene 
        and V allele names) to the indices of the V alleles they refer to.        
    cutJ_genomic_CDR3_segs : list of strings
        List of the J germline nucleotide sequences, trimmed to end at the
        CDR3 region (includes the conserved F or W residue) with the maximum 
        number of reverse complementary palindromic insertions appended.        
    PJdelJ_nt_pos_vec : list
        For each J allele, format P(J)*P(delJ|J) into the correct form for a Pi 
        array or J(D)^{x_4}. This is only done for the first and last position 
        in each codon.        
    PJdelJ_2nd_nt_pos_per_aa_vec : list
        For each J allele, and each 'amino acid', format P(J)*P(delJ|J) for 
        positions in the middle of a codon into the correct form for a Pi array 
        or J(D)^{x_4} given the 'amino acid'.
    PD_given_J : ndarray
        Probability distribution of D conditioned on J, i.e. P(D|J). 
         
    PinsVD : ndarray
        Probability distribution of the VD (N1) insertion sequence length        
    PinsDJ : ndarray
        Probability distribution of the DJ (N2) insertion sequence length        
    Rvd : ndarray
        Markov transition matrix for the VD insertion junction.
    Rdj : ndarray
        Markov transition matrix for the DJ insertion junction.        
    first_nt_bias_insVD : ndarray
        (4,) array of the probability distribution of the indentity of the 
        first nucleotide insertion for the VD junction.        
    zero_nt_bias_insVD : ndarray
        (4,) array of the probability distribution of the indentity of the 
        the nucleotide BEFORE the VD insertion.
        zero_nt_bias_insVD = Rvd^{-1}first_nt_bias_insVD        
    first_nt_bias_insDJ : ndarray
        (4,) array of the probability distribution of the indentity of the 
        first nucleotide insertion for the DJ junction.        
    zero_nt_bias_insDJ : ndarray
        (4,) array of the probability distribution of the indentity of the 
        the nucleotide BEFORE the DJ insertion. Note, as the Markov model
        at the DJ junction goes 3' to 5' this is the position AFTER the
        insertions reading left to right.
        zero_nt_bias_insVD = Rvd^{-1}first_nt_bias_insVD        
    Tvd, Svd, Dvd, lTvd, lDvd : dicts
        Dictionaries (by 'amino acid') of insertion transfer matrices 
        ((4, 4) ndarrays) for the VD junction.        
    Tdj, Sdj, Ddj, rTdj, rDdj : dicts
        Dictionaries (by 'amino acid') of insertion transfer matrices 
        ((4, 4) ndarrays) for the DJ junction.
        
    """
    def __init__(self, generative_model, genomic_data, alphabet_file = None):
        """Initialize GenerationProbabilityVDJ.
        
        This intialization inherits all of the attributes of 
        PreprocessedParametersVDJ (which include all of the processed 
        parameters needed for Pgen computation) and the methods of 
        GenerationProbability which include some wrappers/formatting of
        sequences to make Pgen computation of nucleotide and regular expression
        sequences easier (etc).
        
        Parameters
        ----------
        generative_model : GenerativeModelVDJ
            VDJ generative model class containing the model parameters.            
        genomic_data : GenomicDataVDJ
            VDJ genomic data class containing the V, D, and J germline 
            sequences and info.            
        alphabet_file : str, optional
            File name (full pathing from current directory) for a custom alphabet
            definition. If no file is provided, the default alphabet is used, i.e. 
            standard amino acids, undetermined amino acids (B, J, X, and Z), and
            single codon symbols.
        
        """
        GenerationProbability.__init__(self)
        PreprocessedParametersVDJ.__init__(self, generative_model, genomic_data, alphabet_file)
        
    def compute_CDR3_pgen(self, CDR3_seq, V_usage_mask, J_usage_mask):
        """Compute Pgen for CDR3 'amino acid' sequence CDR3_seq from VDJ model.
    
        Conditioned on the already formatted V genes/alleles indicated in 
        V_usage_mask and the J genes/alleles in J_usage_mask. 
        (Examples are TCRB sequences/model)
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' (single character symbols
            each corresponding to a collection of codons as given by codons_dict).
        V_usage_mask : list
            Indices of the V alleles to be considered in the Pgen computation
        J_usage_mask : list
            Indices of the J alleles to be considered in the Pgen computation
    
        Returns
        -------
        pgen : float
            The generation probability (Pgen) of the sequence
        
        Examples
        --------
        >>> compute_CDR3_pgen('CAWSVAPDRGGYTF', ppp, [42], [1])
        1.203646865765782e-10
        >>> compute_CDR3_pgen(nt2codon_rep('TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC'), ppp, [42], [1])
        2.3986503758867323e-12
        >>> compute_CDR3_pgen('\xbb\x96\xab\xb8\x8e\xb6\xa5\x92\xa8\xba\x9a\x93\x94\x9f', ppp, [42], [1])
        2.3986503758867323e-12
        
        """
        
        
        #Genomic V alignment/matching (contribution from P(V, delV)), return Pi_V
        Pi_V, max_V_align = self.compute_Pi_V(CDR3_seq, V_usage_mask)
        
        #Include VD insertions (Rvd and PinsVD) to get the total contribution from the left (3') side. Return Pi_L
        Pi_L = self.compute_Pi_L(CDR3_seq, Pi_V, max_V_align)
        
        #Genomic J alignment/matching (contribution from P(D, J, delJ)), return Pi_J_given_D
        Pi_J_given_D, max_J_align = self.compute_Pi_J_given_D(CDR3_seq, J_usage_mask)
        
        #Include DJ insertions (Rdj and PinsDJ), return Pi_JinsDJ_given_D
        Pi_JinsDJ_given_D = self.compute_Pi_JinsDJ_given_D(CDR3_seq, Pi_J_given_D, max_J_align)
        
        #Include D genomic contribution (P(delDl, delDr | D)) to complete the contribution from the right (5') side. Return Pi_R
        Pi_R = self.compute_Pi_R(CDR3_seq, Pi_JinsDJ_given_D)
        
        pgen = 0
        
        #zip Pi_L and Pi_R together to get total pgen
        for pos in range(len(CDR3_seq)*3 - 1):
            pgen += np.dot(Pi_L[:, pos], Pi_R[:, pos+1])
            
        return pgen
    
    #Genomic V alignment/matching (contribution from P(V, delV)), return Pi_V
    def compute_Pi_V(self, CDR3_seq, V_usage_mask):
        """Compute Pi_V.
    
        This function returns the Pi array from the model factors of the V genomic 
        contributions, P(V)*P(delV|V). This corresponds to V_{x_1}. 
        
        For clarity in parsing the algorithm implementation, we include which 
        instance attributes are used in the method as 'parameters.'
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' (single character symbols
            each corresponding to a collection of codons as given by codons_dict).            
        V_usage_mask : list
            Indices of the V alleles to be considered in the Pgen computation
            
        self.cutV_genomic_CDR3_segs : list of strings
            List of all the V genomic nucleotide sequences trimmed to begin at the 
            conserved C residue and with the maximum number of palindromic 
            insertions appended.
        self.PVdelV_nt_pos_vec : list of ndarrays
            For each V allele, format P(V)*P(delV|V) into the correct form for 
            a Pi array or V_{x_1}. This is only done for the first and last 
            position in each codon.    
        self.PVdelV_2nd_nt_pos_per_aa_vec : list of dicts
            For each V allele, and each 'amino acid', format P(V)*P(delV|V) for 
            positions in the middle of a codon into the correct form for a Pi 
            array or V_{x_1} given the 'amino acid'.
    
        Returns
        -------
        Pi_V : ndarray
            (4, 3L) array corresponding to V_{x_1}.
        max_V_align: int
            Maximum alignment of the CDR3_seq to any genomic V allele allowed by
            V_usage_mask.
            
        """
        #Note, the cutV_genomic_CDR3_segs INCLUDE the palindromic insertions and thus are max_palindrome nts longer than the template.
        #furthermore, the genomic sequence should be pruned to start at the conserved C
        
        Pi_V = np.zeros((4, len(CDR3_seq)*3)) #Holds the aggregate weight for each nt possiblity and position
        alignment_lengths = []
        for V_in in V_usage_mask:
            try:
                cutV_gen_seg = self.cutV_genomic_CDR3_segs[V_in]
            except IndexError:
                print('Check provided V usage mask. Contains indicies out of allowed range.')
                continue
            current_alignment_length = self.max_nt_to_aa_alignment_left(CDR3_seq, cutV_gen_seg)
            alignment_lengths += [current_alignment_length]
            current_Pi_V = np.zeros((4, len(CDR3_seq)*3))
            
            if current_alignment_length > 0:
                #For first and last nt in a codon use PVdelV_nt_pos_vec
                current_Pi_V[:, :current_alignment_length] = self.PVdelV_nt_pos_vec[V_in][:, :current_alignment_length]
                for pos in range(1, current_alignment_length, 3): #for middle nt use PVdelV_2nd_nt_pos_per_aa_vec
                    current_Pi_V[:, pos] = self.PVdelV_2nd_nt_pos_per_aa_vec[V_in][CDR3_seq[pos//3]][:, pos]
                Pi_V[:, :current_alignment_length] += current_Pi_V[:, :current_alignment_length]
        
        return Pi_V, max(alignment_lengths)
    
    #Include VD insertions (Rvd and PinsVD) to get the total contribution from the left (5') side. Return Pi_L
    def compute_Pi_L(self, CDR3_seq, Pi_V, max_V_align):
        """Compute Pi_L.
    
        This function returns the Pi array from the model factors of the V genomic 
        contributions, P(V)*P(delV|V), and the VD (N1) insertions,
        first_nt_bias_insVD(m_1)PinsVD(\ell_{VD})\prod_{i=2}^{\ell_{VD}}Rvd(m_i|m_{i-1}). 
        This corresponds to V_{x_1}{M^{x_1}}_{x_2}.
        
        For clarity in parsing the algorithm implementation, we include which 
        instance attributes are used in the method as 'parameters.'
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' (single character symbols
            each corresponding to a collection of codons as given by codons_dict).
        Pi_V : ndarray
            (4, 3L) array corresponding to V_{x_1}.
        max_V_align : int
            Maximum alignment of the CDR3_seq to any genomic V allele allowed by
            V_usage_mask.
            
        self.PinsVD : ndarray
            Probability distribution of the VD (N1) insertion sequence length            
        self.first_nt_bias_insVD : ndarray
            (4,) array of the probability distribution of the indentity of the 
            first nucleotide insertion for the VD junction.        
        self.zero_nt_bias_insVD : ndarray
            (4,) array of the probability distribution of the indentity of the 
            the nucleotide BEFORE the VD insertion.
            zero_nt_bias_insVD = Rvd^{-1}first_nt_bias_insVD       
        self.Tvd : dict
            Dictionary of full codon transfer matrices ((4, 4) ndarrays) by 
            'amino acid'.
        self.Svd : dict
            Dictionary of transfer matrices ((4, 4) ndarrays) by 'amino acid' for 
            the VD insertion ending in the first position.
        self.Dvd : dict
            Dictionary of transfer matrices ((4, 4) ndarrays) by 'amino acid' for 
            the VD insertion ending in the second position.
        self.lTvd : dict
            Dictionary of transfer matrices ((4, 4) ndarrays) by 'amino acid' for 
            the VD insertion starting in the first position.
        self.lDvd : dict
            Dictionary of transfer matrices ((4, 4) ndarrays) by 'amino acid' for
            VD insertion starting in the first position and ending in the second 
            position of the same codon.
    
        Returns
        -------
        Pi_L : ndarray
            (4, 3L) array corresponding to V_{x_1}{M^{x_1}}_{x_2}.
            
        """
        #max_insertions = 30 #len(PinsVD) - 1 should zeropad the last few spots
        max_insertions = len(self.PinsVD) - 1
        
        Pi_L = np.zeros((4, len(CDR3_seq)*3))
        
        #start position is first nt in a codon
        for init_pos in range(0, max_V_align, 3):
            #Zero insertions
            Pi_L[:, init_pos] += self.PinsVD[0]*Pi_V[:, init_pos]
            
            #One insertion
            Pi_L[:, init_pos+1] += self.PinsVD[1]*np.dot(self.lDvd[CDR3_seq[init_pos//3]], Pi_V[:, init_pos])
    
            #Two insertions and compute the base nt vec for the standard loop        
            current_base_nt_vec = np.dot(self.lTvd[CDR3_seq[init_pos//3]], Pi_V[:, init_pos])
            Pi_L[0, init_pos+2] += self.PinsVD[2]*np.sum(current_base_nt_vec)
            
            base_ins = 2
            
            #Loop over all other insertions using base_nt_vec
            for aa in CDR3_seq[init_pos//3 + 1: init_pos//3 + max_insertions//3]:
                Pi_L[:, init_pos+base_ins+1] += self.PinsVD[base_ins + 1]*np.dot(self.Svd[aa], current_base_nt_vec)
                Pi_L[:, init_pos+base_ins+2] += self.PinsVD[base_ins + 2]*np.dot(self.Dvd[aa], current_base_nt_vec)
                current_base_nt_vec = np.dot(self.Tvd[aa], current_base_nt_vec)
                Pi_L[0, init_pos+base_ins+3] += self.PinsVD[base_ins + 3]*np.sum(current_base_nt_vec)
                base_ins +=3
            
        
        #start position is second nt in a codon
        for init_pos in range(1, max_V_align, 3):
            #Zero insertions
            Pi_L[:, init_pos] += self.PinsVD[0]*Pi_V[:, init_pos]
            #One insertion --- we first compute our p vec by pairwise mult with the ss distr
            current_base_nt_vec = np.multiply(Pi_V[:, init_pos], self.first_nt_bias_insVD)
            Pi_L[0, init_pos+1] += self.PinsVD[1]*np.sum(current_base_nt_vec)
            base_ins = 1
            
            #Loop over all other insertions using base_nt_vec
            for aa in CDR3_seq[init_pos//3 + 1: init_pos//3 + max_insertions//3]:
                Pi_L[:, init_pos+base_ins+1] += self.PinsVD[base_ins + 1]*np.dot(self.Svd[aa], current_base_nt_vec)
                Pi_L[:, init_pos+base_ins+2] += self.PinsVD[base_ins + 2]*np.dot(self.Dvd[aa], current_base_nt_vec)
                current_base_nt_vec = np.dot(self.Tvd[aa], current_base_nt_vec)
                Pi_L[0, init_pos+base_ins+3] += self.PinsVD[base_ins + 3]*np.sum(current_base_nt_vec)
                base_ins +=3
            
        #start position is last nt in a codon   
        for init_pos in range(2, max_V_align, 3):
            
            #Zero insertions
            Pi_L[0, init_pos] += self.PinsVD[0]*Pi_V[0, init_pos]
            #current_base_nt_vec = first_nt_bias_insVD*Pi_V[0, init_pos] #Okay for steady state
            current_base_nt_vec = self.zero_nt_bias_insVD*Pi_V[0, init_pos]
            base_ins = 0
            
            #Loop over all other insertions using base_nt_vec
            for aa in CDR3_seq[init_pos//3 + 1: init_pos//3 + max_insertions//3]:
                Pi_L[:, init_pos+base_ins+1] += self.PinsVD[base_ins + 1]*np.dot(self.Svd[aa], current_base_nt_vec)
                Pi_L[:, init_pos+base_ins+2] += self.PinsVD[base_ins + 2]*np.dot(self.Dvd[aa], current_base_nt_vec)
                current_base_nt_vec = np.dot(self.Tvd[aa], current_base_nt_vec)
                Pi_L[0, init_pos+base_ins+3] += self.PinsVD[base_ins + 3]*np.sum(current_base_nt_vec)
                base_ins +=3
    
         
        return Pi_L
    
    #Genomic J alignment/matching (contribution from P(D, J, delJ)), return Pi_J_given_D
    def compute_Pi_J_given_D(self, CDR3_seq, J_usage_mask):
        """Compute Pi_J conditioned on D.
    
        This function returns the Pi array from the model factors of the D and J 
        genomic contributions, P(D, J)*P(delJ|J) = P(D|J)P(J)P(delJ|J). This 
        corresponds to J(D)^{x_4}.
        
        For clarity in parsing the algorithm implementation, we include which 
        instance attributes are used in the method as 'parameters.'
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' (single character symbols
            each corresponding to a collection of codons as given by codons_dict).   
        J_usage_mask : list
            Indices of the J alleles to be considered in the Pgen computation.
            
        self.cutJ_genomic_CDR3_segs : list
            List of all the J genomic nucleotide sequences trimmed to begin at the 
            conserved 3' residue (F/W) and with the maximum number of palindromic 
            insertions appended.
        self.PD_given_J : ndarray
            Probability distribution of D conditioned on J, i.e. P(D|J).
        self.PJdelJ_nt_pos_vec : list of ndarrays
            For each J allele, format P(J)*P(delJ|J) into the correct form for 
            a Pi array or J(D)^{x_4}. This is only done for the first and last 
            position in each codon.    
        self.PJdelJ_2nd_nt_pos_per_aa_vec : list of dicts
            For each J allele, and each 'amino acid', format P(J)*P(delJ|J) for 
            positions in the middle of a codon into the correct form for a Pi 
            array or J(D)^{x_4} given the 'amino acid'.
        
        Returns
        -------
        Pi_J_given_D : list
            List of (4, 3L) ndarrays corresponding to J(D)^{x_4}.
        max_J_align: int
            Maximum alignment of the CDR3_seq to any genomic J allele allowed by
            J_usage_mask.
        
            
        """
        
        #Note, the cutJ_genomic_CDR3_segs INCLUDE the palindromic insertions and thus are max_palindrome nts longer than the template.
        #furthermore, the genomic sequence should be pruned to start at a conserved region on the J side
        num_D_genes = self.PD_given_J.shape[0]
        Pi_J_given_D = [np.zeros((4, len(CDR3_seq)*3)) for i in range(num_D_genes)] #Holds the aggregate weight for each nt possiblity and position
        alignment_lengths = []
        for J_in in J_usage_mask:
            try:
                cutJ_gen_seg = self.cutJ_genomic_CDR3_segs[J_in]
            except IndexError:
                print('Check provided V usage mask. Contains indicies out of allowed range.')
                continue
            current_alignment_length = self.max_nt_to_aa_alignment_right(CDR3_seq, cutJ_gen_seg)
            alignment_lengths += [current_alignment_length]
            current_Pi_J = np.zeros((4, len(CDR3_seq)*3))
    
            if current_alignment_length > 0:
                #For first and last nt in a codon use PJdelJ_nt_pos_vec
                current_Pi_J[:, -current_alignment_length:] = self.PJdelJ_nt_pos_vec[J_in][:, -current_alignment_length:]          
                for pos in range(-2, -current_alignment_length-1, -3): #for middle nt use PJdelJ_2nd_nt_pos_per_aa_vec
                    current_Pi_J[:, pos] = self.PJdelJ_2nd_nt_pos_per_aa_vec[J_in][CDR3_seq[pos//3]][:, pos]
    
            for D_in, pd_given_j in enumerate(self.PD_given_J[:, J_in]):
                Pi_J_given_D[D_in][:, -current_alignment_length:] += pd_given_j*current_Pi_J[:, -current_alignment_length:]
        
        return Pi_J_given_D, max(alignment_lengths)
    
    #Include DJ insertions (Rdj and PinsDJ), return Pi_JinsDJ_given_D
    def compute_Pi_JinsDJ_given_D(self, CDR3_seq, Pi_J_given_D, max_J_align):
        """Compute Pi_JinsDJ conditioned on D.
    
        This function returns the Pi array from the model factors of the J genomic 
        contributions, P(D,J)*P(delJ|J), and the DJ (N2) insertions,
        first_nt_bias_insDJ(n_1)PinsDJ(\ell_{DJ})\prod_{i=2}^{\ell_{DJ}}Rdj(n_i|n_{i-1}) 
        conditioned on D identity. This corresponds to {N^{x_3}}_{x_4}J(D)^{x_4}.
        
        For clarity in parsing the algorithm implementation, we include which 
        instance attributes are used in the method as 'parameters.'
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' (single character symbols
            each corresponding to a collection of codons as given by codons_dict).
        Pi_J_given_D : ndarray
            List of (4, 3L) ndarrays corresponding to J(D)^{x_4}.
        max_J_align : int
            Maximum alignment of the CDR3_seq to any genomic J allele allowed by
            J_usage_mask.
            
        self.PinsDJ : ndarray
            Probability distribution of the DJ (N2) insertion sequence length    
        self.first_nt_bias_insDJ : ndarray
            (4,) array of the probability distribution of the indentity of the 
            first nucleotide insertion for the DJ junction.        
        self.zero_nt_bias_insDJ : ndarray
            (4,) array of the probability distribution of the indentity of the 
            the nucleotide BEFORE the DJ insertion. Note, as the Markov model
            at the DJ junction goes 3' to 5' this is the position AFTER the
            insertions reading left to right.
        self.Tdj : dict
            Dictionary of full codon transfer matrices ((4, 4) ndarrays) by 
            'amino acid'.
        self.Sdj : dict
            Dictionary of transfer matrices ((4, 4) ndarrays) by 'amino acid' for 
            the DJ insertion ending in the first position.
        self.Ddj : dict
            Dictionary of transfer matrices ((4, 4) ndarrays) by 'amino acid' for 
            the VD insertion ending in the second position.
        self.rTdj : dict
            Dictionary of transfer matrices ((4, 4) ndarrays) by 'amino acid' for 
            the DJ insertion starting in the first position.
        self.rDdj : dict
            Dictionary of transfer matrices ((4, 4) ndarrays) by 'amino acid' for
            DJ insertion starting in the first position and ending in the second 
            position of the same codon.
    
        Returns
        -------
        Pi_JinsDJ_given_D : list
            List of (4, 3L) ndarrays corresponding to {N^{x_3}}_{x_4}J(D)^{x_4}.
            
        """
        #max_insertions = 30 #len(PinsVD) - 1 should zeropad the last few spots
        max_insertions = len(self.PinsDJ) - 1
        
        
        Pi_JinsDJ_given_D = [np.zeros((4, len(CDR3_seq)*3)) for i in range(len(Pi_J_given_D))]
        
        for D_in in range(len(Pi_J_given_D)):
            #start position is first nt in a codon
            for init_pos in range(-1, -(max_J_align+1), -3):
                #Zero insertions
                Pi_JinsDJ_given_D[D_in][:, init_pos] += self.PinsDJ[0]*Pi_J_given_D[D_in][:, init_pos]
                
                #One insertion
                Pi_JinsDJ_given_D[D_in][:, init_pos-1] += self.PinsDJ[1]*np.dot(self.rDdj[CDR3_seq[init_pos//3]], Pi_J_given_D[D_in][:, init_pos])
        
                #Two insertions and compute the base nt vec for the standard loop        
                current_base_nt_vec = np.dot(self.rTdj[CDR3_seq[init_pos//3]], Pi_J_given_D[D_in][:, init_pos])
                Pi_JinsDJ_given_D[D_in][0, init_pos-2] += self.PinsDJ[2]*np.sum(current_base_nt_vec)
                
                base_ins = 2
                
                #Loop over all other insertions using base_nt_vec
                for aa in CDR3_seq[init_pos//3 - 1: init_pos//3 - max_insertions//3:-1]:
                    Pi_JinsDJ_given_D[D_in][:, init_pos-base_ins-1] += self.PinsDJ[base_ins + 1]*np.dot(self.Sdj[aa], current_base_nt_vec)
                    Pi_JinsDJ_given_D[D_in][:, init_pos-base_ins-2] += self.PinsDJ[base_ins + 2]*np.dot(self.Ddj[aa], current_base_nt_vec)
                    current_base_nt_vec = np.dot(self.Tdj[aa], current_base_nt_vec)
                    Pi_JinsDJ_given_D[D_in][0, init_pos-base_ins-3] += self.PinsDJ[base_ins + 3]*np.sum(current_base_nt_vec)
                    base_ins +=3
                
            
            #start position is second nt in a codon
            for init_pos in range(-2, -(max_J_align+1), -3):
                #Zero insertions
                Pi_JinsDJ_given_D[D_in][:, init_pos] += self.PinsDJ[0]*Pi_J_given_D[D_in][:, init_pos]
                #One insertion --- we first compute our p vec by pairwise mult with the ss distr
                current_base_nt_vec = np.multiply(Pi_J_given_D[D_in][:, init_pos], self.first_nt_bias_insDJ)
                Pi_JinsDJ_given_D[D_in][0, init_pos-1] += self.PinsDJ[1]*np.sum(current_base_nt_vec)
                base_ins = 1
                
                #Loop over all other insertions using base_nt_vec
                for aa in CDR3_seq[init_pos//3 - 1: init_pos//3 - max_insertions//3:-1]:
                    Pi_JinsDJ_given_D[D_in][:, init_pos-base_ins-1] += self.PinsDJ[base_ins + 1]*np.dot(self.Sdj[aa], current_base_nt_vec)
                    Pi_JinsDJ_given_D[D_in][:, init_pos-base_ins-2] += self.PinsDJ[base_ins + 2]*np.dot(self.Ddj[aa], current_base_nt_vec)
                    current_base_nt_vec = np.dot(self.Tdj[aa], current_base_nt_vec)
                    Pi_JinsDJ_given_D[D_in][0, init_pos-base_ins-3] += self.PinsDJ[base_ins + 3]*np.sum(current_base_nt_vec)
                    base_ins +=3
                
            #start position is last nt in a codon   
            for init_pos in range(-3, -(max_J_align+1), -3):
                
                #Zero insertions
                Pi_JinsDJ_given_D[D_in][0, init_pos] += self.PinsDJ[0]*Pi_J_given_D[D_in][0, init_pos]
                #current_base_nt_vec = first_nt_bias_insDJ*Pi_J_given_D[D_in][0, init_pos] #Okay for steady state
                current_base_nt_vec = self.zero_nt_bias_insDJ*Pi_J_given_D[D_in][0, init_pos]
                base_ins = 0
                
                #Loop over all other insertions using base_nt_vec
                for aa in CDR3_seq[init_pos//3 - 1: init_pos//3 - max_insertions//3:-1]:
                    Pi_JinsDJ_given_D[D_in][:, init_pos-base_ins-1] += self.PinsDJ[base_ins + 1]*np.dot(self.Sdj[aa], current_base_nt_vec)
                    Pi_JinsDJ_given_D[D_in][:, init_pos-base_ins-2] += self.PinsDJ[base_ins + 2]*np.dot(self.Ddj[aa], current_base_nt_vec)
                    current_base_nt_vec = np.dot(self.Tdj[aa], current_base_nt_vec)
                    Pi_JinsDJ_given_D[D_in][0, init_pos-base_ins-3] += self.PinsDJ[base_ins + 3]*np.sum(current_base_nt_vec)
                    base_ins +=3
    
         
        return Pi_JinsDJ_given_D
    
    #Include D genomic contribution (P(delDl, delDr | D)) to complete the contribution from the right (5') side. Return Pi_R
    def compute_Pi_R(self, CDR3_seq, Pi_JinsDJ_given_D):
        """Compute Pi_R.
    
        This function returns the Pi array from the model factors of the D and J 
        genomic contributions, P(D, J)*P(delJ|J)P(delDl, delDr |D) and
        the DJ (N2) insertions,
        first_nt_bias_insDJ(n_1)PinsDJ(\ell_{DJ})\prod_{i=2}^{\ell_{DJ}}Rdj(n_i|n_{i-1}).
        This corresponds to \sum_D {D^{x_2}}_{x_3}{N^{x_3}}_{x_4}J(D)^{x_4}.
        
        For clarity in parsing the algorithm implementation, we include which 
        instance attributes are used in the method as 'parameters.'
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' (single character symbols
            each corresponding to a collection of codons as given by codons_dict).
        Pi_JinsDJ_given_D : list
            List of (4, 3L) ndarrays corresponding to {N^{x_3}}_{x_4}J(D)^{x_4}.
            
        self.cutD_genomic_CDR3_segs : list of strings
            List of all the D genomic nucleotide sequences with the maximum number 
            of palindromic insertions appended on both ends.
        self.PD_given_J : ndarray
            Probability distribution of D conditioned on J, i.e. P(D|J).
        self.PD_nt_pos_vec : list of ndarrays
            For each D allele, format P(delDl, delDr|D) into the correct form 
            for a Pi array as if each position were the first in a codon.
        self.PD_2nd_nt_pos_per_aa_vec : list of dicts
            For each D allele, and each 'amino acid', format P(delDl, delDr|D) 
            for positions in the middle of a codon into the correct form for a 
            Pi array as if each position were the middle of a codon 
            corresponding to the 'amino acid'.
        self.min_delDl_given_DdelDr : list of lists
            minimum delDl for each delDr, D combination.
        self.max_delDl_given_DdelDr : list of lists
            maximum delDl for each delDr, D combination.
        self.PdelDldelDr_given_D : ndarray
            Joint probability distribution of the D deletions given the D allele,
            i.e. P(delDl, delDr |D)
        self.zeroD_given_D : list of floats
            The probability that a given D allele is fully deleted away.
        self.codons_dict : dict
            Dictionary, keyed by the allowed 'amino acid' symbols with the values 
            being lists of codons corresponding to the symbol.
        self.sub_codons_right : dict
            Dictionary of the 1 and 2 nucleotide suffixes (read from 5') for 
            each codon in an 'amino acid' grouping
    
        Returns
        -------
        Pi_L : ndarray
            (4, 3L) array corresponding to 
            \sum_D {D^{x_2}}_{x_3}{N^{x_3}}_{x_4}J(D)^{x_4}.
        
        """
        
        
        #Need to consider all D alignments from all possible positions and right deletions.
        nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        
        #n_aaseq = [aa_dict[aa] for aa in CDR3_seq]
        Pi_R = np.zeros((4, len(CDR3_seq)*3))
        min_pos = -len(CDR3_seq)*3
        
        num_dell_pos, num_delr_pos, num_D_genes = self.PdelDldelDr_given_D.shape
        
        for D_in, cutD_gen_seg  in enumerate(self.cutD_genomic_CDR3_segs):
            l_D_seg = len(cutD_gen_seg)
            
            #start position is first nt in a codon
            for init_pos in range(-1,-len(CDR3_seq)*3-1,-3):
                Pi_R[:, init_pos] += Pi_JinsDJ_given_D[D_in][:, init_pos]*self.zeroD_given_D[D_in]
                second_pos_dict = {'A': np.zeros(4), 'C': np.zeros(4), 'G': np.zeros(4), 'T': np.zeros(4)}
                codon_prefix_dict = {}
    
                for last_nt in 'ACGT':
                    for second_nt in 'ACGT':
                        codon_prefix_dict[last_nt + second_nt] = np.zeros(4)
                        #for first_nt in ['ACGT'[nt] for nt in range(4) if Pi_JinsDJ_given_D[D_in][nt, init_pos] > 0]:
                        for first_nt in 'ACGT':
                            if last_nt + second_nt + first_nt in self.codons_dict[CDR3_seq[init_pos//3]]: #possible allowed codon
                                second_pos_dict[second_nt][nt2num[last_nt]] += Pi_JinsDJ_given_D[D_in][nt2num[first_nt], init_pos] #base weight for middle pos nt
                                codon_prefix_dict[last_nt + second_nt][0] += Pi_JinsDJ_given_D[D_in][nt2num[first_nt], init_pos] #base weight for last pos nt
                for nt1 in 'ACGT':
                    if np.sum(second_pos_dict[nt1]) == 0:
                        second_pos_dict.pop(nt1, None)
                    for nt2 in 'ACGT':
                        if np.sum(codon_prefix_dict[nt1+nt2])== 0:
                            codon_prefix_dict.pop(nt1+nt2, None)
                            
    #            if len(second_pos_dict)> 0:
    #                print second_pos_dict
    #                return -1
                
                for delDr in range(num_delr_pos):
                    if self.min_delDl_given_DdelDr[D_in][delDr] == -1: # P(delDr | D) = 0 for this delDr --> move to next
                        continue
                    #Check if first nt from the D segment is okay
                    if cutD_gen_seg[l_D_seg - delDr - 1] in second_pos_dict.keys():
                        #The dell pos may be out of range of the PdelDldelDr_given_D -- check!
                        if l_D_seg - delDr - 1 <= self.max_delDl_given_DdelDr[D_in][delDr]:
                            Pi_R[:, init_pos - 1] += self.PdelDldelDr_given_D[l_D_seg - delDr - 1, delDr, D_in]*second_pos_dict[cutD_gen_seg[l_D_seg - delDr - 1]]
                    else:
                        continue #not okay, reject the alignment
                    
                    #Check if the second nt from the D segment is okay
                    if cutD_gen_seg[l_D_seg - delDr - 2:l_D_seg - delDr] in codon_prefix_dict.keys():
                        #The dell pos may be out of range of the PdelDldelDr_given_D -- check!
                        if l_D_seg - delDr - 2 <= self.max_delDl_given_DdelDr[D_in][delDr]:
                            Pi_R[0, init_pos - 2] += self.PdelDldelDr_given_D[l_D_seg - delDr - 2, delDr, D_in]*codon_prefix_dict[cutD_gen_seg[l_D_seg - delDr - 2:l_D_seg - delDr]][0]
                        base_prob = codon_prefix_dict[cutD_gen_seg[l_D_seg - delDr - 2:l_D_seg - delDr]][0]
                    else:
                        continue #no longer aligned, move to next delDr
                        
                    #Enter main loop
                    for pos in range(init_pos - 3, max(init_pos - l_D_seg + delDr, min_pos)-1, -1):
                        #note delDl = D_pos
                        D_pos = l_D_seg - delDr - 1 - ((init_pos - 1) - pos)
                        
                        #The dell pos may be out of range of the PdelDldelDr_given_D -- check!
                        if D_pos > self.max_delDl_given_DdelDr[D_in][delDr]:
                            current_PdelDldelDr = 0
                        else:
                            current_PdelDldelDr = self.PdelDldelDr_given_D[D_pos, delDr, D_in]
                        #Position is the first nt in codon
                        if pos%3 == 2:
                            #check alignment
                            if cutD_gen_seg[D_pos] in self.sub_codons_right[CDR3_seq[pos//3]]:
                                Pi_R[:, pos] += current_PdelDldelDr*base_prob*self.PD_nt_pos_vec[D_in][:, D_pos]
                            else:
                                break #no longer aligned -- exit loop
                        #Position is the second nt in codon
                        elif pos%3 == 1:
                            #check alignment
                            if cutD_gen_seg[D_pos:D_pos + 2] in self.sub_codons_right[CDR3_seq[pos//3]]:
                                Pi_R[:, pos] += current_PdelDldelDr*base_prob*self.PD_2nd_nt_pos_per_aa_vec[D_in][CDR3_seq[pos//3]][ :, D_pos]
                            else:
                                break #no longer aligned --- exit loop
                        #Position is the last nt in codon
                        else:
                            #check alignment
                            if cutD_gen_seg[D_pos:D_pos + 3] in self.codons_dict[CDR3_seq[pos//3]]:
                                Pi_R[0, pos] += current_PdelDldelDr*base_prob
                            else:
                                break #no longer aligned --- exit loop
                        
                    
                        
            #start position is second nt in a codon
            for init_pos in range(-2,-len(CDR3_seq)*3-1,-3):
                Pi_R[:, init_pos] += Pi_JinsDJ_given_D[D_in][:, init_pos]*self.zeroD_given_D[D_in]
                allowed_final_nts = ['ACGT'[nt] for nt in range(4) if Pi_JinsDJ_given_D[D_in][nt, init_pos] > 0]
                
                for delDr in range(num_delr_pos):
                    if self.min_delDl_given_DdelDr[D_in][delDr] == -1: # P(delDr | D) = 0 for this delDr --> move to next
                        continue
                    #check first nt of the D region (last in the codon)
                    if cutD_gen_seg[l_D_seg - delDr - 1] in allowed_final_nts: #first nt match
                        base_prob = Pi_JinsDJ_given_D[D_in][nt2num[cutD_gen_seg[l_D_seg - delDr - 1]], init_pos]
                        #The dell pos may be out of range of the PdelDldelDr_given_D -- check!
                        if l_D_seg - delDr - 1 <= self.max_delDl_given_DdelDr[D_in][delDr]:
                            Pi_R[0, init_pos-1] += self.PdelDldelDr_given_D[l_D_seg - delDr - 1, delDr, D_in]*base_prob
                    else:
                        continue #no alignment
                        
                    #Enter main loop
                    for pos in range(init_pos - 2, max(init_pos - l_D_seg + delDr, min_pos)-1, -1):
                        #note delDl = D_pos
                        D_pos = l_D_seg - delDr - 1 - ((init_pos - 1) - pos)
                        
                        #The dell pos may be out of range of the PdelDldelDr_given_D -- check!
                        if D_pos > self.max_delDl_given_DdelDr[D_in][delDr]:
                            current_PdelDldelDr = 0
                        else:
                            current_PdelDldelDr = self.PdelDldelDr_given_D[D_pos, delDr, D_in]
                        #Position is the first nt in codon
                        if pos%3 == 2:
                            #check alignment
                            if cutD_gen_seg[D_pos] in self.sub_codons_right[CDR3_seq[pos//3]]:
                                Pi_R[:, pos] += current_PdelDldelDr*base_prob*self.PD_nt_pos_vec[D_in][:, D_pos]
                            else:
                                break #no longer aligned -- exit loop
                        #Position is the second nt in codon
                        elif pos%3 == 1:
                            #check alignment
                            if cutD_gen_seg[D_pos:D_pos + 2] in self.sub_codons_right[CDR3_seq[pos//3]]:
                                Pi_R[:, pos] += current_PdelDldelDr*base_prob*self.PD_2nd_nt_pos_per_aa_vec[D_in][CDR3_seq[pos//3]][ :, D_pos]
                            else:
                                break #no longer aligned --- exit loop
                        #Position is the last nt in codon
                        else:
                            #check alignment
                            if cutD_gen_seg[D_pos:D_pos + 3] in self.codons_dict[CDR3_seq[pos//3]]:
                                Pi_R[0, pos] += current_PdelDldelDr*base_prob
                            else:
                                break #no longer aligned --- exit loop
                
            #start position is last nt in a codon
            for init_pos in range(-3,-len(CDR3_seq)*3-1,-3):
                Pi_R[0, init_pos] += Pi_JinsDJ_given_D[D_in][0, init_pos]*self.zeroD_given_D[D_in] 
                for delDr in range(num_delr_pos):
                    if self.min_delDl_given_DdelDr[D_in][delDr] == -1: # P(delDr | D) = 0 for this delDr --> move to next
                        continue
                    base_prob = Pi_JinsDJ_given_D[D_in][0, init_pos]
                    for pos in range(init_pos - 1, max(init_pos - l_D_seg + delDr, min_pos)-1, -1):
                        #note delDl = D_pos
                        D_pos = l_D_seg - delDr - 1 - ((init_pos - 1) - pos)
                        
                        #The dell pos may be out of range of the PdelDldelDr_given_D -- check!
                        if D_pos > self.max_delDl_given_DdelDr[D_in][delDr]:
                            current_PdelDldelDr = 0
                        else:
                            current_PdelDldelDr = self.PdelDldelDr_given_D[D_pos, delDr, D_in]
                        #Position is the first nt in codon
                        if pos%3 == 2:
                            #check alignment
                            if cutD_gen_seg[D_pos] in self.sub_codons_right[CDR3_seq[pos//3]]:
                                Pi_R[:, pos] += current_PdelDldelDr*base_prob*self.PD_nt_pos_vec[D_in][:, D_pos]
                            else:
                                break #no longer aligned -- exit loop
                        #Position is the second nt in codon
                        elif pos%3 == 1:
                            #check alignment
                            if cutD_gen_seg[D_pos:D_pos + 2] in self.sub_codons_right[CDR3_seq[pos//3]]:
                                Pi_R[:, pos] += current_PdelDldelDr*base_prob*self.PD_2nd_nt_pos_per_aa_vec[D_in][CDR3_seq[pos//3]][ :, D_pos]
                            else:
                                break #no longer aligned --- exit loop
                        #Position is the last nt in codon
                        else:
                            #check alignment
                            if cutD_gen_seg[D_pos:D_pos + 3] in self.codons_dict[CDR3_seq[pos//3]]:
                                Pi_R[0, pos] += current_PdelDldelDr*base_prob
                            else:
                                break #no longer aligned --- exit loop
                                
        return Pi_R



#%%    
class GenerationProbabilityVJ(GenerationProbability, PreprocessedParametersVJ):
    """Class used to compute the Pgen of CDR3 sequences from a  VJ model.
    
    All of the attributes of GenerationProbabilityVJ are inherited from the
    class PreprocessedParametersVJ. The methods of the class are used to 
    compute the Pgen of an 'amino acid' sequence from a VJ generative model.

    Attributes
    ----------
    codons_dict : dict
        Dictionary, keyed by the allowed 'amino acid' symbols with the values 
        being lists of codons corresponding to the symbol.        
    sub_codons_left : dict
        Dictionary of the 1 and 2 nucleotide prefixes (read from 5') for 
        each codon in an 'amino acid' grouping        
    sub_codons_right : dict
        Dictionary of the 1 and 2 nucleotide suffixes (read from 5') for 
        each codon in an 'amino acid' grouping
        
    V_allele_names : list of strings
        List of V allele names in genomic_data        
    d_V_usage_mask : list of int
        Default V usage mask of indices of all productive V genes/alleles.        
    V_mask_mapping : dict
        Dictionary mapping allowed keywords (strings corresponding to V gene 
        and V allele names) to the indices of the V alleles they refer to.        
    cutV_genomic_CDR3_segs : list of strings
        List of the V germline nucleotide sequences, trimmed to begin at the
        CDR3 region (includes the conserved C residue) with the maximum number
        of reverse complementary palindromic insertions appended. 
    PVJ : ndarray
        Joint probability distribution of V and J, P(V, J).
    PVdelV_nt_pos_vec : list of ndarrays
        For each V allele, format P(delV|V) into the correct form for a Pi 
        array or V(J)_{x_1}. This is only done for the first and last position 
        in each codon.    
    PVdelV_2nd_nt_pos_per_aa_vec : list of dicts
        For each V allele, and each 'amino acid', format P(V)*P(delV|V) for 
        positions in the middle of a codon into the correct form for a Pi array 
        or V(J)_{x_1} given the 'amino acid'.
        
    J_allele_names : list of strings
        List of J allele names in genomic_data        
    d_J_usage_mask : list of int
        Default J usage mask of indices of all productive J genes/alleles.        
    J_mask_mapping : dict
        Dictionary mapping allowed keywords (strings corresponding to V gene 
        and V allele names) to the indices of the V alleles they refer to.        
    cutJ_genomic_CDR3_segs : list of strings
        List of the J germline nucleotide sequences, trimmed to end at the
        CDR3 region (includes the conserved F or W residue) with the maximum 
        number of reverse complementary palindromic insertions appended.    
    PJdelJ_nt_pos_vec : list of ndarrays
        For each J allele, format P(delJ|J) into the correct form for a Pi 
        array or J^{x_2}. This is only done for the first and last position in 
        each codon.           
    PJdelJ_2nd_nt_pos_per_aa_vec : list of dicts
        For each J allele, and each 'amino acid', format P(delJ|J) for 
        positions in the middle of a codon into the correct form for a Pi array 
        or J^{x_2} given the 'amino acid'.
    
    PinsVJ : ndarray
        Probability distribution of the VJ (N) insertion sequence length        
    Rvj : ndarray
        Markov transition matrix for the VJ insertion junction.        
    first_nt_bias_insVJ : ndarray
        (4,) array of the probability distribution of the indentity of the 
        first nucleotide insertion for the VJ junction.        
    zero_nt_bias_insVJ : ndarray
        (4,) array of the probability distribution of the indentity of the 
        the nucleotide BEFORE the VJ insertion.
        zero_nt_bias_insVJ = Rvj^{-1}first_nt_bias_insVJ        
    Tvj, Svj, Dvj, lTvj, lDvj : dicts
        Dictionaries (by 'amino acid') of insertion transfer matrices 
        ((4, 4) ndarrays) for the VJ junction.
    
    """
    def __init__(self, generative_model, genomic_data, alphabet_file = None):
        """Initialize GenerationProbabilityVJ.
        
        This intialization inherits all of the attributes of 
        PreprocessedParametersVJ (which include all of the processed 
        parameters needed for Pgen computation) and the methods of 
        GenerationProbability which include some wrappers/formatting of
        sequences to make Pgen computation of nucleotide and regular expression
        sequences easier (etc).
        
        Parameters
        ----------
        generative_model : GenerativeModelVJ
            VJ generative model class containing the model parameters.            
        genomic_data : GenomicDataVJ
            VJ genomic data class containing the V and J germline sequences and 
            info.            
        alphabet_file : str, optional
            File name (full pathing from current directory) for a custom alphabet
            definition. If no file is provided, the default alphabet is used, i.e. 
            standard amino acids, undetermined amino acids (B, J, X, and Z), and
            single codon symbols.
        
        """
        GenerationProbability.__init__(self)
        PreprocessedParametersVJ.__init__(self, generative_model, genomic_data, alphabet_file)
        
    def compute_CDR3_pgen(self, CDR3_seq, V_usage_mask, J_usage_mask):
        """Compute Pgen for CDR3 'amino acid' sequence CDR3_seq from VJ model.
    
        Conditioned on the already formatted V genes/alleles indicated in 
        V_usage_mask and the J genes/alleles in J_usage_mask.
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' (single character symbols
            each corresponding to a collection of codons as given by codons_dict).
        V_usage_mask : list
            Indices of the V alleles to be considered in the Pgen computation
        J_usage_mask : list
            Indices of the J alleles to be considered in the Pgen computation
    
        Returns
        -------
        pgen : float
            The generation probability (Pgen) of the sequence
        
        Examples
        --------
        >>> compute_CDR3_pgen('CAVKIQGAQKLVF', ppp, [72], [56])
        4.1818202431143785e-07
        >>> compute_CDR3_pgen(nt2codon_rep('TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC'), ppp, [42], [1])
        1.3971676613008565e-08
        >>> compute_CDR3_pgen('\xbb\xb6\xbe\x80\xbc\xa1\x8a\x96\xa1\xa0\xad\x8e\xbf', ppp, [72], [56])
        1.3971676613008565e-08
        
        """
        
        #Genomic J alignment/matching (contribution from P(delJ | J)), return Pi_J and reduced J_usage_mask
        Pi_J, r_J_usage_mask = self.compute_Pi_J(CDR3_seq, J_usage_mask)
        
        #Genomic V alignment/matching conditioned on J gene (contribution from P(V, J, delV)), return Pi_V_given_J
        Pi_V_given_J, max_V_align = self.compute_Pi_V_given_J(CDR3_seq, V_usage_mask, r_J_usage_mask)
        
        #Include insertions (R and PinsVJ) to get the total contribution from the left (3') side conditioned on J gene. Return Pi_V_insVJ_given_J
        Pi_V_insVJ_given_J = self.compute_Pi_V_insVJ_given_J(CDR3_seq, Pi_V_given_J, max_V_align)
        
        pgen = 0
        #zip Pi_V_insVJ_given_J and Pi_J together for each J gene to get total pgen
        for j in range(len(r_J_usage_mask)):
            for pos in range(len(CDR3_seq)*3 - 1):
                pgen += np.dot(Pi_V_insVJ_given_J[j][:, pos], Pi_J[j][:, pos+1])
        return pgen
    
    #Genomic V alignment/matching (contribution from P(V, delV)), return Pi_V
    def compute_Pi_V_given_J(self, CDR3_seq, V_usage_mask, J_usage_mask):
        """Compute Pi_V conditioned on J.
    
        This function returns the Pi array from the model factors of the V genomic 
        contributions, P(V, J)*P(delV|V). This corresponds to V(J)_{x_1}.
        
        For clarity in parsing the algorithm implementation, we include which 
        instance attributes are used in the method as 'parameters.'
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' (single character symbols
            each corresponding to a collection of codons as given by codons_dict).
        V_usage_mask : list
            Indices of the V alleles to be considered in the Pgen computation
        J_usage_mask : list
            Indices of the J alleles to be considered in the Pgen computation
            
        self.cutV_genomic_CDR3_segs : list of strings
            List of all the V genomic nucleotide sequences trimmed to begin at the 
            conserved C residue and with the maximum number of palindromic 
            insertions appended.
        self.PVdelV_nt_pos_vec : list of ndarrays
            For each V allele, format P(delV|V) into the correct form for a Pi 
            array or V(J)_{x_1}. This is only done for the first and last 
            position in each codon.    
        self.PVdelV_2nd_nt_pos_per_aa_vec : list of dicts
            For each V allele, and each 'amino acid', format P(V)*P(delV|V) for 
            positions in the middle of a codon into the correct form for a Pi 
            array or V(J)_{x_1} given the 'amino acid'.
        self.PVJ : ndarray
            Joint probability distribution of V and J, P(V, J).
    
        Returns
        -------
        Pi_V_given_J : list
            List of (4, 3L) ndarrays corresponding to V(J)_{x_1}.
        max_V_align: int
            Maximum alignment of the CDR3_seq to any genomic V allele allowed by
            V_usage_mask.
            
        """
        
        #Note, the cutV_genomic_CDR3_segs INCLUDE the palindromic insertions and thus are max_palindrome nts longer than the template.
        #furthermore, the genomic sequence should be pruned to start at the conserved C
        
        Pi_V_given_J = [np.zeros((4, len(CDR3_seq)*3)) for i in J_usage_mask] #Holds the aggregate weight for each nt possiblity and position
        alignment_lengths = []
        for V_in in V_usage_mask:
            try:
                cutV_gen_seg = self.cutV_genomic_CDR3_segs[V_in]
            except IndexError:
                print('Check provided V usage mask. Contains indicies out of allowed range.')
                continue
            current_alignment_length = self.max_nt_to_aa_alignment_left(CDR3_seq, cutV_gen_seg)
            alignment_lengths += [current_alignment_length]
            current_Pi_V = np.zeros((4, len(CDR3_seq)*3))
            
            if current_alignment_length > 0:
                #For first and last nt in a codon use PVdelV_nt_pos_vec
                current_Pi_V[:, :current_alignment_length] = self.PVdelV_nt_pos_vec[V_in][:, :current_alignment_length]
                for pos in range(1, current_alignment_length, 3): #for middle nt use PVdelV_2nd_nt_pos_per_aa_vec
                    current_Pi_V[:, pos] = self.PVdelV_2nd_nt_pos_per_aa_vec[V_in][CDR3_seq[pos//3]][:, pos]
                for j, J_in in enumerate(J_usage_mask):
                    Pi_V_given_J[j][:, :current_alignment_length] += self.PVJ[V_in, J_in]*current_Pi_V[:, :current_alignment_length]
        
        return Pi_V_given_J, max(alignment_lengths)
    
    #Include insertions (R and PinsVJ) to get the total contribution from the the V and insertions conditioned on J identity. Return Pi_V_insVJ_given_J
    def compute_Pi_V_insVJ_given_J(self, CDR3_seq, Pi_V_given_J, max_V_align):
        """Compute Pi_V_insVJ conditioned on J.
    
        This function returns the Pi array from the model factors of the V genomic 
        contributions, P(V, J)*P(delV|V), and the VJ (N) insertions,
        first_nt_bias_insVJ(m_1)PinsVJ(\ell_{VJ})\prod_{i=2}^{\ell_{VJ}}Rvj(m_i|m_{i-1}). 
        This corresponds to V(J)_{x_1}{M^{x_1}}_{x_2}.
        
        For clarity in parsing the algorithm implementation, we include which 
        instance attributes are used in the method as 'parameters.'
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' (single character symbols
            each corresponding to a collection of codons as given by codons_dict).
        Pi_V_given_J : ndarray
            List of (4, 3L) ndarrays corresponding to V(J)_{x_1}.
        max_V_align : int
            Maximum alignment of the CDR3_seq to any genomic V allele allowed by
            V_usage_mask.
            
        self.PinsVJ : ndarray
            Probability distribution of the VJ insertion sequence length
        self.first_nt_bias_insVJ : ndarray
            (4,) array of the probability distribution of the indentity of the 
            first nucleotide insertion for the VJ junction.        
        self.zero_nt_bias_insVJ : ndarray
            (4,) array of the probability distribution of the indentity of the 
            the nucleotide BEFORE the VJ insertion.
            zero_nt_bias_insVJ = Rvj^{-1}first_nt_bias_insVJ 
        self.Tvj : dict
            Dictionary of full codon transfer matrices ((4, 4) ndarrays) by 
            'amino acid'.
        self.Svj : dict
            Dictionary of transfer matrices ((4, 4) ndarrays) by 'amino acid' for 
            the VD insertion ending in the first position.
        self.Dvj : dict
            Dictionary of transfer matrices ((4, 4) ndarrays) by 'amino acid' for 
            the VD insertion ending in the second position.
        self.lTvj : dict
            Dictionary of transfer matrices ((4, 4) ndarrays) by 'amino acid' for 
            the VD insertion starting in the first position.
        self.lDvj : dict
            Dictionary of transfer matrices ((4, 4) ndarrays) by 'amino acid' for
            VD insertion starting in the first position and ending in the second 
            position of the same codon.
    
        Returns
        -------
        Pi_V_insVJ_given_J : list
            List of (4, 3L) ndarrays corresponding to V(J)_{x_1}{M^{x_1}}_{x_2}.
            
        """
        #max_insertions = 30 #len(PinsVJ) - 1 should zeropad the last few spots
        max_insertions = len(self.PinsVJ) - 1
        
        Pi_V_insVJ_given_J = [np.zeros((4, len(CDR3_seq)*3)) for i in range(len(Pi_V_given_J))]
        
        
        for j in range(len(Pi_V_given_J)):
            #start position is first nt in a codon
            for init_pos in range(0, max_V_align, 3):
                #Zero insertions
                Pi_V_insVJ_given_J[j][:, init_pos] += self.PinsVJ[0]*Pi_V_given_J[j][:, init_pos]
                
                #One insertion
                Pi_V_insVJ_given_J[j][:, init_pos+1] += self.PinsVJ[1]*np.dot(self.lDvj[CDR3_seq[init_pos//3]], Pi_V_given_J[j][:, init_pos])
        
                #Two insertions and compute the base nt vec for the standard loop        
                current_base_nt_vec = np.dot(self.lTvj[CDR3_seq[init_pos//3]], Pi_V_given_J[j][:, init_pos])
                Pi_V_insVJ_given_J[j][0, init_pos+2] += self.PinsVJ[2]*np.sum(current_base_nt_vec)
                
                base_ins = 2
                
                #Loop over all other insertions using base_nt_vec
                for aa in CDR3_seq[init_pos//3 + 1: init_pos//3 + max_insertions//3]:
                    Pi_V_insVJ_given_J[j][:, init_pos+base_ins+1] += self.PinsVJ[base_ins + 1]*np.dot(self.Svj[aa], current_base_nt_vec)
                    Pi_V_insVJ_given_J[j][:, init_pos+base_ins+2] += self.PinsVJ[base_ins + 2]*np.dot(self.Dvj[aa], current_base_nt_vec)
                    current_base_nt_vec = np.dot(self.Tvj[aa], current_base_nt_vec)
                    Pi_V_insVJ_given_J[j][0, init_pos+base_ins+3] += self.PinsVJ[base_ins + 3]*np.sum(current_base_nt_vec)
                    base_ins +=3
                
            
            #start position is second nt in a codon
            for init_pos in range(1, max_V_align, 3):
                #Zero insertions
                Pi_V_insVJ_given_J[j][:, init_pos] += self.PinsVJ[0]*Pi_V_given_J[j][:, init_pos]
                #One insertion --- we first compute our p vec by pairwise mult with the ss distr
                current_base_nt_vec = np.multiply(Pi_V_given_J[j][:, init_pos], self.first_nt_bias_insVJ)
                Pi_V_insVJ_given_J[j][0, init_pos+1] += self.PinsVJ[1]*np.sum(current_base_nt_vec)
                base_ins = 1
                
                #Loop over all other insertions using base_nt_vec
                for aa in CDR3_seq[init_pos//3 + 1: init_pos//3 + max_insertions//3]:
                    Pi_V_insVJ_given_J[j][:, init_pos+base_ins+1] += self.PinsVJ[base_ins + 1]*np.dot(self.Svj[aa], current_base_nt_vec)
                    Pi_V_insVJ_given_J[j][:, init_pos+base_ins+2] += self.PinsVJ[base_ins + 2]*np.dot(self.Dvj[aa], current_base_nt_vec)
                    current_base_nt_vec = np.dot(self.Tvj[aa], current_base_nt_vec)
                    Pi_V_insVJ_given_J[j][0, init_pos+base_ins+3] += self.PinsVJ[base_ins + 3]*np.sum(current_base_nt_vec)
                    base_ins +=3
                
            #start position is last nt in a codon   
            for init_pos in range(2, max_V_align, 3):
                
                #Zero insertions
                Pi_V_insVJ_given_J[j][0, init_pos] += self.PinsVJ[0]*Pi_V_given_J[j][0, init_pos]
                #current_base_nt_vec = first_nt_bias_insVJ*Pi_V_given_J[j][0, init_pos] #Okay for steady state
                current_base_nt_vec = self.zero_nt_bias_insVJ*Pi_V_given_J[j][0, init_pos]
                base_ins = 0
                
                #Loop over all other insertions using base_nt_vec
                for aa in CDR3_seq[init_pos//3 + 1: init_pos//3 + max_insertions//3]:
                    Pi_V_insVJ_given_J[j][:, init_pos+base_ins+1] += self.PinsVJ[base_ins + 1]*np.dot(self.Svj[aa], current_base_nt_vec)
                    Pi_V_insVJ_given_J[j][:, init_pos+base_ins+2] += self.PinsVJ[base_ins + 2]*np.dot(self.Dvj[aa], current_base_nt_vec)
                    current_base_nt_vec = np.dot(self.Tvj[aa], current_base_nt_vec)
                    Pi_V_insVJ_given_J[j][0, init_pos+base_ins+3] += self.PinsVJ[base_ins + 3]*np.sum(current_base_nt_vec)
                    base_ins +=3
    
         
        return Pi_V_insVJ_given_J
    
    #Genomic J alignment/matching (contribution from P(delJ | J)), return Pi_J and the reduced J_usage_mask (reduced based on non-zero alignment)
    def compute_Pi_J(self, CDR3_seq, J_usage_mask):
        """Compute Pi_J.
    
        This function returns the Pi array from the model factors of the J genomic 
        contributions, P(delJ|J). This corresponds to J(D)^{x_4}.
        
        For clarity in parsing the algorithm implementation, we include which 
        instance attributes are used in the method as 'parameters.'
    
        Parameters
        ----------
        CDR3_seq : str
            CDR3 sequence composed of 'amino acids' (single character symbols
            each corresponding to a collection of codons as given by codons_dict).
        J_usage_mask : list
            Indices of the J alleles to be considered in the Pgen computation
        
        self.cutJ_genomic_CDR3_segs : list of strings
            List of all the J genomic nucleotide sequences trimmed to begin at the 
            conserved 3' residue (F/W) and with the maximum number of palindromic 
            insertions appended.
        self.PJdelJ_nt_pos_vec : list of ndarrays
            For each J allele, format P(delJ|J) into the correct form for a Pi 
            array or J^{x_2}. This is only done for the first and last position 
            in each codon.    
        self.PJdelJ_2nd_nt_pos_per_aa_vec : list of dicts
            For each J allele, and each 'amino acid', format P(delJ|J) for 
            positions in the middle of a codon into the correct form for a Pi 
            array or J^{x_2} given the 'amino acid'.
        
        Returns
        -------
        Pi_J : ndarray
            (4, 3L) array corresponding to J^{x_4}.
        r_J_usage_mask: list
            Reduced J_usage mask. J genes/alleles with no contribution (bad 
            alignment) are removed from the mask. This is done to speed up the 
            computation on the V side (which must be done conditioned on the J).
        
        """
        #Note, the cutJ_genomic_CDR3_segs INCLUDE the palindromic insertions and thus are max_palindrome nts longer than the template.
        #furthermore, the genomic sequence should be pruned to start at a conserved region on the J side
        
        Pi_J = [] #Holds the aggregate weight for each nt possiblity and position
        r_J_usage_mask = []
        for j, J_in in enumerate(J_usage_mask):
            try:
                cutJ_gen_seg = self.cutJ_genomic_CDR3_segs[J_in]
            except IndexError:
                print('Check provided J usage mask. Contains indicies out of allowed range.')
                continue
            current_alignment_length = self.max_nt_to_aa_alignment_right(CDR3_seq, cutJ_gen_seg)
            #alignment_lengths += [current_alignment_length]
            current_Pi_J = np.zeros((4, len(CDR3_seq)*3))
            
            if current_alignment_length > 0:
                #For first and last nt in a codon use PJdelJ_nt_pos_vec
                current_Pi_J[:, -current_alignment_length:] = self.PJdelJ_nt_pos_vec[J_in][:, -current_alignment_length:]          
                for pos in range(-2, -current_alignment_length-1, -3): #for middle nt use PJdelJ_2nd_nt_pos_per_aa_vec
                    current_Pi_J[:, pos] = self.PJdelJ_2nd_nt_pos_per_aa_vec[J_in][CDR3_seq[pos//3]][:, pos]
                if np.sum(current_Pi_J) > 0:
                    Pi_J.append(current_Pi_J)
                    r_J_usage_mask.append(J_in)
        
        return Pi_J, r_J_usage_mask