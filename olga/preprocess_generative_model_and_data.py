#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process classes GenerativeModelV(D)J GenomicModelV(D)J before computing Pgen.

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
    

As many of the arrays that will be used to compute the Pgen of a specific 
sequence depend only on the generative model and the germline sequences (as 
opposed to the CDR3 sequence itself) much of the computation and organization
can be done in advance and does not need to be repeated for each individual 
sequence. This module prepares and formats all of the sequence-independent 
components which are stored as attributes of classes PreprocessedParametersVDJ 
and PreprocessedParametersVJ which will be inherited by the 
GenerationProbabilityV(D)J classes.

There is also an option of defining a custom 'amino acid' alphabet. Please 
refer to the example alphabet file and the README for more info on how to 
define customized alphabets for faster Pgen computation of CDR3 motifs.

This module and the classes PreprocessedParametersV(D)J are unlikely to be
used directly, but are used in the instantiation of the classes
GenerationProbabilityV(D)J.

@author: zacharysethna

"""
from __future__ import division
import numpy as np
from olga.utils import construct_codons_dict, generate_sub_codons_left, generate_sub_codons_right, calc_steady_state_dist, gene_to_num_str

class PreprocessedParameters(object):
    """Class used to preprocess the parameters that both VDJ and VJ models have.

    The class is the parent of classes PreprocessedParametersVDJ and 
    PreprocessedParametersVJ.

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
        
    """
    def __init__(self, generative_model, genomic_data, alphabet_file = None):
        """Initialize PreprocessedParameters
        
        This initialization is used to compute the preprocess parameters that 
        are in common with both the VDJ and VJ model types.
        
        Parameters
        ----------
        generative_model : GenerativeModelVDJ or GenerativeModelVJ
            V(D)J generative model classes containing the model parameters.            
        genomic_data : GenomicDataVDJ GenomicDataVJ
            V(D)J genomic data classes containing the V(,D,) and J germline 
            sequences and info.            
        alphabet_file : str, optional
            File name (full pathing from current directory) for a custom alphabet
            definition. If no file is provided, the default alphabet is used, i.e. 
            standard amino acids, undetermined amino acids (B, J, X, and Z), and
            single codon symbols.
        
        """
        
        #determine if generative_model and genomic_data are VDJ or VJ
        if generative_model.__class__.__name__.endswith('VDJ') and genomic_data.__class__.__name__.endswith('VDJ'):
            recomb_type = 'VDJ'
        elif generative_model.__class__.__name__.endswith('VJ') and genomic_data.__class__.__name__.endswith('VJ'):
            recomb_type = 'VJ'
        else:
            raise ValueError #recomb types must match
        
        
        self.codons_dict = construct_codons_dict(alphabet_file)
        self.sub_codons_left = generate_sub_codons_left(self.codons_dict)
        self.sub_codons_right = generate_sub_codons_right(self.codons_dict)
        
        self.V_allele_names = None
        self.V_mask_mapping = None
        self.J_allele_names = None
        self.J_mask_mapping = None
        self.make_V_and_J_mask_mapping(genomic_data.genV, genomic_data.genJ)
        
        self.cutV_genomic_CDR3_segs = genomic_data.cutV_genomic_CDR3_segs
        
        #construct the default V_usage_mask which represents all functional V genes with non-zero probability
        d_V_usage_mask = []
        for V_in, cutVseg in enumerate(genomic_data.cutV_genomic_CDR3_segs):
            if recomb_type == 'VDJ':
                if len(cutVseg) > 0 and generative_model.PV[V_in] > 0:
                    d_V_usage_mask.append(V_in)
            elif recomb_type == 'VJ':
                if len(cutVseg) > 0 and np.sum(generative_model.PVJ[V_in, :]) > 0:
                    d_V_usage_mask.append(V_in)
        self.d_V_usage_mask = d_V_usage_mask
        
        
        
        self.cutJ_genomic_CDR3_segs = genomic_data.cutJ_genomic_CDR3_segs
        
        #construct the default J_usage_mask which represents all functional J genes with non-zero probability
        d_J_usage_mask = []
        for J_in, cutJseg in enumerate(genomic_data.cutJ_genomic_CDR3_segs):
            if recomb_type == 'VDJ':
                if len(cutJseg) > 0 and np.sum(generative_model.PDJ[:, J_in]) > 0:
                    d_J_usage_mask.append(J_in)
            elif recomb_type == 'VJ':
                if len(cutJseg) > 0 and np.sum(generative_model.PVJ[:, J_in]) > 0:
                    d_J_usage_mask.append(J_in)
        self.d_J_usage_mask = d_J_usage_mask
        
        
    def make_V_and_J_mask_mapping(self, genV, genJ):
        """Constructs the V and J mask mapping dictionaries.
        
        Parameters
        ----------
        genV : list
            List of genomic V information.            
        genJ : list
            List of genomic J information.
        
        """
        #construct mapping between allele/gene names and index for custom V_usage_masks
        V_allele_names = [V[0] for V in genV]
        V_reduced_allele_names = [gene_to_num_str(v, 'v') for v in V_allele_names]
        V_mask_mapping = {}
        for v in [x.split('*')[0] for x in V_reduced_allele_names]:
            V_mask_mapping[v] = []
            if '-' not in v:
                V_mask_mapping[v+'-1'] = []
            else:
                V_mask_mapping[v.split('-')[0]] = []

        for i, v in enumerate(V_reduced_allele_names):
            V_mask_mapping[v] = [i]
            V_mask_mapping[v.split('*')[0]].append(i)
            if '-' not in v:
                V_mask_mapping[v.replace('*', '-1*')] = [i]
                V_mask_mapping[v.split('*')[0] + '-1'].append(i)
            else:
                V_mask_mapping[v.split('*')[0].split('-')[0]].append(i)

        #construct mapping between allele/gene names and index for custom J_usage_masks
        J_allele_names = [J[0] for J in genJ]
        J_reduced_allele_names = [gene_to_num_str(j, 'j') for j in J_allele_names]
        J_mask_mapping = {}
        for j in [x.split('*')[0] for x in J_reduced_allele_names]:
            J_mask_mapping[j] = []
            if '-' not in j:
                J_mask_mapping[j+'-1'] = []
            else:
                J_mask_mapping[j.split('-')[0]] = []
        for i, j in enumerate(J_reduced_allele_names):
            J_mask_mapping[j] = [i]
            J_mask_mapping[j.split('*')[0]].append(i)
            if '-' not in j:
                J_mask_mapping[j.replace('*', '-1*')] = [i]
                J_mask_mapping[j.split('*')[0] + '-1'].append(i)
            else:
                J_mask_mapping[j.split('*')[0].split('-')[0]].append(i)
            
        self.V_allele_names = V_allele_names
        self.V_mask_mapping = V_mask_mapping
        self.J_allele_names = J_allele_names
        self.J_mask_mapping = J_mask_mapping

#%% PreprocessedParametersVJ Class def
class PreprocessedParametersVDJ(PreprocessedParameters):
    """Class used to preprocess the parameters of VDJ recombination models.

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
        """Initialize PreprocessedParametersVDJ
        
        This intialization computes all of the attributes that will be needed
        for Pgen computation of CDR3 sequences generated by VDJ recombination.
        
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
        PreprocessedParameters.__init__(self, generative_model, genomic_data, alphabet_file)
        
        #Check that the generative_model and genomic_data are VDJ
        if not all([generative_model.__class__.__name__.endswith('VDJ'),genomic_data.__class__.__name__.endswith('VDJ')]):
            raise ValueError #Need VDJ model and data
        
        self.cutD_genomic_CDR3_segs = genomic_data.cutD_genomic_CDR3_segs
        self.D_allele_names = [D[0] for D in genomic_data.genD]
        
        #Set the V genomic Pi arrays
        self.PVdelV_nt_pos_vec = None
        self.PVdelV_2nd_nt_pos_per_aa_vec = None
        self.generate_PVdelV_nt_pos_vecs(generative_model, genomic_data) #Computes and sets the above attributes
        
        #Set the D genomic Pi arrays and info
        self.PD_nt_pos_vec = None
        self.PD_2nd_nt_pos_per_aa_vec = None
        self.min_delDl_given_DdelDr = None
        self.max_delDl_given_DdelDr = None
        self.zeroD_given_D = None
        self.preprocess_D_segs(generative_model, genomic_data) #Computes and sets the above attributes
        
        self.PdelDldelDr_given_D = generative_model.PdelDldelDr_given_D
        
        #Set the J genomic Pi arrays and info
        self.PJdelJ_nt_pos_vec = None
        self.PJdelJ_2nd_nt_pos_per_aa_vec = None
        self.generate_PJdelJ_nt_pos_vecs(generative_model, genomic_data) #Computes and sets the above attributes
        
        #allow for 0 prob J genes
        self.PD_given_J = np.zeros(generative_model.PDJ.shape)
        self.PD_given_J[:, np.sum(generative_model.PDJ, axis = 0) > 0] = np.multiply(generative_model.PDJ[:, np.sum(generative_model.PDJ, axis = 0) > 0], 1/np.sum(generative_model.PDJ, axis = 0)[np.sum(generative_model.PDJ, axis = 0) > 0])

        #Trim, then zeropad PinsVD and PinsDJ
        self.PinsVD = np.append(np.trim_zeros(generative_model.PinsVD, 'b'),  [0., 0., 0., 0.])
        self.PinsDJ = np.append(np.trim_zeros(generative_model.PinsDJ, 'b'),  [0., 0., 0., 0.])
        
        
        #insVD and insDJ Dinucleotide bias transition matrices
        #Note, these used to be called 'RnucleotideVD_per_nucleotideVD_5prime'
        #and 'RnucleotideDJ_per_nucleotideDJ_3prime'
        self.Rvd = generative_model.Rvd
        self.Rdj = generative_model.Rdj
        
        
        #Check if first insertion nt probability distributions are given. 
        #Note, this is not normally inferred by IGoR. If these nt bias dists
        #aren't present, use the steady-state distributions defined from Rvd
        #and Rdj. When not using the steady-state distributions we need to
        #compute the nt bias for the position previous to first_nt_bias 
        #(R^{-1}p_0)
        
        if generative_model.first_nt_bias_insVD == None:
            self.first_nt_bias_insVD = calc_steady_state_dist(self.Rvd)
            #In steady-state zero_nt_bias_insVD = first_nt_bias_insVD so no need to recalculate
            self.zero_nt_bias_insVD = self.first_nt_bias_insVD
        else:
            self.first_nt_bias_insVD = generative_model.first_nt_bias_insVD
            #Require Rvd*zero_nt_bias_insVD = first_nt_bias_insVD  --- can use pseudo-inverse
            self.zero_nt_bias_insVD = np.dot(np.linalg.pinv(self.Rvd), generative_model.first_nt_bias_insVD)
        
        if generative_model.first_nt_bias_insDJ == None:
            self.first_nt_bias_insDJ = calc_steady_state_dist(self.Rdj)
            #In steady-state zero_nt_bias_insDJ = first_nt_bias_insDJ so no need to recalculate
            self.zero_nt_bias_insDJ = self.first_nt_bias_insDJ
        else:
            self.first_nt_bias_insDJ = generative_model.first_nt_bias_insDJ
            #Require Rdj*zero_nt_bias_insDJ = first_nt_bias_insDJ  --- can use pseudo-inverse
            self.zero_nt_bias_insDJ = np.dot(np.linalg.pinv(self.Rdj), generative_model.first_nt_bias_insDJ)
            
        
        #Compute VD insertion transfer matrices
        self.Tvd = None
        self.Svd = None
        self.Dvd = None
        self.lTvd = None
        self.lDvd = None
        self.generate_VD_junction_transfer_matrices() #Computes and sets the transfer matrices
        
        #Compute DJ insertion transfer matrices
        self.Tdj = None
        self.Sdj = None
        self.Ddj = None
        self.rTdj = None
        self.rDdj = None
        self.generate_DJ_junction_transfer_matrices() #Computes and sets the transfer matrices

    def generate_PVdelV_nt_pos_vecs(self, generative_model, genomic_data):
        """Process P(V)*P(delV|V) into Pi arrays.
        
        Sets the attributes PVdelV_nt_pos_vec and PVdelV_2nd_nt_pos_per_aa_vec.
    
        Parameters
        ----------
        generative_model : GenerativeModelVDJ
            VDJ generative model class containing the model parameters.            
        genomic_data : GenomicDataVDJ
            VDJ genomic data class containing the V, D, and J germline 
            sequences and info.
        
        """
    
        cutV_genomic_CDR3_segs = genomic_data.cutV_genomic_CDR3_segs
        nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        num_del_pos = generative_model.PdelV_given_V.shape[0]
        
        PVdelV_nt_pos_vec = [[]]*len(generative_model.PV)
        PVdelV_2nd_nt_pos_per_aa_vec = [[]]*len(generative_model.PV)
        for V_in, pv in enumerate(generative_model.PV):
            current_PVdelV_nt_pos_vec = np.zeros((4, len(cutV_genomic_CDR3_segs[V_in])))
            current_PVdelV_2nd_nt_pos_per_aa_vec = {}
            for aa in self.codons_dict.keys():
                current_PVdelV_2nd_nt_pos_per_aa_vec[aa] = np.zeros((4, len(cutV_genomic_CDR3_segs[V_in])))
            for pos, nt in enumerate(cutV_genomic_CDR3_segs[V_in]):
                if len(cutV_genomic_CDR3_segs[V_in]) - pos >  num_del_pos:
                    continue
                if pos%3 == 0: #Start of a codon
                    current_PVdelV_nt_pos_vec[nt2num[nt], pos] = pv*generative_model.PdelV_given_V[len(cutV_genomic_CDR3_segs[V_in])-pos-1, V_in]    
                elif pos%3 == 1: #Mid codon position
                    for ins_nt in 'ACGT':
                        #We need to find what possible codons are allowed for any aa (or motif)
                        for aa in self.codons_dict.keys():
                            if cutV_genomic_CDR3_segs[V_in][pos-1:pos+1]+ ins_nt in self.codons_dict[aa]:
                                current_PVdelV_2nd_nt_pos_per_aa_vec[aa][nt2num[ins_nt], pos] = pv*generative_model.PdelV_given_V[len(cutV_genomic_CDR3_segs[V_in])-pos-1, V_in]            
                elif pos%3 == 2: #End of codon
                    current_PVdelV_nt_pos_vec[0, pos] = pv*generative_model.PdelV_given_V[len(cutV_genomic_CDR3_segs[V_in])-pos-1, V_in]
            PVdelV_nt_pos_vec[V_in] = current_PVdelV_nt_pos_vec
            PVdelV_2nd_nt_pos_per_aa_vec[V_in] = current_PVdelV_2nd_nt_pos_per_aa_vec
    
        self.PVdelV_nt_pos_vec = PVdelV_nt_pos_vec
        self.PVdelV_2nd_nt_pos_per_aa_vec = PVdelV_2nd_nt_pos_per_aa_vec
        
    def preprocess_D_segs(self, generative_model, genomic_data):
        """Process P(delDl, delDr|D) into Pi arrays.
        
        Sets the attributes PD_nt_pos_vec, PD_2nd_nt_pos_per_aa_vec, 
        min_delDl_given_DdelDr, max_delDl_given_DdelDr, and zeroD_given_D.
    
        Parameters
        ----------
        generative_model : GenerativeModelVDJ
            VDJ generative model class containing the model parameters.            
        genomic_data : GenomicDataVDJ
            VDJ genomic data class containing the V, D, and J germline 
            sequences and info.
        
        """
    
        cutD_genomic_CDR3_segs = genomic_data.cutD_genomic_CDR3_segs
        nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        num_dell_pos, num_delr_pos, num_D_genes = generative_model.PdelDldelDr_given_D.shape
        
        #These arrays only include the nt identity information, not the PdelDldelDr_given_D info
        PD_nt_pos_vec = [[]]*num_D_genes
        PD_2nd_nt_pos_per_aa_vec = [[]]*num_D_genes
        for D_in in range(num_D_genes):
           
            current_PD_nt_pos_vec = np.zeros((4, len(cutD_genomic_CDR3_segs[D_in])))
            current_PD_2nd_nt_pos_per_aa_vec = {}
            for aa in self.codons_dict.keys():
                current_PD_2nd_nt_pos_per_aa_vec[aa] = np.zeros((4, len(cutD_genomic_CDR3_segs[D_in])))
            
            for pos, nt in enumerate(cutD_genomic_CDR3_segs[D_in]):
                current_PD_nt_pos_vec[nt2num[nt], pos] = 1
                for ins_nt in 'ACGT':
                    for aa in self.codons_dict.keys():
                        if ins_nt + cutD_genomic_CDR3_segs[D_in][pos:pos+2] in self.codons_dict[aa]:
                            current_PD_2nd_nt_pos_per_aa_vec[aa][nt2num[ins_nt], pos] = 1
                            
            PD_nt_pos_vec[D_in] = current_PD_nt_pos_vec
            PD_2nd_nt_pos_per_aa_vec[D_in] = current_PD_2nd_nt_pos_per_aa_vec
        
        min_delDl_given_DdelDr = [[]]*num_D_genes
        max_delDl_given_DdelDr = [[]]*num_D_genes
        zeroD_given_D = [[]]*num_D_genes
        for D_in in range(num_D_genes):
            current_min_delDl_given_delDr = [0]*num_delr_pos
            current_max_delDl_given_delDr = [0]*num_delr_pos
            current_zeroD = 0
            for delr in range(num_delr_pos):
                
                if num_dell_pos > len(cutD_genomic_CDR3_segs[D_in])-delr:
                    current_zeroD += generative_model.PdelDldelDr_given_D[len(cutD_genomic_CDR3_segs[D_in])-delr, delr, D_in]
                
                dell = 0
                while generative_model.PdelDldelDr_given_D[dell, delr, D_in]==0 and dell<num_dell_pos-1:
                    dell+=1
                if generative_model.PdelDldelDr_given_D[dell, delr, D_in] == 0:
                    current_min_delDl_given_delDr[delr] = -1
                else:
                    current_min_delDl_given_delDr[delr] = dell
                if current_min_delDl_given_delDr[delr] == -1:
                    current_max_delDl_given_delDr[delr] = -1
                else:
                    dell = num_dell_pos-1
                    while generative_model.PdelDldelDr_given_D[dell, delr, D_in]==0 and dell>=0:
                        dell -= 1
                    if generative_model.PdelDldelDr_given_D[dell, delr, D_in] == 0:
                        current_max_delDl_given_delDr[delr] = -1
                    else:
                        current_max_delDl_given_delDr[delr] = dell
                
            min_delDl_given_DdelDr[D_in] = current_min_delDl_given_delDr
            max_delDl_given_DdelDr[D_in] = current_max_delDl_given_delDr
            zeroD_given_D[D_in] = current_zeroD
        
        self.PD_nt_pos_vec = PD_nt_pos_vec
        self.PD_2nd_nt_pos_per_aa_vec = PD_2nd_nt_pos_per_aa_vec
        self.min_delDl_given_DdelDr = min_delDl_given_DdelDr 
        self.max_delDl_given_DdelDr = max_delDl_given_DdelDr
        self.zeroD_given_D = zeroD_given_D
        
    def generate_PJdelJ_nt_pos_vecs(self, generative_model, genomic_data):
        """Process P(J)*P(delJ|J) into Pi arrays.
        
        Sets the attributes PJdelJ_nt_pos_vec and PJdelJ_2nd_nt_pos_per_aa_vec.
    
        Parameters
        ----------
        generative_model : GenerativeModelVDJ
            VDJ generative model class containing the model parameters.            
        genomic_data : GenomicDataVDJ
            VDJ genomic data class containing the V, D, and J germline 
            sequences and info.
        
        """
    
        cutJ_genomic_CDR3_segs = genomic_data.cutJ_genomic_CDR3_segs
        nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        num_del_pos = generative_model.PdelJ_given_J.shape[0]
    
        num_D_genes, num_J_genes = generative_model.PDJ.shape
        PJ = np.sum(generative_model.PDJ, axis = 0)
    
        PJdelJ_nt_pos_vec = [[]]*num_J_genes
        PJdelJ_2nd_nt_pos_per_aa_vec = [[]]*num_J_genes
        for J_in, pj in enumerate(PJ):
            #We include the marginal PJ here
            current_PJdelJ_nt_pos_vec = np.zeros((4, len(cutJ_genomic_CDR3_segs[J_in])))
            current_PJdelJ_2nd_nt_pos_per_aa_vec  = {}
            for aa in self.codons_dict.keys():
                current_PJdelJ_2nd_nt_pos_per_aa_vec[aa] = np.zeros((4, len(cutJ_genomic_CDR3_segs[J_in])))
    
            for pos, nt in enumerate(cutJ_genomic_CDR3_segs[J_in]):
                if pos >=  num_del_pos:
                    continue
                if (len(cutJ_genomic_CDR3_segs[J_in]) - pos)%3 == 1: #Start of a codon
                    current_PJdelJ_nt_pos_vec[nt2num[nt], pos] = pj*generative_model.PdelJ_given_J[pos, J_in]
                elif (len(cutJ_genomic_CDR3_segs[J_in]) - pos)%3 == 2: #Mid codon position
                    for ins_nt in 'ACGT':
                        #We need to find what possible codons are allowed for any aa (or motif)
                        for aa in self.codons_dict.keys():
                            if ins_nt + cutJ_genomic_CDR3_segs[J_in][pos:pos+2] in self.codons_dict[aa]:
                                current_PJdelJ_2nd_nt_pos_per_aa_vec[aa][nt2num[ins_nt], pos] = pj*generative_model.PdelJ_given_J[pos, J_in]
                                
                elif (len(cutJ_genomic_CDR3_segs[J_in]) - pos)%3 == 0: #End  of codon
                    current_PJdelJ_nt_pos_vec[0, pos] = pj*generative_model.PdelJ_given_J[pos, J_in]
            PJdelJ_nt_pos_vec[J_in] = current_PJdelJ_nt_pos_vec
            PJdelJ_2nd_nt_pos_per_aa_vec[J_in] = current_PJdelJ_2nd_nt_pos_per_aa_vec
        
        self.PJdelJ_nt_pos_vec = PJdelJ_nt_pos_vec
        self.PJdelJ_2nd_nt_pos_per_aa_vec = PJdelJ_2nd_nt_pos_per_aa_vec
        
    def generate_VD_junction_transfer_matrices(self):
        """Compute the transfer matrices for the VD junction.
        
        Sets the attributes Tvd, Svd, Dvd, lTvd, and lDvd.
        
        """  
        
        nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}                
        
        #Compute Tvd
        Tvd = {}
        for aa in self.codons_dict.keys():
            current_Tvd = np.zeros((4, 4))
            for init_nt in 'ACGT':
                for codon in self.codons_dict[aa]:
                    current_Tvd[nt2num[codon[2]], nt2num[init_nt]] += self.Rvd[nt2num[codon[2]],nt2num[codon[1]]]*self.Rvd[nt2num[codon[1]],nt2num[codon[0]]] * self.Rvd[nt2num[codon[0]],nt2num[init_nt]]
            Tvd[aa] = current_Tvd
            
        #Compute Svd
        Svd = {}
        for aa in self.codons_dict.keys():
            current_Svd = np.zeros((4, 4))
            for ins_nt in 'ACGT':
                if any([codon.startswith(ins_nt) for codon in self.codons_dict[aa]]):
                    current_Svd[nt2num[ins_nt], :] = self.Rvd[nt2num[ins_nt], :]
                
            Svd[aa] = current_Svd
        
        #Compute Dvd                
        Dvd = {}
        for aa in self.codons_dict.keys():
            current_Dvd = np.zeros((4, 4))
            for init_nt in 'ACGT':
                for codon in self.codons_dict[aa]:
                    current_Dvd[nt2num[codon[2]], nt2num[init_nt]] += self.Rvd[nt2num[codon[1]],nt2num[codon[0]]] * self.Rvd[nt2num[codon[0]],nt2num[init_nt]]
            Dvd[aa] = current_Dvd
     

        #Compute lTvd
        lTvd = {}
        for aa in self.codons_dict.keys():
            current_lTvd = np.zeros((4, 4))
            for codon in self.codons_dict[aa]:
                current_lTvd[nt2num[codon[2]], nt2num[codon[0]]] += self.Rvd[nt2num[codon[2]],nt2num[codon[1]]]*self.first_nt_bias_insVD[nt2num[codon[1]]]
            lTvd[aa] = current_lTvd

        
        #Compute lDvd
        lDvd = {}
        for aa in self.codons_dict.keys():
            current_lDvd = np.zeros((4, 4))
            for codon in self.codons_dict[aa]:
                current_lDvd[nt2num[codon[2]], nt2num[codon[0]]] += self.first_nt_bias_insVD[nt2num[codon[1]]]
            lDvd[aa] = current_lDvd
        
        #Set the attributes
        self.Tvd = Tvd
        self.Svd = Svd
        self.Dvd = Dvd
        self.lTvd = lTvd
        self.lDvd = lDvd
        
    def generate_DJ_junction_transfer_matrices(self):
        """Compute the transfer matrices for the VD junction.
        
        Sets the attributes Tdj, Sdj, Ddj, rTdj, and rDdj.
        
        """   
        
        nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

        #Compute Tdj    
        Tdj = {}
        for aa in self.codons_dict.keys():
            current_Tdj = np.zeros((4, 4))
            for init_nt in 'ACGT':
                for codon in self.codons_dict[aa]:
                    current_Tdj[nt2num[codon[0]], nt2num[init_nt]] += self.Rdj[nt2num[codon[0]],nt2num[codon[1]]]*self.Rdj[nt2num[codon[1]],nt2num[codon[2]]] * self.Rdj[nt2num[codon[2]],nt2num[init_nt]]
            Tdj[aa] = current_Tdj
        
        #Compute Sdj
        Sdj = {}
        for aa in self.codons_dict.keys():
            current_Sdj = np.zeros((4, 4))
            for ins_nt in 'ACGT':
                if any([codon.endswith(ins_nt) for codon in self.codons_dict[aa]]):
                    current_Sdj[nt2num[ins_nt], :] = self.Rdj[nt2num[ins_nt], :]    
            Sdj[aa] = current_Sdj
        
        #Compute Ddj
        Ddj = {}
        for aa in self.codons_dict.keys():
            current_Ddj = np.zeros((4, 4))
            for init_nt in 'ACGT':
                for codon in self.codons_dict[aa]:
                    current_Ddj[nt2num[codon[0]], nt2num[init_nt]] += self.Rdj[nt2num[codon[1]],nt2num[codon[2]]] * self.Rdj[nt2num[codon[2]],nt2num[init_nt]]
            Ddj[aa] = current_Ddj
        
        #Compute rTdj
        rTdj = {}
        for aa in self.codons_dict.keys():
            current_lTdj = np.zeros((4, 4))
            for codon in self.codons_dict[aa]:
                current_lTdj[nt2num[codon[0]], nt2num[codon[2]]] += self.Rdj[nt2num[codon[0]],nt2num[codon[1]]]*self.first_nt_bias_insDJ[nt2num[codon[1]]]
            rTdj[aa] = current_lTdj
        
        #Compute rDdj
        rDdj = {}
        for aa in self.codons_dict.keys():
            current_rDdj = np.zeros((4, 4))
            for codon in self.codons_dict[aa]:
                current_rDdj[nt2num[codon[0]], nt2num[codon[2]]] += self.first_nt_bias_insDJ[nt2num[codon[1]]]
            rDdj[aa] = current_rDdj
    
        #Set the attributes
        self.Tdj = Tdj
        self.Sdj = Sdj
        self.Ddj = Ddj
        self.rTdj = rTdj
        self.rDdj = rDdj


#%% PreprocessedParametersVJ Class def
class PreprocessedParametersVJ(PreprocessedParameters):
    """Class used to preprocess the parameters of VJ recombination models.

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
        """Initialize PreprocessedParametersVJ
        
        This intialization computes all of the attributes that will be needed
        for Pgen computation of CDR3 sequences generated by VJ recombination.
        
        Parameters
        ----------
        generative_model : GenerativeModelVJ
            VJ generative model class containing the model parameters.            
        genomic_data : GenomicDataVJ
            VJ genomic data class containing the V and J germline 
            sequences and info.            
        alphabet_file : str, optional
            File name (full pathing from current directory) for a custom alphabet
            definition. If no file is provided, the default alphabet is used, i.e. 
            standard amino acids, undetermined amino acids (B, J, X, and Z), and
            single codon symbols.
        
        """
        PreprocessedParameters.__init__(self, generative_model, genomic_data, alphabet_file)
        
        #Check that the generative_model and genomic_data are VDJ
        if not all([generative_model.__class__.__name__.endswith('VJ'),genomic_data.__class__.__name__.endswith('VJ')]):
            raise ValueError #Need VJ model and data
            

        self.PVJ = generative_model.PVJ
        
        #Set the V genomic Pi arrays
        self.PVdelV_nt_pos_vec = None
        self.PVdelV_2nd_nt_pos_per_aa_vec = None
        self.generate_PVdelV_nt_pos_vecs(generative_model, genomic_data) #Computes and sets the above attributes
        
        #Set the J genomic Pi arrays
        self.PJdelJ_nt_pos_vec = None
        self.PJdelJ_2nd_nt_pos_per_aa_vec = None
        self.generate_PJdelJ_nt_pos_vecs(generative_model, genomic_data) #Computes and sets the above attributes
    
        #Trim, then zeropad PinsVJ
        self.PinsVJ = np.append(np.trim_zeros(generative_model.PinsVJ, 'b'),  [0., 0., 0., 0.])
        
        
        #insVJ Dinucleotide bias transition matrix
        #Note, this used to be called 'RnucleotideVJ_per_nucleotideVJ_5prime'
        self.Rvj = generative_model.Rvj
        
        #Check if first insertion nt probability distribution is given. 
        #Note, this is not normally inferred by IGoR. If the nt bias dist isn't
        #present, use the steady-state distributions defined from Rvj.
        #When not using the steady-state distributions we need to compute the
        #nt bias for the position previous to first_nt_bias (R^{-1}p_0)
        
        if generative_model.first_nt_bias_insVJ == None:
            self.first_nt_bias_insVJ = calc_steady_state_dist(self.Rvj)
            #In steady-state zero_nt_bias_insVJ = first_nt_bias_insVJ so no need to recalculate
            self.zero_nt_bias_insVJ = self.first_nt_bias_insVJ
        else:
            self.first_nt_bias_insVJ = generative_model.first_nt_bias_insVJ
            #Require Rvj*zero_nt_bias_insVJ = first_nt_bias_insVJ  --- can use pseudo-inverse
            self.zero_nt_bias_insVJ = np.dot(np.linalg.pinv(self.Rvj), generative_model.first_nt_bias_insVJ)
            
        #Compute VJ insertion transfer matrices
        self.Tvj = None
        self.Svj = None
        self.Dvj = None
        self.lTvj = None
        self.lDvj = None
        self.generate_VJ_junction_transfer_matrices() #Computes and sets the transfer matrices
    
    def generate_PVdelV_nt_pos_vecs(self, generative_model, genomic_data):
        """Process P(delV|V) into Pi arrays.
        
        Set the attributes PVdelV_nt_pos_vec and PVdelV_2nd_nt_pos_per_aa_vec.
    
        Parameters
        ----------
        generative_model : GenerativeModelVJ
            VJ generative model class containing the model parameters.            
        genomic_data : GenomicDataVJ
            VJ genomic data class containing the V and J germline 
            sequences and info.
        
        """
    
        cutV_genomic_CDR3_segs = genomic_data.cutV_genomic_CDR3_segs
        nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        num_del_pos = generative_model.PdelV_given_V.shape[0]
        num_V_genes = generative_model.PdelV_given_V.shape[1]
        PVdelV_nt_pos_vec = [[]]*num_V_genes
        PVdelV_2nd_nt_pos_per_aa_vec = [[]]*num_V_genes
        for V_in in range(num_V_genes):
            current_PVdelV_nt_pos_vec = np.zeros((4, len(cutV_genomic_CDR3_segs[V_in])))
            current_PVdelV_2nd_nt_pos_per_aa_vec = {}
            for aa in self.codons_dict.keys():
                current_PVdelV_2nd_nt_pos_per_aa_vec[aa] = np.zeros((4, len(cutV_genomic_CDR3_segs[V_in])))
            for pos, nt in enumerate(cutV_genomic_CDR3_segs[V_in]):
                if len(cutV_genomic_CDR3_segs[V_in]) - pos >  num_del_pos:
                    continue
                if pos%3 == 0: #Start of a codon
                    current_PVdelV_nt_pos_vec[nt2num[nt], pos] = generative_model.PdelV_given_V[len(cutV_genomic_CDR3_segs[V_in])-pos-1, V_in]    
                elif pos%3 == 1: #Mid codon position
                    for ins_nt in 'ACGT':
                        #We need to find what possible codons are allowed for any aa (or motif)
                        for aa in self.codons_dict.keys():
                            if cutV_genomic_CDR3_segs[V_in][pos-1:pos+1]+ ins_nt in self.codons_dict[aa]:
                                current_PVdelV_2nd_nt_pos_per_aa_vec[aa][nt2num[ins_nt], pos] = generative_model.PdelV_given_V[len(cutV_genomic_CDR3_segs[V_in])-pos-1, V_in]            
                elif pos%3 == 2: #End of codon
                    current_PVdelV_nt_pos_vec[0, pos] = generative_model.PdelV_given_V[len(cutV_genomic_CDR3_segs[V_in])-pos-1, V_in]
            PVdelV_nt_pos_vec[V_in] = current_PVdelV_nt_pos_vec
            PVdelV_2nd_nt_pos_per_aa_vec[V_in] = current_PVdelV_2nd_nt_pos_per_aa_vec
    
        
        self.PVdelV_nt_pos_vec = PVdelV_nt_pos_vec
        self.PVdelV_2nd_nt_pos_per_aa_vec = PVdelV_2nd_nt_pos_per_aa_vec

    def generate_PJdelJ_nt_pos_vecs(self, generative_model, genomic_data):
        """Process P(delJ|J) into Pi arrays.
        
        Set the attributes PJdelJ_nt_pos_vec and PJdelJ_2nd_nt_pos_per_aa_vec.
    
        Parameters
        ----------
        generative_model : GenerativeModelVJ
            VJ generative model class containing the model parameters.            
        genomic_data : GenomicDataVJ
            VJ genomic data class containing the V and J germline 
            sequences and info.
        
        """
        
        cutJ_genomic_CDR3_segs = genomic_data.cutJ_genomic_CDR3_segs
        nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
        num_del_pos = generative_model.PdelJ_given_J.shape[0]
    
        num_D_genes, num_J_genes = generative_model.PVJ.shape
    
        PJdelJ_nt_pos_vec = [[]]*num_J_genes
        PJdelJ_2nd_nt_pos_per_aa_vec = [[]]*num_J_genes
        for J_in in range(num_J_genes):
            current_PJdelJ_nt_pos_vec = np.zeros((4, len(cutJ_genomic_CDR3_segs[J_in])))
            current_PJdelJ_2nd_nt_pos_per_aa_vec  = {}
            for aa in self.codons_dict.keys():
                current_PJdelJ_2nd_nt_pos_per_aa_vec[aa] = np.zeros((4, len(cutJ_genomic_CDR3_segs[J_in])))
    
            for pos, nt in enumerate(cutJ_genomic_CDR3_segs[J_in]):
                if pos >=  num_del_pos:
                    continue
                if (len(cutJ_genomic_CDR3_segs[J_in]) - pos)%3 == 1: #Start of a codon
                    current_PJdelJ_nt_pos_vec[nt2num[nt], pos] = generative_model.PdelJ_given_J[pos, J_in]
                elif (len(cutJ_genomic_CDR3_segs[J_in]) - pos)%3 == 2: #Mid codon position
                    for ins_nt in 'ACGT':
                        #We need to find what possible codons are allowed for any aa (or motif)
                        for aa in self.codons_dict.keys():
                            if ins_nt + cutJ_genomic_CDR3_segs[J_in][pos:pos+2] in self.codons_dict[aa]:
                                current_PJdelJ_2nd_nt_pos_per_aa_vec[aa][nt2num[ins_nt], pos] = generative_model.PdelJ_given_J[pos, J_in]
                                
                elif (len(cutJ_genomic_CDR3_segs[J_in]) - pos)%3 == 0: #End  of codon
                    current_PJdelJ_nt_pos_vec[0, pos] = generative_model.PdelJ_given_J[pos, J_in]
            PJdelJ_nt_pos_vec[J_in] = current_PJdelJ_nt_pos_vec
            PJdelJ_2nd_nt_pos_per_aa_vec[J_in] = current_PJdelJ_2nd_nt_pos_per_aa_vec
            
        self.PJdelJ_nt_pos_vec = PJdelJ_nt_pos_vec
        self.PJdelJ_2nd_nt_pos_per_aa_vec = PJdelJ_2nd_nt_pos_per_aa_vec
        
        
    def generate_VJ_junction_transfer_matrices(self):
        """Compute the transfer matrices for the VJ junction.
        
        Sets the attributes Tvj, Svj, Dvj, lTvj, and lDvj.
        
        """    
        
        nt2num = {'A': 0, 'C': 1, 'G': 2, 'T': 3}                
        
        #Compute Tvj
        Tvj = {}
        for aa in self.codons_dict.keys():
            current_Tvj = np.zeros((4, 4))
            for init_nt in 'ACGT':
                for codon in self.codons_dict[aa]:
                    current_Tvj[nt2num[codon[2]], nt2num[init_nt]] += self.Rvj[nt2num[codon[2]],nt2num[codon[1]]]*self.Rvj[nt2num[codon[1]],nt2num[codon[0]]] * self.Rvj[nt2num[codon[0]],nt2num[init_nt]]
            Tvj[aa] = current_Tvj
    
        #Compute Svj
        Svj = {}
        for aa in self.codons_dict.keys():
            current_Svj = np.zeros((4, 4))
            for ins_nt in 'ACGT':
                if any([codon.startswith(ins_nt) for codon in self.codons_dict[aa]]):
                    current_Svj[nt2num[ins_nt], :] = self.Rvj[nt2num[ins_nt], :]             
            Svj[aa] = current_Svj
        
        #Compute Dvj               
        Dvj = {}    
        for aa in self.codons_dict.keys():
            current_Dvj = np.zeros((4, 4))
            for init_nt in 'ACGT':
                for codon in self.codons_dict[aa]:
                    current_Dvj[nt2num[codon[2]], nt2num[init_nt]] += self.Rvj[nt2num[codon[1]],nt2num[codon[0]]] * self.Rvj[nt2num[codon[0]],nt2num[init_nt]]
            Dvj[aa] = current_Dvj

        #Compute lTvj
        lTvj = {}
        for aa in self.codons_dict.keys():
            current_lTvj = np.zeros((4, 4))
            for codon in self.codons_dict[aa]:
                current_lTvj[nt2num[codon[2]], nt2num[codon[0]]] += self.Rvj[nt2num[codon[2]],nt2num[codon[1]]]*self.first_nt_bias_insVJ[nt2num[codon[1]]]
            lTvj[aa] = current_lTvj

        #Compute lDvj        
        lDvj = {}    
        for aa in self.codons_dict.keys():
            current_lDvj = np.zeros((4, 4))
            for codon in self.codons_dict[aa]:
                current_lDvj[nt2num[codon[2]], nt2num[codon[0]]] += self.first_nt_bias_insVJ[nt2num[codon[1]]]
            lDvj[aa] = current_lDvj
    

        #Set the attributes
        self.Tvj = Tvj
        self.Svj = Svj
        self.Dvj = Dvj
        self.lTvj = lTvj
        self.lDvj = lDvj