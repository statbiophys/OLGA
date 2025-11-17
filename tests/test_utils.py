"""Tests for utility functions."""
import pytest
import numpy as np
from olga.utils import (
    nt2aa, nt2codon_rep, cutR_seq, cutL_seq,
    construct_codons_dict, determine_seq_type,
    calc_steady_state_dist, gene_to_num_str,
    calc_S, calc_S_single_gene, calc_S_joint_genes, calc_Sins
)


class TestNt2Aa:
    """Tests for nucleotide to amino acid translation."""
    
    def test_basic_translation(self):
        """Test basic nucleotide to amino acid translation."""
        ntseq = 'TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC'
        expected = 'CAWSVAPDRGGYTF'
        assert nt2aa(ntseq) == expected
    
    def test_lowercase(self):
        """Test that lowercase nucleotides work."""
        ntseq = 'tgtgcctggagtgtagctccggacaggggtggctacaccttc'
        expected = 'CAWSVAPDRGGYTF'
        assert nt2aa(ntseq) == expected
    
    def test_mixed_case(self):
        """Test mixed case nucleotides."""
        ntseq = 'TGTgccTGGagtGTAGCTCCGGACAGGGGTGGCTACACCTTC'
        expected = 'CAWSVAPDRGGYTF'
        assert nt2aa(ntseq) == expected
    
    def test_empty_sequence(self):
        """Test empty sequence."""
        assert nt2aa('') == ''
    
    def test_incomplete_codon(self):
        """Test sequence with incomplete codon at end."""
        ntseq = 'TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTCAA'  # Extra AA
        # Should ignore incomplete codon
        result = nt2aa(ntseq)
        assert len(result) == len(ntseq) // 3


class TestNt2CodonRep:
    """Tests for nucleotide to codon representation."""
    
    def test_basic_conversion(self):
        """Test basic conversion."""
        ntseq = 'TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC'
        result = nt2codon_rep(ntseq)
        assert len(result) == len(ntseq) // 3
        assert isinstance(result, str)
    
    def test_reversibility(self):
        """Test that conversion is consistent."""
        ntseq = 'TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC'
        result1 = nt2codon_rep(ntseq)
        result2 = nt2codon_rep(ntseq)
        assert result1 == result2


class TestCutSequences:
    """Tests for sequence cutting functions."""
    
    def test_cutR_seq_no_cut(self):
        """Test cutR_seq with no cutting."""
        seq = 'TGCGCCAGCAGTGAGTC'
        result = cutR_seq(seq, 0, 4)
        assert len(result) > len(seq)  # Should add palindrome
    
    def test_cutR_seq_with_cut(self):
        """Test cutR_seq with cutting."""
        seq = 'TGCGCCAGCAGTGAGTC'
        result = cutR_seq(seq, 8, 4)
        assert len(result) < len(seq)  # Should cut
    
    def test_cutL_seq_no_cut(self):
        """Test cutL_seq with no cutting."""
        seq = 'TGAACACTGAAGCTTTCTTT'
        result = cutL_seq(seq, 0, 4)
        assert len(result) > len(seq)  # Should add palindrome
    
    def test_cutL_seq_with_cut(self):
        """Test cutL_seq with cutting."""
        seq = 'TGAACACTGAAGCTTTCTTT'
        result = cutL_seq(seq, 8, 4)
        assert len(result) < len(seq)  # Should cut


class TestConstructCodonsDict:
    """Tests for codon dictionary construction."""
    
    def test_default_dict(self):
        """Test default codon dictionary."""
        codons_dict = construct_codons_dict()
        assert 'A' in codons_dict
        assert 'C' in codons_dict
        assert 'X' in codons_dict  # Ambiguous amino acid
        assert len(codons_dict) > 20  # Should have standard AAs + ambiguous + codon symbols
    
    def test_custom_alphabet(self, tmp_path):
        """Test custom alphabet file."""
        alphabet_file = tmp_path / 'alphabet.txt'
        alphabet_file.write_text('^: A, G, R\n')
        
        codons_dict = construct_codons_dict(str(alphabet_file))
        assert '^' in codons_dict


class TestDetermineSeqType:
    """Tests for sequence type determination."""
    
    def test_ntseq(self):
        """Test nucleotide sequence detection."""
        aa_alphabet = ''.join(construct_codons_dict().keys())
        seq = 'TGTGCCTGGAGTGTAGCTCCGGACAGGGGTGGCTACACCTTC'
        assert determine_seq_type(seq, aa_alphabet) == 'ntseq'
    
    def test_aaseq(self):
        """Test amino acid sequence detection."""
        aa_alphabet = ''.join(construct_codons_dict().keys())
        seq = 'CAWSVAPDRGGYTF'
        assert determine_seq_type(seq, aa_alphabet) == 'aaseq'
    
    def test_regex(self):
        """Test regex sequence detection."""
        aa_alphabet = ''.join(construct_codons_dict().keys())
        seq = 'CASSAX{0,5}SARPEQFF'
        assert determine_seq_type(seq, aa_alphabet) == 'regex'


class TestSteadyState:
    """Tests for steady state distribution calculation."""
    
    def test_steady_state(self):
        """Test steady state calculation."""
        R = np.array([[0.25, 0.25, 0.25, 0.25],
                      [0.25, 0.25, 0.25, 0.25],
                      [0.25, 0.25, 0.25, 0.25],
                      [0.25, 0.25, 0.25, 0.25]])
        p_ss = calc_steady_state_dist(R)
        assert len(p_ss) == 4
        assert np.allclose(np.sum(p_ss), 1.0)
        assert np.allclose(p_ss, 0.25)


class TestGeneToNumStr:
    """Tests for gene name conversion."""
    
    def test_basic_conversion(self):
        """Test basic gene name conversion."""
        result = gene_to_num_str('TRBV30*01', 'V')
        assert isinstance(result, str)
        assert 'v' in result.lower()
    
    def test_different_gene_types(self):
        """Test different gene types."""
        v_result = gene_to_num_str('TRBV30*01', 'V')
        j_result = gene_to_num_str('TRBJ1-2*01', 'J')
        assert v_result != j_result


class TestEntropyFunctions:
    """Tests for entropy calculation functions."""
    
    def test_calc_S(self):
        """Test basic entropy calculation."""
        P = np.array([0.5, 0.5])
        S = calc_S(P, 2)
        assert S == 1.0  # Entropy of fair coin flip
    
    def test_calc_S_single_gene(self):
        """Test single gene entropy calculation."""
        PG = np.array([0.5, 0.5])
        PdelG_given_G = np.array([[0.5, 0.5], [0.5, 0.5]])
        SG, SdelG = calc_S_single_gene(PG, PdelG_given_G, base=2)
        assert SG > 0
        assert SdelG > 0
    
    def test_calc_S_joint_genes(self):
        """Test joint gene entropy calculation."""
        PG1G2 = np.array([[0.25, 0.25], [0.25, 0.25]])
        PdelG1 = np.array([[0.5, 0.5], [0.5, 0.5]])
        PdelG2 = np.array([[0.5, 0.5], [0.5, 0.5]])
        SG1G2, SdelG1, SdelG2 = calc_S_joint_genes(PG1G2, PdelG1, PdelG2, base=2)
        assert SG1G2 > 0
        assert SdelG1 > 0
        assert SdelG2 > 0
    
    def test_calc_Sins(self):
        """Test insertion entropy calculation."""
        Pins = np.array([0.1, 0.2, 0.3, 0.4])
        Pins = Pins / np.sum(Pins)  # Normalize
        R = np.array([[0.25, 0.25, 0.25, 0.25],
                      [0.25, 0.25, 0.25, 0.25],
                      [0.25, 0.25, 0.25, 0.25],
                      [0.25, 0.25, 0.25, 0.25]])
        Sins, Smarkov = calc_Sins(Pins, R, base=2)
        assert Sins > 0
        assert Smarkov > 0

