"""Tests for FastPgen (numba-accelerated) functionality."""
import pytest
import numpy as np
from olga.performance.fast_pgen import FastPgen


class TestFastPgenVDJ:
    """Tests for FastPgen with VDJ models."""
    
    def test_fast_pgen_initialization(self, human_t_beta_model):
        """Test FastPgen initialization with VDJ model."""
        pgen_model = human_t_beta_model['pgen_model']
        
        fast_model = FastPgen(pgen_model)
        
        assert fast_model._impl == pgen_model
        assert hasattr(fast_model, 'codons_dict')
    
    def test_fast_pgen_compute_aa(self, human_t_beta_model):
        """Test FastPgen amino acid computation."""
        pgen_model = human_t_beta_model['pgen_model']
        fast_model = FastPgen(pgen_model)
        
        cdr3_seq = 'CASSTGQANYGYTF'
        
        pgen_slow = pgen_model.compute_aa_CDR3_pgen(cdr3_seq)
        pgen_fast = fast_model.compute_aa_CDR3_pgen(cdr3_seq)
        
        # Results should be very close (within numerical precision)
        assert np.allclose(pgen_slow, pgen_fast, rtol=1e-6)
        assert pgen_fast > 0
    
    def test_fast_pgen_compute_nt(self, human_t_beta_model):
        """Test FastPgen nucleotide computation."""
        pgen_model = human_t_beta_model['pgen_model']
        fast_model = FastPgen(pgen_model)
        
        ntseq = 'TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT'
        
        pgen_slow = pgen_model.compute_nt_CDR3_pgen(ntseq, print_warnings=False)
        pgen_fast = fast_model.compute_nt_CDR3_pgen(ntseq, print_warnings=False)
        
        # Results should be very close
        assert np.allclose(pgen_slow, pgen_fast, rtol=1e-6)
    
    def test_fast_pgen_with_masks(self, human_t_beta_model):
        """Test FastPgen with V and J masks."""
        pgen_model = human_t_beta_model['pgen_model']
        fast_model = FastPgen(pgen_model)
        
        cdr3_seq = 'CASSTGQANYGYTF'
        
        # Get some V and J genes
        v_genes = list(pgen_model.V_allele_names[:2]) if len(pgen_model.V_allele_names) >= 2 else pgen_model.V_allele_names
        j_genes = list(pgen_model.J_allele_names[:2]) if len(pgen_model.J_allele_names) >= 2 else pgen_model.J_allele_names
        
        pgen_slow = pgen_model.compute_aa_CDR3_pgen(cdr3_seq, v_genes, j_genes)
        pgen_fast = fast_model.compute_aa_CDR3_pgen(cdr3_seq, v_genes, j_genes)
        
        assert np.allclose(pgen_slow, pgen_fast, rtol=1e-6)


class TestFastPgenVJ:
    """Tests for FastPgen with VJ models."""
    
    def test_fast_pgen_initialization(self, human_t_alpha_model):
        """Test FastPgen initialization with VJ model."""
        pgen_model = human_t_alpha_model['pgen_model']
        
        fast_model = FastPgen(pgen_model)
        
        assert fast_model._impl == pgen_model
    
    def test_fast_pgen_compute_aa(self, human_t_alpha_model):
        """Test FastPgen amino acid computation for VJ model."""
        pgen_model = human_t_alpha_model['pgen_model']
        fast_model = FastPgen(pgen_model)
        
        cdr3_seq = 'CAVKIQGAQKLVF'
        
        pgen_slow = pgen_model.compute_aa_CDR3_pgen(cdr3_seq)
        pgen_fast = fast_model.compute_aa_CDR3_pgen(cdr3_seq)
        
        assert np.allclose(pgen_slow, pgen_fast, rtol=1e-6)
        assert pgen_fast > 0


class TestFastPgenPickling:
    """Tests for FastPgen pickling (for multiprocessing)."""
    
    def test_fast_pgen_pickling(self, human_t_beta_model):
        """Test that FastPgen can be pickled and unpickled."""
        import pickle
        
        pgen_model = human_t_beta_model['pgen_model']
        fast_model = FastPgen(pgen_model)
        
        # Test that it can be pickled
        pickled = pickle.dumps(fast_model)
        unpickled = pickle.loads(pickled)
        
        # Test that unpickled model still works
        cdr3_seq = 'CASSTGQANYGYTF'
        pgen_original = fast_model.compute_aa_CDR3_pgen(cdr3_seq)
        pgen_unpickled = unpickled.compute_aa_CDR3_pgen(cdr3_seq)
        
        assert np.allclose(pgen_original, pgen_unpickled, rtol=1e-6)


class TestFastPgenConsistency:
    """Tests for FastPgen consistency with standard implementation."""
    
    def test_multiple_sequences_consistency(self, human_t_beta_model):
        """Test that FastPgen gives consistent results for multiple sequences."""
        pgen_model = human_t_beta_model['pgen_model']
        fast_model = FastPgen(pgen_model)
        
        test_sequences = [
            'CASSTGQANYGYTF',
            'CAWSVAPDRGGYTF',
            'CASSDAQGRNRGTEAFF'
        ]
        
        for seq in test_sequences:
            pgen_slow = pgen_model.compute_aa_CDR3_pgen(seq, print_warnings=False)
            pgen_fast = fast_model.compute_aa_CDR3_pgen(seq, print_warnings=False)
            
            if pgen_slow > 0:  # Only check if slow version gives valid result
                assert np.allclose(pgen_slow, pgen_fast, rtol=1e-6)

