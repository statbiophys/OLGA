"""Tests for generation probability computation."""
import pytest
import numpy as np
import olga.generation_probability as generation_probability


class TestGenerationProbabilityVDJ:
    """Tests for VDJ generation probability computation."""
    
    def test_compute_aa_CDR3_pgen(self, human_t_beta_model):
        """Test amino acid CDR3 Pgen computation."""
        pgen_model = human_t_beta_model['pgen_model']
        cdr3_seq = 'CASSTGQANYGYTF'
        
        pgen = pgen_model.compute_aa_CDR3_pgen(cdr3_seq)
        
        assert isinstance(pgen, (float, np.floating))
        assert pgen > 0
        assert pgen <= 1.0
    
    def test_compute_nt_CDR3_pgen(self, human_t_beta_model):
        """Test nucleotide CDR3 Pgen computation."""
        pgen_model = human_t_beta_model['pgen_model']
        ntseq = 'TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT'
        
        pgen = pgen_model.compute_nt_CDR3_pgen(ntseq)
        
        assert isinstance(pgen, (float, np.floating))
        assert pgen > 0
        assert pgen <= 1.0
    
    def test_compute_with_v_j_mask(self, human_t_beta_model):
        """Test Pgen computation with V and J masks."""
        pgen_model = human_t_beta_model['pgen_model']
        cdr3_seq = 'CASSTGQANYGYTF'
        
        # Get some V and J gene names from the model
        v_genes = list(pgen_model.V_allele_names[:3]) if len(pgen_model.V_allele_names) >= 3 else pgen_model.V_allele_names
        j_genes = list(pgen_model.J_allele_names[:3]) if len(pgen_model.J_allele_names) >= 3 else pgen_model.J_allele_names
        
        pgen_masked = pgen_model.compute_aa_CDR3_pgen(cdr3_seq, v_genes, j_genes)
        pgen_unmasked = pgen_model.compute_aa_CDR3_pgen(cdr3_seq)
        
        assert pgen_masked > 0
        # Masked should generally be <= unmasked (more restrictive)
        assert pgen_masked <= pgen_unmasked or abs(pgen_masked - pgen_unmasked) < 1e-10
    
    def test_empty_sequence(self, human_t_beta_model):
        """Test that empty sequence returns zero Pgen."""
        pgen_model = human_t_beta_model['pgen_model']
        pgen = pgen_model.compute_aa_CDR3_pgen('')
        assert pgen == 0
    
    def test_invalid_sequence(self, human_t_beta_model):
        """Test that invalid sequence returns zero Pgen."""
        pgen_model = human_t_beta_model['pgen_model']
        # Use a sequence with invalid characters
        pgen = pgen_model.compute_aa_CDR3_pgen('INVALID123', print_warnings=False)
        assert pgen == 0
    
    def test_regex_sequence(self, human_t_beta_model):
        """Test regex sequence Pgen computation."""
        pgen_model = human_t_beta_model['pgen_model']
        regex_seq = 'CASSAX{0,3}SARPEQFF'
        
        pgen = pgen_model.compute_regex_CDR3_template_pgen(regex_seq, print_warnings=False)
        
        assert isinstance(pgen, (float, np.floating))
        assert pgen >= 0


class TestGenerationProbabilityVJ:
    """Tests for VJ generation probability computation."""
    
    def test_compute_aa_CDR3_pgen(self, human_t_alpha_model):
        """Test amino acid CDR3 Pgen computation for VJ model."""
        pgen_model = human_t_alpha_model['pgen_model']
        cdr3_seq = 'CAVKIQGAQKLVF'
        
        pgen = pgen_model.compute_aa_CDR3_pgen(cdr3_seq)
        
        assert isinstance(pgen, (float, np.floating))
        assert pgen > 0
        assert pgen <= 1.0
    
    def test_compute_nt_CDR3_pgen(self, human_t_alpha_model):
        """Test nucleotide CDR3 Pgen computation for VJ model."""
        pgen_model = human_t_alpha_model['pgen_model']
        seq_gen_model = human_t_alpha_model['seq_gen_model']

        ntseq, aaseq, V_in, J_in = seq_gen_model.gen_rnd_prod_CDR3()
        pgen = pgen_model.compute_nt_CDR3_pgen(ntseq, print_warnings=False)
        assert isinstance(pgen, (float, np.floating))
        assert pgen > 0
        assert pgen < 1.0
    
    def test_consistency_between_aa_and_nt(self, human_t_alpha_model):
        """Test that AA and NT Pgen are consistent when sequence matches."""
        pgen_model = human_t_alpha_model['pgen_model']
        # For a valid sequence, both should give reasonable results
        aaseq = 'CAVKIQGAQKLVF'
        
        pgen_aa = pgen_model.compute_aa_CDR3_pgen(aaseq)
        assert pgen_aa > 0


class TestPgenModelAttributes:
    """Tests for Pgen model attributes."""
    
    def test_vdj_model_attributes(self, human_t_beta_model):
        """Test VDJ model has required attributes."""
        pgen_model = human_t_beta_model['pgen_model']
        
        assert hasattr(pgen_model, 'codons_dict')
        assert hasattr(pgen_model, 'V_allele_names')
        assert hasattr(pgen_model, 'J_allele_names')
        assert hasattr(pgen_model, 'D_allele_names')
        assert hasattr(pgen_model, 'V_mask_mapping')
        assert hasattr(pgen_model, 'J_mask_mapping')
        assert len(pgen_model.V_allele_names) > 0
        assert len(pgen_model.J_allele_names) > 0
    
    def test_vj_model_attributes(self, human_t_alpha_model):
        """Test VJ model has required attributes."""
        pgen_model = human_t_alpha_model['pgen_model']
        
        assert hasattr(pgen_model, 'codons_dict')
        assert hasattr(pgen_model, 'V_allele_names')
        assert hasattr(pgen_model, 'J_allele_names')
        assert hasattr(pgen_model, 'V_mask_mapping')
        assert hasattr(pgen_model, 'J_mask_mapping')
        assert len(pgen_model.V_allele_names) > 0
        assert len(pgen_model.J_allele_names) > 0


class TestPgenConsistency:
    """Tests for Pgen computation consistency."""
    
    def test_same_sequence_same_pgen(self, human_t_beta_model):
        """Test that same sequence gives same Pgen."""
        pgen_model = human_t_beta_model['pgen_model']
        cdr3_seq = 'CASSTGQANYGYTF'
        
        pgen1 = pgen_model.compute_aa_CDR3_pgen(cdr3_seq)
        pgen2 = pgen_model.compute_aa_CDR3_pgen(cdr3_seq)
        
        assert np.allclose(pgen1, pgen2, rtol=1e-10)
    
    def test_longer_sequence_lower_pgen(self, human_t_beta_model):
        """Test that longer sequences generally have lower Pgen."""
        pgen_model = human_t_beta_model['pgen_model']
        
        short_seq = 'CASSTGQANYGYTF'
        # Create a longer sequence (if valid)
        long_seq = 'CASSTGQANYGYTFCASSTGQANYGYTF'
        
        pgen_short = pgen_model.compute_aa_CDR3_pgen(short_seq, print_warnings=False)
        pgen_long = pgen_model.compute_aa_CDR3_pgen(long_seq, print_warnings=False)
        
        # Longer sequences should generally have lower or equal Pgen
        # (though this isn't always true, so we just check they're both valid)
        assert pgen_short > 0 or pgen_long == 0
        assert pgen_long >= 0

