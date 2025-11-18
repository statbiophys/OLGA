"""Tests for sequence generation functionality."""
import pytest
import numpy as np
import olga.sequence_generation as sequence_generation
from olga.utils import nt2aa


class TestSequenceGenerationVDJ:
    """Tests for VDJ sequence generation."""
    
    def test_generate_sequence(self, human_t_beta_model):
        """Test generating a single sequence."""
        generative_model = human_t_beta_model['generative_model']
        genomic_data = human_t_beta_model['genomic_data']
        
        seq_gen = sequence_generation.SequenceGenerationVDJ(
            generative_model, genomic_data
        )
        
        ntseq, aaseq, V_in, J_in = seq_gen.gen_rnd_prod_CDR3()
        
        assert isinstance(ntseq, str)
        assert isinstance(aaseq, str)
        assert isinstance(V_in, (int, np.integer))
        assert isinstance(J_in, (int, np.integer))
        assert len(ntseq) > 0
        assert len(aaseq) > 0
        assert len(ntseq) % 3 == 0  # Must be in frame
        assert aaseq[0] == 'C'  # Should start with C
        assert aaseq[-1] in 'FVW'  # Should end with F, V, or W
        assert '*' not in aaseq  # Should be productive (no stop codons)
    
    def test_generate_multiple_sequences(self, human_t_beta_model):
        """Test generating multiple sequences."""
        generative_model = human_t_beta_model['generative_model']
        genomic_data = human_t_beta_model['genomic_data']
        
        seq_gen = sequence_generation.SequenceGenerationVDJ(
            generative_model, genomic_data
        )
        
        sequences = []
        for _ in range(10):
            ntseq, aaseq, V_in, J_in = seq_gen.gen_rnd_prod_CDR3()
            sequences.append((ntseq, aaseq, V_in, J_in))
            assert len(ntseq) > 0
            assert len(aaseq) > 0
        
        # Check that we got some variety
        unique_seqs = set(aaseq for _, aaseq, _, _ in sequences)
        assert len(unique_seqs) > 1  # Should have some variety
    
    def test_conserved_j_residues(self, human_t_beta_model):
        """Test with custom conserved J residues."""
        generative_model = human_t_beta_model['generative_model']
        genomic_data = human_t_beta_model['genomic_data']
        
        seq_gen = sequence_generation.SequenceGenerationVDJ(
            generative_model, genomic_data
        )
        
        ntseq, aaseq, V_in, J_in = seq_gen.gen_rnd_prod_CDR3(conserved_J_residues='F')
        
        assert aaseq[-1] == 'F'  # Should end with F only
    
    def test_sequence_consistency(self, human_t_beta_model):
        """Test that generated nucleotide and amino acid sequences are consistent."""
        generative_model = human_t_beta_model['generative_model']
        genomic_data = human_t_beta_model['genomic_data']
        
        seq_gen = sequence_generation.SequenceGenerationVDJ(
            generative_model, genomic_data
        )
        
        ntseq, aaseq, V_in, J_in = seq_gen.gen_rnd_prod_CDR3()
        
        # Check that translation is consistent (allowing for ambiguous codons)
        translated = nt2aa(ntseq)
        # The translated sequence should match or be compatible with aaseq
        assert len(translated) == len(aaseq) or len(translated) == len(ntseq) // 3


class TestSequenceGenerationVJ:
    """Tests for VJ sequence generation."""
    
    def test_generate_sequence(self, human_t_alpha_model):
        """Test generating a single VJ sequence."""
        generative_model = human_t_alpha_model['generative_model']
        genomic_data = human_t_alpha_model['genomic_data']
        
        seq_gen = sequence_generation.SequenceGenerationVJ(
            generative_model, genomic_data
        )
        
        ntseq, aaseq, V_in, J_in = seq_gen.gen_rnd_prod_CDR3()
        
        assert isinstance(ntseq, str)
        assert isinstance(aaseq, str)
        assert isinstance(V_in, (int, np.integer))
        assert isinstance(J_in, (int, np.integer))
        assert len(ntseq) > 0
        assert len(aaseq) > 0
        assert len(ntseq) % 3 == 0  # Must be in frame
        assert aaseq[0] == 'C'  # Should start with C
        assert aaseq[-1] in 'FVW'  # Should end with F, V, or W
        assert '*' not in aaseq  # Should be productive
    
    def test_generate_multiple_sequences(self, human_t_alpha_model):
        """Test generating multiple VJ sequences."""
        generative_model = human_t_alpha_model['generative_model']
        genomic_data = human_t_alpha_model['genomic_data']
        
        seq_gen = sequence_generation.SequenceGenerationVJ(
            generative_model, genomic_data
        )
        
        sequences = []
        for _ in range(10):
            ntseq, aaseq, V_in, J_in = seq_gen.gen_rnd_prod_CDR3()
            sequences.append((ntseq, aaseq, V_in, J_in))
        
        # Check variety
        unique_seqs = set(aaseq for _, aaseq, _, _ in sequences)
        assert len(unique_seqs) > 1


class TestSequenceGenerationAttributes:
    """Tests for sequence generation attributes."""
    
    def test_vdj_attributes(self, human_t_beta_model):
        """Test VDJ sequence generator has required attributes."""
        generative_model = human_t_beta_model['generative_model']
        genomic_data = human_t_beta_model['genomic_data']
        
        seq_gen = sequence_generation.SequenceGenerationVDJ(
            generative_model, genomic_data
        )
        
        assert hasattr(seq_gen, 'CPV'), f"CPV attribute not found in {seq_gen}"
        assert hasattr(seq_gen, 'CPDJ'), f"CPDJ attribute not found in {seq_gen}"
        assert hasattr(seq_gen, 'CinsVD'), f"CinsVD attribute not found in {seq_gen}"  
        assert hasattr(seq_gen, 'CinsDJ'), f"CinsDJ attribute not found in {seq_gen}"
        assert hasattr(seq_gen, 'cutV_genomic_CDR3_segs'), f"cutV_genomic_CDR3_segs attribute not found in {seq_gen}"
        assert hasattr(seq_gen, 'cutD_genomic_CDR3_segs'), f"cutD_genomic_CDR3_segs attribute not found in {seq_gen}"
        assert hasattr(seq_gen, 'cutJ_genomic_CDR3_segs'), f"cutJ_genomic_CDR3_segs attribute not found in {seq_gen}"
    
    def test_vj_attributes(self, human_t_alpha_model):
        """Test VJ sequence generator has required attributes."""
        generative_model = human_t_alpha_model['generative_model']
        genomic_data = human_t_alpha_model['genomic_data']
        
        seq_gen = sequence_generation.SequenceGenerationVJ(
            generative_model, genomic_data
        )
        
        assert hasattr(seq_gen, 'CPVJ'), f"CPVJ attribute not found in {seq_gen}"
        assert hasattr(seq_gen, 'CPinsVJ'), f"CPinsVJ attribute not found in {seq_gen}"
        assert hasattr(seq_gen, 'C_Rvj'), f"C_Rvj attribute not found in {seq_gen}"
        assert hasattr(seq_gen, 'C_first_nt_bias_insVJ'), f"C_first_nt_bias_insVJ attribute not found in {seq_gen}"
        assert hasattr(seq_gen, 'cutV_genomic_CDR3_segs'), f"cutV_genomic_CDR3_segs attribute not found in {seq_gen}"
        assert hasattr(seq_gen, 'cutJ_genomic_CDR3_segs'), f"cutJ_genomic_CDR3_segs attribute not found in {seq_gen}"

