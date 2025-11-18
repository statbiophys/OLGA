"""Tests for model loading functionality."""
import os
import pytest
import numpy as np
import olga.load_model as load_model


class TestGenomicDataVDJ:
    """Tests for VDJ genomic data loading."""
    
    def test_load_human_t_beta(self, model_dir):
        """Test loading human T beta genomic data."""
        model_folder = os.path.join(model_dir, 'human_T_beta')
        params_file = os.path.join(model_folder, 'model_params.txt')
        V_anchor_file = os.path.join(model_folder, 'V_gene_CDR3_anchors.csv')
        J_anchor_file = os.path.join(model_folder, 'J_gene_CDR3_anchors.csv')
        
        genomic_data = load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(params_file, V_anchor_file, J_anchor_file)
        
        assert genomic_data.genV is not None
        assert genomic_data.genD is not None
        assert genomic_data.genJ is not None
        assert len(genomic_data.genV) > 0
        assert len(genomic_data.genD) > 0
        assert len(genomic_data.genJ) > 0
        assert genomic_data.cutV_genomic_CDR3_segs is not None
        assert genomic_data.cutD_genomic_CDR3_segs is not None
        assert genomic_data.cutJ_genomic_CDR3_segs is not None


class TestGenomicDataVJ:
    """Tests for VJ genomic data loading."""
    
    def test_load_human_t_alpha(self, model_dir):
        """Test loading human T alpha genomic data."""
        model_folder = os.path.join(model_dir, 'human_T_alpha')
        params_file = os.path.join(model_folder, 'model_params.txt')
        V_anchor_file = os.path.join(model_folder, 'V_gene_CDR3_anchors.csv')
        J_anchor_file = os.path.join(model_folder, 'J_gene_CDR3_anchors.csv')
        
        genomic_data = load_model.GenomicDataVJ()
        genomic_data.load_igor_genomic_data(params_file, V_anchor_file, J_anchor_file)
        
        assert genomic_data.genV is not None
        assert genomic_data.genJ is not None
        assert len(genomic_data.genV) > 0
        assert len(genomic_data.genJ) > 0
        assert genomic_data.cutV_genomic_CDR3_segs is not None
        assert genomic_data.cutJ_genomic_CDR3_segs is not None


class TestGenerativeModelVDJ:
    """Tests for VDJ generative model loading."""
    
    def test_load_human_t_beta_model(self, model_dir):
        """Test loading human T beta generative model."""
        model_folder = os.path.join(model_dir, 'human_T_beta')
        marginals_file = os.path.join(model_folder, 'model_marginals.txt')
        
        generative_model = load_model.GenerativeModelVDJ()
        generative_model.load_and_process_igor_model(marginals_file)
        
        assert generative_model.PV is not None
        assert generative_model.PDJ is not None
        assert generative_model.PinsVD is not None
        assert generative_model.PinsDJ is not None
        assert generative_model.Rvd is not None
        assert generative_model.Rdj is not None
        assert generative_model.PdelV_given_V is not None
        assert generative_model.PdelJ_given_J is not None
        assert generative_model.PdelDldelDr_given_D is not None
        
        # Check probabilities sum to 1
        assert np.allclose(np.sum(generative_model.PV), 1.0, atol=1e-6)
        assert np.allclose(np.sum(generative_model.PDJ), 1.0, atol=1e-6)
        assert np.allclose(np.sum(generative_model.PinsVD), 1.0, atol=1e-6)
        assert np.allclose(np.sum(generative_model.PinsDJ), 1.0, atol=1e-6)


class TestGenerativeModelVJ:
    """Tests for VJ generative model loading."""
    
    def test_load_human_t_alpha_model(self, model_dir):
        """Test loading human T alpha generative model."""
        model_folder = os.path.join(model_dir, 'human_T_alpha')
        marginals_file = os.path.join(model_folder, 'model_marginals.txt')
        
        generative_model = load_model.GenerativeModelVJ()
        generative_model.load_and_process_igor_model(marginals_file)
        
        assert generative_model.PVJ is not None
        assert generative_model.PinsVJ is not None
        assert generative_model.Rvj is not None
        assert generative_model.PdelV_given_V is not None
        assert generative_model.PdelJ_given_J is not None
        
        # Check probabilities sum to 1
        assert np.allclose(np.sum(generative_model.PVJ), 1.0, atol=1e-6)
        assert np.allclose(np.sum(generative_model.PinsVJ), 1.0, atol=1e-6)


class TestModelAttributes:
    """Tests for model attribute consistency."""
    
    def test_vdj_model_entropy(self, model_dir):
        """Test that VDJ model has entropy attributes."""
        model_folder = os.path.join(model_dir, 'human_T_beta')
        marginals_file = os.path.join(model_folder, 'model_marginals.txt')
        
        generative_model = load_model.GenerativeModelVDJ()
        generative_model.load_and_process_igor_model(marginals_file)
        
        assert generative_model.SV is not None
        assert generative_model.SDJ is not None
        assert generative_model.Sscenario is not None
        assert generative_model.Sscenario > 0

