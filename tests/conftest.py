"""Pytest configuration and shared fixtures."""
import os
import pytest
import olga.load_model as load_model
import olga.generation_probability as generation_probability
import olga.sequence_generation as sequence_generation

@pytest.fixture
def model_dir():
    """Get the path to default models directory."""
    return os.path.join(os.path.dirname(__file__), '..', 'olga', 'default_models')


@pytest.fixture
def human_t_beta_model(model_dir):
    """Load human T beta (VDJ) model for testing."""
    model_folder = os.path.join(model_dir, 'human_T_beta')
    params_file = os.path.join(model_folder, 'model_params.txt')
    marginals_file = os.path.join(model_folder, 'model_marginals.txt')
    V_anchor_file = os.path.join(model_folder, 'V_gene_CDR3_anchors.csv')
    J_anchor_file = os.path.join(model_folder, 'J_gene_CDR3_anchors.csv')
    
    genomic_data = load_model.GenomicDataVDJ()
    genomic_data.load_igor_genomic_data(params_file, V_anchor_file, J_anchor_file)
    
    generative_model = load_model.GenerativeModelVDJ()
    generative_model.load_and_process_igor_model(marginals_file)
    
    pgen_model = generation_probability.GenerationProbabilityVDJ(
        generative_model, genomic_data
    )
    seq_gen_model = sequence_generation.SequenceGenerationVDJ(
        generative_model, genomic_data
    )
    
    return {
        'pgen_model': pgen_model,
        'generative_model': generative_model,
        'genomic_data': genomic_data,
        'model_folder': model_folder,
        'seq_gen_model': seq_gen_model
    }


@pytest.fixture
def human_t_alpha_model(model_dir):
    """Load human T alpha (VJ) model for testing."""
    model_folder = os.path.join(model_dir, 'human_T_alpha')
    params_file = os.path.join(model_folder, 'model_params.txt')
    marginals_file = os.path.join(model_folder, 'model_marginals.txt')
    V_anchor_file = os.path.join(model_folder, 'V_gene_CDR3_anchors.csv')
    J_anchor_file = os.path.join(model_folder, 'J_gene_CDR3_anchors.csv')
    
    genomic_data = load_model.GenomicDataVJ()
    genomic_data.load_igor_genomic_data(params_file, V_anchor_file, J_anchor_file)
    
    generative_model = load_model.GenerativeModelVJ()
    generative_model.load_and_process_igor_model(marginals_file)
    
    pgen_model = generation_probability.GenerationProbabilityVJ(
        generative_model, genomic_data
    )
    seq_gen_model = sequence_generation.SequenceGenerationVJ(
        generative_model, genomic_data
    )
    return {
        'pgen_model': pgen_model,
        'generative_model': generative_model,
        'genomic_data': genomic_data,
        'model_folder': model_folder,
        'seq_gen_model': seq_gen_model
    }

