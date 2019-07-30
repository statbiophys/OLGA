
"""
July 30, 2019 
@kmayerb 

To use olga with Python 3, we updated the repo using by the following:

```
2to3 -f all -f idioms -f ws_comma -f set_literal -v -w -n ./
```

- Then code was then manually debugged for additional issues. 

- Manual changes involved fixing integer division in e.g x/y in Python 2.7 to a//b in Python 3.6

- To locally test that the scripts run without compilation, 
default_models, compute_pgen.py, and generate_sequences.py 
were coppied one directory up from olga/

>>> cp ./olga/compute_pgen.py ./compute_pgen.py
>>> cp ./olga/generate_sequences.py ./generate_sequences.py
>>> rsync -r ./olga/default_models/ ./default_models/

Then using Python 3.6.8

Python 3.6.8 |Anaconda, Inc.| (default, Dec 29 2018, 19:04:46) 
[GCC 4.2.1 Compatible Clang 4.0.1 (tags/RELEASE_401/final)] on darwin

>>> python --version
Python 3.6.8 :: Anaconda, Inc.

>>> python compute_pgen.py --humanTRB CASSLGRDGGHEQYF
Pgen of the amino acid sequence CASSLGRDGGHEQYF: 7.253421763151436e-10
Completed pgen computation in: 0.02 seconds.

>>> python generate_sequences.py --humanTRB -n 5
TGTGCCAGCAGCAACTCTGGCGAGCAGTACTTC   CASSNSGEQYF TRBV9   TRBJ2-7
TGTGCCAGCAGCTCCCAGGGACAGATAGACGAGCAGTACTTC  CASSSQGQIDEQYF  TRBV7-9 TRBJ2-7
TGTGCCTGGAAAGGACAGGATACACTGTTTTTT   CAWKGQDTLFF TRBV30  TRBJ2-2
TGTGCCACGGCCGGGGGCGGGGATGGAAACACCATATATTTT  CATAGGGDGNTIYF  TRBV24-1    TRBJ1-3
TGTGCCACCAGTGATCCCCTTTACTGCGGAAATGAAAAACTGTTTTTT    CATSDPLYCGNEKLFF    TRBV24-1    TRBJ1-4

>>> python generate_sequences.py --humanTRB -o example_seqs.tsv -n 1e2
Starting sequence generation... 
Completed generating all 100 sequences in 0.01 seconds.

>>> python compute_pgen.py --humanTRB -i example_seqs.tsv -m 5
python compute_pgen.py --humanTRB -i example_seqs.tsv -m 5
TGTGCCAGCAGTTACGCAAGGCTCCCGGGTGATTTGAACACTGAAGCTTTCTTT  2.5025349112903632e-17  CASSYARLPGDLNTEAFF  7.004131042191115e-13
TGTGCCAGCGAAAGAGGGGTCATGGCCAAAAACATTCAGTACTTC   4.118971361321766e-14   CASERGVMAKNIQYF 8.058582236907881e-12
TGTGCCAGCAGCTGGGGGTACGGAAACACCATATATTTT 1.1214063616323165e-09  CASSWGYGNTIYF   1.8359637976984517e-08
TGCGCCAGCAGCTTGGACGGCGGCCAGGGCTATGGCTACACCTTC   1.124706004150688e-12   CASSLDGGQGYGYTF 1.0116145533668686e-08
TGTGCCAGCAGCCAGCCGCAACTCGATTTACAGGGAGAACATTCACCCCTCCACTTT   6.757203697069712e-22   CASSQPQLDLQGEHSPLHF 5.225492688284462e-17

>>> python olga-compute_pgen --humanTRB -i example_seqs.tsv -o example_pgens.tsv
Successfully read in and formatted 100 sequences and any V or J usages.
...
Completed Pgen computation for 100 sequences: in 4.27 seconds.


Comparison to Pip Installed Version of Olga on Python 2.7.11
------------------------------------------------------------
>>> python --version
Python 2.7.11 :: Continuum Analytics, Inc.

>>>olga-compute_pgen --humanTRB CASSLGRDGGHEQYF
Pgen of the amino acid sequence CASSLGRDGGHEQYF: 7.253421763151438e-10
Completed pgen computation in: 0.03 seconds.

>>>olga-generate_sequences --humanTRB -n 5
TGTGCCTCAAGAAATGAAAAACTGTTTTTT  CASRNEKLFF  TRBV24-1    TRBJ1-4
TGCAGTAAACCAATCTCAGCCGTCCTGTCCACTGAAGCTTTCTTT   CSKPISAVLSTEAFF TRBV20-1    TRBJ1-1
TGTGCCAGCAGCACCAGGACTAGATCTTTTAGGGAGACCCAGTACTTC    CASSTRTRSFRETQYF    TRBV7-9 TRBJ2-5
TGCAGTGAAGATAGGGAGGGGGACAATGAGCAGTTCTTC CSEDREGDNEQFF   TRBV20-1    TRBJ2-1
TGTGCCAGCAGCAACCATGGACCACCATCCGGGACAGGGTCCCTTACCGAGCAGTACTTC    CASSNHGPPSGTGSLTEQYF    TRBV13  TRBJ2-7

Some additional tests:
"""

import unittest
import pandas as pd
import numpy as np

import olga.load_model as load_model
import olga.generation_probability as generation_probability
import olga.generation_probability as pgen
import olga.sequence_generation as seq_gen
from olga.utils import nt2aa, determine_seq_type
import os.path as op

from olga.paths import path_to_olga_default_models


class test_olga_in_python3(unittest.TestCase):

    def test_olga_pgen(self, chain_folder = 'human_T_beta', cdr3 = 'CAWSVAPDRGGYTF', v_b = 'TRBV30*01', v_j ='TRBJ1-2*01'):
        """
        NOT AN ACTUAL UNIT TEST, JUST CODE USED TO TEST IF OLGA WORKED IN PYTHON 3
        """

        params_file_name = op.join(path_to_olga_default_models,
                                   chain_folder,
                                   'model_params.txt')
        marginals_file_name = op.join(path_to_olga_default_models,
                                      chain_folder,
                                      'model_marginals.txt')
        V_anchor_pos_file = op.join(path_to_olga_default_models,
                                    chain_folder,
                                    'V_gene_CDR3_anchors.csv')
        J_anchor_pos_file = op.join(path_to_olga_default_models,
                                    chain_folder,
                                     'J_gene_CDR3_anchors.csv')

        #Load data
        genomic_data = load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        #Load model
        generative_model = load_model.GenerativeModelVDJ()
        generative_model.load_and_process_igor_model(marginals_file_name)

        #Process model/data for pgen computation by instantiating GenerationProbabilityVDJ
        pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)
        #Compute some sequence pgens
        x = pgen_model.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF', 'TRBV30*01', 'TRBJ1-2*01')
        self.assertTrue(np.isclose(x,1.203646865765782e-10, atol =1e-15))


    def test_olga_sgen(self, chain_folder = 'human_T_beta'):
        """
        NOT AN ACTUAL UNIT TEST, JUST CODE USED TO TEST IF OLGA WORKED AT ALL IN PYTHON 3
        """

        params_file_name = op.join(path_to_olga_default_models,
                                   chain_folder,
                                   'model_params.txt')
        marginals_file_name = op.join(path_to_olga_default_models,
                                      chain_folder,
                                      'model_marginals.txt')
        V_anchor_pos_file = op.join(path_to_olga_default_models,
                                    chain_folder,
                                    'V_gene_CDR3_anchors.csv')
        J_anchor_pos_file = op.join(path_to_olga_default_models,
                                    chain_folder,
                                     'J_gene_CDR3_anchors.csv')

        #Load data
        genomic_data = load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        #Load model
        generative_model = load_model.GenerativeModelVDJ()
        generative_model.load_and_process_igor_model(marginals_file_name)

        #Process model/data for sequence generation by instantiating SequenceGenerationVDJ
        seq_gen_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)

        #Generate some random sequences
        x = seq_gen_model.gen_rnd_prod_CDR3()
        self.assertTrue(x)


    def test_olga(self,    chain_folder = 'human_T_beta'):
        """
        NOT AN ACTUAL UNIT TEST, JUST CODE USED TO TEST IF OLGA WORKED AT ALL IN PYTHON 3
        """
        params_file_name = op.join(path_to_olga_default_models,
                                   chain_folder,
                                   'model_params.txt')
        marginals_file_name = op.join(path_to_olga_default_models,
                                      chain_folder,
                                      'model_marginals.txt')
        V_anchor_pos_file = op.join(path_to_olga_default_models,
                                    chain_folder,
                                    'V_gene_CDR3_anchors.csv')
        J_anchor_pos_file = op.join(path_to_olga_default_models,
                                    chain_folder,
                                     'J_gene_CDR3_anchors.csv')

        #Load data
        genomic_data = load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        #Load model
        generative_model = load_model.GenerativeModelVDJ()
        generative_model.load_and_process_igor_model(marginals_file_name)

        #Process model/data for pgen computation by instantiating GenerationProbabilityVDJ
        pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)

        #Compute some sequence pgens
        x = pgen_model.compute_regex_CDR3_template_pgen('CASSAX{0,5}SARPEQFF')
        assert(np.isclose(x,6.846877804096558e-10,atol =1e-15))
        x = pgen_model.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF', 'TRBV30*01', 'TRBJ1-2*01')
        assert(np.isclose(x,1.203646865765782e-10,atol =1e-15))
        x = pgen_model.compute_nt_CDR3_pgen('TGTGCCAGTAGTATAACAACCCAGGGCTTGTACGAGCAGTACTTC')
        assert(np.isclose(x,3.9945642868171824e-14,atol =1e-17))


        #Process model/data for sequence generation by instantiating SequenceGenerationVDJ
        seq_gen_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)

        #Generate some random sequences
        x = seq_gen_model.gen_rnd_prod_CDR3()
        #assert(x == ('TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT', 'CASSEKRQWESGELFF', 27, 8))
        assert(isinstance(x[0], str))# == 'TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT' )
        assert(isinstance(x[1], str))# == 'CASSEKRQWESGELFF')

        #('TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT', 'CASSEKRQWESGELFF', 27, 8)
        x = seq_gen_model.gen_rnd_prod_CDR3()
        assert(isinstance(x[0], str))# == 'TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT' )
        assert(isinstance(x[1], str))# == 'CASSEKRQWESGELFF')

        #assert(x == ('TGTGCCAGCAGTTTAGTGGGAAGGGCGGGGCCCTATGGCTACACCTTC', 'CASSLVGRAGPYGYTF', 14, 1))
        #('TGTGCCAGCAGTTTAGTGGGAAGGGCGGGGCCCTATGGCTACACCTTC', 'CASSLVGRAGPYGYTF', 14, 1)
        x = seq_gen_model.gen_rnd_prod_CDR3()
        assert(isinstance(x[0], str))# == 'TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT' )
        assert(isinstance(x[1], str))# == 'CASSEKRQWESGELFF')



if __name__ == '__main__':
    unittest.main()
