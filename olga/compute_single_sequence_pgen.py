#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Command line script for Pgen computation of a single sequence.

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


DON'T USE THIS FUNCTION TO LOOP OVER A LOT OF SEQUENCES. This function loads
up and processes a full model to compute Pgen of a single sequence, and thus
pays a nontrivial overhead time cost. If many sequences need to be run it is
better to use run_pgen.py or to write your own wrapper.

The sequence which is read in must be TRIMMED TO ONLY THE CDR3 region as
defined by the V and J anchor files (default is to INCLUDE the conserved
residues of the C in the V region and the F/W in the J region).

It is also possible to condition the Pgen computation on V and J identity by
specifying the V or J usages. See documentation below for an explanation of the
syntax and examples of calls.

Required input:

1) Sequence
THIS IS ASSUMED TO BE THE FIRST ARGUMENT (see example calls)

2) Generative model used to define the generation probability  of a sequence.
Flags for default models:

--humanTCRA or --human_T_alpha (default Human T cell alpha chain model)
--humanTCRB or --human_T_beta (default Human T cell beta chain model)
--mouseTCRB or --mouse_T_beta (default Mouse T cell beta chain model)
--humanIGH or --human_B_heavy (default Human B cell heavy chain model)

In order to use these default model flags, the program must be executed from
the same directory where models/ is. For example, models/human_T_beta/ should
be an okay pathing to the humanTCRB model folder.

To specify a custom model folder use:

--VJ_model_folder (a generative model of VJ recombination, e.g. T alpha chain)
--VDJ_model_folder (a generative model of VDJ recombination, e.g. T beta chain)

Note, if specifying a custom model folder for either a VJ recombination model
(e.g. an alpha or light chain model) or a VDJ recombination model
(e.g. a beta or heavy chain model), the folder must contain the following files
with the exact naming convention:

model_params.txt (IGoR inference param file)
model_marginals.txt (IGoR inference marginal file)
V_gene_CDR3_anchors.csv (V residue anchor and functionality file)
J_gene_CDR3_anchors.csv (J residue anchor and functionality file)
-------------------------------------------------------------------------------
Example call (example calls are formatted as executed functions instead of the
console script entry point. The arguments are identical in either case.)

```
$ ./compute_single_sequence_pgen.py CARQGALYEQYF --humanTCRB

------------------------------------------------------------------------------------------
Pgen of the amino acid sequence CARQGALYEQYF: 1.73865590343e-08
------------------------------------------------------------------------------------------
```


-------------------------------------------------------------------------------
Options:
  -h, --help            show this help message and exit
  --humanTCRA, --human_T_alpha
                        use default human TCRA model (T cell alpha chain)
  --humanTCRB, --human_T_beta
                        use default human TCRB model (T cell beta chain)
  --mouseTCRB, --mouse_T_beta
                        use default mouse TCRB model (T cell beta chain)
  --humanIGH, --human_B_heavy
                        use default human IGH model (B cell heavy chain)
  --VDJ_model_folder=PATH/TO/FOLDER/
                        specify PATH/TO/FOLDER/ for a custom VDJ generative
                        model
  --VJ_model_folder=PATH/TO/FOLDER/
                        specify PATH/TO/FOLDER/ for a custom VJ generative
                        model
  -v V_MASK, --v_mask=V_MASK
                        specify V usage to condition Pgen on
  -j J_MASK, --j_mask=J_MASK
                        specify J usage to condition Pgen on
  -a PATH/TO/FILE, --alphabet_filename=PATH/TO/FILE
                        specify PATH/TO/FILE defining a custom 'amino acid'
                        alphabet
  -t, --time            print time to compute Pgen.
  --seq_type=SEQ_TYPE   declare sequence type. Infers seq_type by default.
                        Choices: 'ntseq', 'nucleotide', 'aaseq', 'amino_acid',
                        'regex', 'regular_expression'
  --print_warnings_off  turn Pgen print warnings off

Example call with optional flags:

```
$ ./compute_single_sequence_pgen.py C[AI]LXXGSNY[KQ]L[TI][FW] --seq_type regex --humanTCRA -v 23/DV6,21 -j 53,33 -t
output:
------------------------------------------------------------------------------------------
Pgen of the regular expression sequence C[AI]LXXGSNY[KQ]L[TI][FW]: 1.38162611621e-06

(Conditioned on the V and J gene/allele usages: ['23/DV6', '21'], ['53', '33'])
------------------------------------------------------------------------------------------
Completed pgen computation in: 0.01 seconds.

```

@author: zacharysethna

"""

#Function assumes that it is in the same directory that the folder app/ is
#in (which should contain all the modules imported).

import os
import sys
import time

import olga.load_model as load_model
import olga.generation_probability as generation_probability
from olga.utils import nt2aa
from optparse import OptionParser

def main():
    """Compute Pgen of a single sequence."""


    parser = OptionParser(conflict_handler="resolve")

    parser.add_option('--humanTCRA', '--human_T_alpha', action='store_true', dest='humanTCRA', default=False, help='use default human TCRA model (T cell alpha chain)')
    parser.add_option('--humanTCRB', '--human_T_beta', action='store_true', dest='humanTCRB', default=False, help='use default human TCRB model (T cell beta chain)')
    parser.add_option('--mouseTCRB', '--mouse_T_beta', action='store_true', dest='mouseTCRB', default=False, help='use default mouse TCRB model (T cell beta chain)')
    parser.add_option('--humanIGH', '--human_B_heavy', action='store_true', dest='humanIGH', default=False, help='use default human IGH model (B cell heavy chain)')
    parser.add_option('--VDJ_model_folder', dest='vdj_model_folder', metavar='PATH/TO/FOLDER/', help='specify PATH/TO/FOLDER/ for a custom VDJ generative model')
    parser.add_option('--VJ_model_folder', dest='vj_model_folder', metavar='PATH/TO/FOLDER/', help='specify PATH/TO/FOLDER/ for a custom VJ generative model')
    parser.add_option('-v', '--v_mask', type='string', dest='V_mask', help='specify V usage to condition Pgen on')
    parser.add_option('-j', '--j_mask', type='string', dest='J_mask', help='specify J usage to condition Pgen on')
    parser.add_option('-a', '--alphabet_filename', dest='alphabet_filename', metavar='PATH/TO/FILE', help="specify PATH/TO/FILE defining a custom 'amino acid' alphabet")
    parser.add_option('-t', '--time', action='store_true', dest='print_c_time', default=False, help='print time to compute Pgen.')
    parser.add_option('--seq_type', type='choice', dest='seq_type',  choices=['ntseq', 'nucleotide', 'aaseq', 'amino_acid', 'regex', 'regular_expression'], help="declare sequence type. Infers seq_type by default. Choices: 'ntseq', 'nucleotide', 'aaseq', 'amino_acid', 'regex', 'regular_expression'")
    parser.add_option('--print_warnings_off', action='store_false', dest="print_warnings", default=True, help='turn Pgen print warnings off')



    (options, args) = parser.parse_args()

    #SEQUENCE IS THE FIRST ARGUMENT
    try:
        seq = args[0]
    except IndexError:
        print 'Need to specify the CDR3 sequence as the first argument!'
        print 'Exiting...'
        return -1

    if len(seq.strip()) == 0:
        print 'Need to specify the CDR3 sequence as the first argument!'
        print 'Exiting...'
        return -1

    main_folder = os.path.dirname(__file__)

    default_models = {}
    default_models['humanTCRA'] = [os.path.join(main_folder, 'default_models', 'human_T_alpha'),  'VJ']
    default_models['humanTCRB'] = [os.path.join(main_folder, 'default_models', 'human_T_beta'), 'VDJ']
    default_models['mouseTCRB'] = [os.path.join(main_folder, 'default_models', 'mouse_T_beta'), 'VDJ']
    default_models['humanIGH'] = [os.path.join(main_folder, 'default_models', 'human_B_heavy'), 'VDJ']

    num_models_specified = sum([1 for x in default_models.keys() + ['vj_model_folder', 'vdj_model_folder'] if getattr(options, x)])

    if num_models_specified == 1: #exactly one model specified
        try:
            d_model = [x for x in default_models.keys() if getattr(options, x)][0]
            model_folder = default_models[d_model][0]
            recomb_type = default_models[d_model][1]
        except IndexError:
            if options.vdj_model_folder: #custom VDJ model specified
                model_folder = options.vdj_model_folder
                recomb_type = 'VDJ'
            elif options.vj_model_folder: #custom VJ model specified
                model_folder = options.vj_model_folder
                recomb_type = 'VJ'
    elif num_models_specified == 0:
        print 'Need to indicate generative model.'
        print 'Exiting...'
        return -1
    elif num_models_specified > 1:
        print 'Only specify one model'
        print 'Exiting...'
        return -1

    #Check that all model and genomic files exist in the indicated model folder
    if not os.path.isdir(model_folder):
        print 'Check pathing... cannot find the model folder: ' + model_folder
        print 'Exiting...'
        return -1

    params_file_name = os.path.join(model_folder,'model_params.txt')
    marginals_file_name = os.path.join(model_folder,'model_marginals.txt')
    V_anchor_pos_file = os.path.join(model_folder,'V_gene_CDR3_anchors.csv')
    J_anchor_pos_file = os.path.join(model_folder,'J_gene_CDR3_anchors.csv')

    for x in [params_file_name, marginals_file_name, V_anchor_pos_file, J_anchor_pos_file]:
        if not os.path.isfile(x):
            print 'Cannot find: ' + x
            print 'Please check the files (and naming conventions) in the model folder ' + model_folder
            print 'Exiting...'
            return -1


    #Optional flags
    alphabet_filename = options.alphabet_filename #used if a custom alphabet is to be specified
    seq_type = options.seq_type
    print_warnings = options.print_warnings
    print_c_time = options.print_c_time

    #Format V and J masks
    try:
        V_mask = options.V_mask.split(',')
    except AttributeError:
        V_mask = options.V_mask #Default is None, i.e. not conditioning on V identity

    try:
        J_mask = options.J_mask.split(',')
    except AttributeError:
         J_mask = options.J_mask #Default is None, i.e. not conditioning on J identity

    if alphabet_filename is not None:
        if not os.path.isfile(alphabet_filename):
            print 'Cannot find extended alphabet file: ' + alphabet_filename
            print 'Exiting...'
            return -1

    #Default --- infer type of sequence.
    if seq_type is None:
        if all([x in 'ACGTacgt' for x in seq]):
            seq_type = 'ntseq'
        elif any([x in '[]{}0123456789,' for x in seq]):
            seq_type = 'regex'
        else:
            seq_type = 'aaseq'
    else:
        #Must be one of these keys
        seq_type = {'ntseq': 'ntseq', 'nucleotide': 'ntseq', 'aaseq': 'aaseq', 'amino_acid': 'aaseq', 'regex': 'regex', 'regular_expression': 'regex'}[seq_type]

        #Vet seq_type
        if seq_type in ['aaseq', 'regex']:
            if all([x in 'ACGTacgt' for x in seq]):
                print 'It looks like the sequence read in is a nucleotide sequence.'
                if raw_input('Change read in seq_type to nucleotide? (y/n)? ') in ['y', 'yes']:
                    seq_type = 'ntseq'
                else:
                    print 'Okay... you are the boss... chances are the Pgen will be 0.'

        if seq_type in ['ntseq', 'aaseq']:
            if any([x in '[]{}0123456789,' for x in seq]):
                print 'It looks like sequence read in includes symbols of a regular expression sequence.'
                if raw_input('Change read in seq_type to regular expression? (y/n)? ') in ['y', 'yes']:
                    seq_type = 'regex'
                else:
                    print 'Okay... you are the boss... cannot read the sequence in that case.'


    #Load up model based on recomb_type
    #VDJ recomb case --- used for TCRB and IGH
    if recomb_type == 'VDJ':
        genomic_data = load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = load_model.GenerativeModelVDJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        pgen_model = generation_probability.GenerationProbabilityVDJ(generative_model, genomic_data, alphabet_filename)
    #VJ recomb case --- used for TCRA and light chain
    elif recomb_type == 'VJ':
        genomic_data = load_model.GenomicDataVJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = load_model.GenerativeModelVJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        pgen_model = generation_probability.GenerationProbabilityVJ(generative_model, genomic_data, alphabet_filename)

    if seq_type == 'ntseq' and not all([x in 'ACGTacgt' for x in seq]):
        if all([x in pgen_model.codons_dict.keys() for x in seq]):
            print 'It looks like sequence read in is an amino acid sequence.'
            if raw_input('Change read in seq_type to amino acid? (y/n)? ') in ['y', 'yes']:
                seq_type = 'aaseq'
            else:
                print 'Okay... you are the boss... cannot read the sequence in that case.'

    print ''
    if seq_type == 'regex' and all([x in pgen_model.codons_dict.keys() for x in seq]):
        seq_type = 'aaseq' #amino acid sequence provided... just change the seq_type for a better printed output

    start_time = time.time()
    if seq_type == 'aaseq':
        pgen = pgen_model.compute_aa_CDR3_pgen(seq, V_mask, J_mask, print_warnings)
    elif seq_type == 'regex':
        pgen = pgen_model.compute_regex_CDR3_template_pgen(seq, V_mask, J_mask, print_warnings)
    elif seq_type == 'ntseq':
        pgen_nt = pgen_model.compute_nt_CDR3_pgen(seq, V_mask, J_mask, print_warnings)
        pgen_aa = pgen_model.compute_aa_CDR3_pgen(nt2aa(seq), V_mask, J_mask, print_warnings)

    c_time = time.time() - start_time
    if c_time > 86400: #more than a day
        c_time_str = '%d days, %d hours, %d minutes, and %.2f seconds.'%(int(c_time)/86400, (int(c_time)/3600)%24, (int(c_time)/60)%60, c_time%60)
    elif c_time > 3600: #more than an hr
        c_time_str = '%d hours, %d minutes, and %.2f seconds.'%((int(c_time)/3600)%24, (int(c_time)/60)%60, c_time%60)
    elif c_time > 60: #more than a min
        c_time_str = '%d minutes and %.2f seconds.'%((int(c_time)/60)%60, c_time%60)
    else:
        c_time_str = '%.2f seconds.'%(c_time)

    print '-'*90
    if seq_type == 'aaseq':
        print 'Pgen of the amino acid sequence ' + seq + ': ' + str(pgen)
    elif seq_type == 'regex':
        print 'Pgen of the regular expression sequence ' + seq + ': ' + str(pgen)
    elif seq_type == 'ntseq':
        print 'Pgen of the nucleotide sequence ' + seq + ': ' + str(pgen_nt)
        print 'Pgen of the amino acid sequence nt2aa(' + seq + ') = ' + nt2aa(seq) + ': ' + str(pgen_aa)

    if V_mask is None:
        V_mask = []
    if J_mask is None:
        J_mask = []

    V_mask = [v for v in V_mask if v in pgen_model.V_mask_mapping.keys()]
    J_mask = [j for j in J_mask if j in pgen_model.J_mask_mapping.keys()]

    if not len(V_mask) == 0 and not len(J_mask) == 0:
        print ''
        print '(Conditioned on the V and J gene/allele usages: ' + str(V_mask) + ', ' + str(J_mask) + ')'
    elif not len(V_mask) == 0:
        print ''
        print '(Conditioned on the V gene/allele usages: ' + str(V_mask) + ')'
    elif not len(J_mask) == 0:
        print ''
        print '(Conditioned on the J gene/allele usages: ' + str(J_mask) + ')'

    print '-'*90
    if print_c_time:
        print 'Completed pgen computation in: ' + c_time_str

    print ''

if __name__ == '__main__': main()
