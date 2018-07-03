#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""Command line script to generate sequences.

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

This program will generate a file of Monte Carlo sampling from a specified
genomic model (either the default models that ship with the scripts, or a
custom IGoR model inference). The sequences generated will have NO ERRORS.
Input arguments are specified by an immediately preceding flag. Capitalization
is irrelevant for the flags and limited option string arguments (but can matter
for arbitrary string arguments).

The following are required inputs:

1) File name for output file: PATH/TO/OUTFILE
THIS IS ASSUMED TO BE THE FIRST ARGUMENT (see example calls)

2) Number of sequences to generate.
-num_seqs or -n : int or float
    If a float number is given, it will be rounded. However, it is convenient
    to be able to specify numbers in scientific notation, e.g. 1e7 is easier than
    10000000.

3) Generative model used to define the generation probability  of a sequence.
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

#Generate 1,000 seqs from the human TCRB default model, delimiter is tab
./generate_synthetic_sequences.py example_seqs.tsv --humanTCRB -n 1e3

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
  -n N, --num_seqs=N    specify the number of sequences to generate.
  --seed=SEED           set seed for pseudorandom number generator. Default is
                        to not set a seed.
  --seqs_per_time_update=SEQS_PER_TIME_UPDATE
                        specify the number of sequences between time updates.
                        Default is 1e5
  --conserved_J_residues=CONSERVED_J_RESIDUES
                        specify conserved J residues. Default is 'FW'.
  --time_updates_off    turn time updates off.
  --seq_type=SEQ_TYPE   declare sequence type for output sequences. Choices:
                        'all' [default], 'ntseq', 'nucleotide', 'aaseq',
                        'amino_acid'
  --record_genes_off    turn off recording V and J gene info.
  -d DELIMITER, --delimiter=DELIMITER
                        declare delimiter choice. Default is tab for .tsv
                        output files, comma for .csv files, and tab for all
                        others. Choices: 'tab', 'space', ',', ';', ':'
  --raw_delimiter=DELIMITER
                        declare delimiter choice as a raw string.


Note about conserved_J_residues:

This specifies a string which must be composed ONLY of amino acids
(i.e. only ACDEFGHIKLMNPQRSTVWY*). The amino acids in that string will
determine functionality of a sequence. Please note that the J genes are
ALREADY ANCHORED at a given residue, thus this string should almost
certainly only include phenylalanine (F) and/or tryptophan (W). If amino
acids are used to define this 'conserved residue' in the string here, but
the J genes are still anchored at a 'F' or 'W' (as the default genomic
files are), this could result in no productive sequences being generated.
Unless the anchor positions are changed, LEAVE THE DEFAULT. The default
string is 'FW'.

-------------------------------------------------------------------------------
Example calls with options

#Output a .csv file with delimiter = comma.
./generate_synthetic_sequences.py example_seqs.csv --VDJ_model_folder models/human_T_beta/ -n 1e3

#Output the amino acid CDR3 seqs and not the nucleotide sequences
./generate_synthetic_sequences.py example_seqs.tsv --humanTCRB -n 1e3 --seq_type amino_acid

#Don't record the V and J genes/alleles used to generate the sequences
./generate_synthetic_sequences.py example_seqs.tsv --humanTCRB -n 1e3 --record_genes_off

#Observe the time updates --- default is 1e5:
./generate_synthetic_sequences.py example_seqs.tsv --humanTCRB -n 1e3 --seqs_per_time_update 100

#Don't give time updates
./generate_synthetic_sequences.py example_seqs.tsv --humanTCRB -n 1e3 --seqs_per_time_update 100 --time_updates_off

#Set the pseudo-random number generator seed to 100 -- should be repeatable.
./generate_synthetic_sequences.py example_seqs.tsv --humanTCRB -n 1e3 --seed 100

#Specify that the conserved J residue is an F
./generate_synthetic_sequences.py example_seqs.tsv --humanTCRB -n 1e3 --conserved_J_residues 'F'

#Change the delimiter of the output file to space ' '.
./generate_synthetic_sequences.py example_seqs.tsv --humanTCRB -n 1e3 --delimiter space

@author: zacharysethna

"""

#Function assumes that it is in the same directory that the folder app/ is
#in (which should contain all the modules imported).

import os
import sys


import olga.load_model as load_model
import olga.sequence_generation as sequence_generation
from optparse import OptionParser
import time
import numpy as np

def main():
    """ Generate sequences."""

    parser = OptionParser(conflict_handler="resolve")

    parser.add_option('--humanTCRA', '--human_T_alpha', action='store_true', dest='humanTCRA', default=False, help='use default human TCRA model (T cell alpha chain)')
    parser.add_option('--humanTCRB', '--human_T_beta', action='store_true', dest='humanTCRB', default=False, help='use default human TCRB model (T cell beta chain)')
    parser.add_option('--mouseTCRB', '--mouse_T_beta', action='store_true', dest='mouseTCRB', default=False, help='use default mouse TCRB model (T cell beta chain)')
    parser.add_option('--humanIGH', '--human_B_heavy', action='store_true', dest='humanIGH', default=False, help='use default human IGH model (B cell heavy chain)')
    parser.add_option('--VDJ_model_folder', dest='vdj_model_folder', metavar='PATH/TO/FOLDER/', help='specify PATH/TO/FOLDER/ for a custom VDJ generative model')
    parser.add_option('--VJ_model_folder', dest='vj_model_folder', metavar='PATH/TO/FOLDER/', help='specify PATH/TO/FOLDER/ for a custom VJ generative model')
    parser.add_option('-n', '--num_seqs', type='float', metavar='N', default = 0, dest='num_seqs_to_generate', help='specify the number of sequences to generate.')
    parser.add_option('--seed', type='int', dest='seed', help='set seed for pseudorandom number generator. Default is to not set a seed.')
    parser.add_option('--seqs_per_time_update', type='float', default = 100000, dest='seqs_per_time_update', help='specify the number of sequences between time updates. Default is 1e5')
    parser.add_option('--conserved_J_residues', type='string', default = 'FW', dest='conserved_J_residues', help="specify conserved J residues. Default is 'FW'.")
    parser.add_option('--time_updates_off', action='store_false', dest='time_updates', default=True, help='turn time updates off.')
    parser.add_option('--seq_type', type='choice', default = 'all', dest='seq_type',  choices=['all', 'ntseq', 'nucleotide', 'aaseq', 'amino_acid'], help="declare sequence type for output sequences. Choices: 'all' [default], 'ntseq', 'nucleotide', 'aaseq', 'amino_acid'")
    parser.add_option('--record_genes_off', action='store_false', dest="record_genes", default=True, help='turn off recording V and J gene info.')
    parser.add_option('-d', '--delimiter', type='choice', dest='delimiter',  choices=['tab', 'space', ',', ';', ':'], help="declare delimiter choice. Default is tab for .tsv output files, comma for .csv files, and tab for all others. Choices: 'tab', 'space', ',', ';', ':'")
    parser.add_option('--raw_delimiter', type='str', dest='delimiter', help="declare delimiter choice as a raw string.")


    (options, args) = parser.parse_args()

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

    #OUTFILE IS THE FIRST ARGUMENT
    try:
        outfile_name = args[0]
    except IndexError:
        print 'Need to specify outfile as the first argument!'
        print 'Exiting...'
        return -1

    if len(outfile_name.strip()) == 0:
        print 'Need to specify outfile as the first argument!'
        print 'Exiting...'
        return -1

    if os.path.isfile(outfile_name):
        if not raw_input(outfile_name + ' already exists. Overwrite (y/n)? ').strip().lower() in ['y', 'yes']:
            print 'Exiting...'
            return -1

    #Parse arguments

    num_seqs_to_generate = int(options.num_seqs_to_generate)

    if num_seqs_to_generate <= 0:
        print 'Need to specify num_seqs (number of sequences to generate).'
        print 'Exiting...'
        return -1


    #Parse default delimiter
    delimiter = options.delimiter
    if delimiter is None:
        delimiter = '\t'
        if outfile_name.endswith('.tsv'):
            delimiter = '\t'
        elif outfile_name.endswith('.csv'):
            delimiter = ','
    else:
        try:
            delimiter = {'tab': '\t', 'space': ' ', ',': ',', ';': ';', ':': ':'}[delimiter]
        except KeyError:
            pass #Other raw string.

    #Optional flags
    seq_type = {'all': 'all', 'ntseq': 'ntseq', 'nucleotide': 'ntseq', 'aaseq': 'aaseq', 'amino_acid': 'aaseq'}[options.seq_type]
    record_genes = options.record_genes
    seqs_per_time_update = int(options.seqs_per_time_update)
    time_updates = options.time_updates
    conserved_J_residues = options.conserved_J_residues


    if options.seed is not None:
        np.random.seed(options.seed)


    #VDJ recomb case --- used for TCRB and IGH
    if recomb_type == 'VDJ':
        genomic_data = load_model.GenomicDataVDJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = load_model.GenerativeModelVDJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        seq_gen = sequence_generation.SequenceGenerationVDJ(generative_model, genomic_data)
    #VJ recomb case --- used for TCRA and light chain
    elif recomb_type == 'VJ':
        genomic_data = load_model.GenomicDataVJ()
        genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
        generative_model = load_model.GenerativeModelVJ()
        generative_model.load_and_process_igor_model(marginals_file_name)
        seq_gen = sequence_generation.SequenceGenerationVJ(generative_model, genomic_data)


    V_gene_names = [V[0].split('*')[0] for V in genomic_data.genV]
    J_gene_names = [J[0].split('*')[0] for J in genomic_data.genJ]


    outfile = open(outfile_name, 'w')

    print 'Starting sequence generation... '
    start_time = time.time()
    for i in range(num_seqs_to_generate):
        ntseq, aaseq, V_in, J_in = seq_gen.gen_rnd_prod_CDR3(conserved_J_residues)
        if seq_type == 'all': #default, include both ntseq and aaseq
            current_line_out = ntseq + delimiter + aaseq
        elif seq_type == 'ntseq': #only record ntseq
            current_line_out = ntseq
        elif seq_type == 'aaseq': #only record aaseq
            current_line_out = aaseq

        if record_genes:
            current_line_out += delimiter + V_gene_names[V_in] + delimiter + J_gene_names[J_in]
        outfile.write(current_line_out + '\n')

        if (i+1)%seqs_per_time_update == 0 and time_updates:
            c_time = time.time() - start_time
            eta = ((num_seqs_to_generate - (i+1))/float(i+1))*c_time
            if c_time > 86400: #more than a day
                c_time_str = '%d days, %d hours, %d minutes, and %.2f seconds.'%(int(c_time)/86400, (int(c_time)/3600)%24, (int(c_time)/60)%60, c_time%60)
            elif c_time > 3600: #more than an hr
                c_time_str = '%d hours, %d minutes, and %.2f seconds.'%((int(c_time)/3600)%24, (int(c_time)/60)%60, c_time%60)
            elif c_time > 60: #more than a min
                c_time_str = '%d minutes and %.2f seconds.'%((int(c_time)/60)%60, c_time%60)
            else:
                c_time_str = '%.2f seconds.'%(c_time)

            if eta > 86400: #more than a day
                eta_str = '%d days, %d hours, %d minutes, and %.2f seconds.'%(int(eta)/86400, (int(eta)/3600)%24, (int(eta)/60)%60, eta%60)
            elif eta > 3600: #more than an hr
                eta_str = '%d hours, %d minutes, and %.2f seconds.'%((int(eta)/3600)%24, (int(eta)/60)%60, eta%60)
            elif eta > 60: #more than a min
                eta_str = '%d minutes and %.2f seconds.'%((int(eta)/60)%60, eta%60)
            else:
                eta_str = '%.2f seconds.'%(eta)

            print '%d sequences generated in %s Estimated time remaining: %s'%(i+1, c_time_str, eta_str)

    c_time = time.time() - start_time
    if c_time > 86400: #more than a day
        c_time_str = '%d days, %d hours, %d minutes, and %.2f seconds.'%(int(c_time)/86400, (int(c_time)/3600)%24, (int(c_time)/60)%60, c_time%60)
    elif c_time > 3600: #more than an hr
        c_time_str = '%d hours, %d minutes, and %.2f seconds.'%((int(c_time)/3600)%24, (int(c_time)/60)%60, c_time%60)
    elif c_time > 60: #more than a min
        c_time_str = '%d minutes and %.2f seconds.'%((int(c_time)/60)%60, c_time%60)
    else:
        c_time_str = '%.2f seconds.'%(c_time)
    print 'Completed generated all %d sequences in %s'%(num_seqs_to_generate, c_time_str)
    outfile.close()

if __name__ == '__main__': main()
