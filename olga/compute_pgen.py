#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Command line script to compute Pgens of CDR3 sequences.

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

This script will compute  the generation probabilities (Pgens) of sequences, as
defined by a specified generative model. The sequences must be TRIMMED TO ONLY
THE CDR3 region as defined by the V and J anchor files (default is to INCLUDE
the conserved residues of the C in the V region and the F/W in the J region).

Each sequence will be determined to be one of:
1) 'Amino acid' sequence including any ambiguous amino acid symbols (see
Options for alphabet_file to specify custom symbols)
2) Nucleotide sequence (in-frame)
3) Regular expression template for 'amino acid sequences'

There are four default generative models that ship with OLGA and can be
specified with a flag:
--humanTRA (Human T cell alpha chain VJ model)
--humanTRB (Human T cell beta chain VDJ model)
--mouseTRB (Mouse T cell beta chain VDJ model)
--humanIGH (Human B cell heavy chain VDJ model)

A custom model can also be specified (see below for details).

This script can read in sequences and output pgens in one of three modes:

1) Pass CDR3 sequences as arguments, output is printed to stdout.
2) Read in CDR3 sequences from a file, output is printed to stdout.
3) Read in CDR3 sequences from a file, output to a file, dynamic display with
time updates printed to stdout.


Mode 1):

This mode is provided as a way to quickly compute pgen of just a single or
couple of sequences. It is not advisable to use this mode to set up a loop
over a lot of sequences as each call of the script demands the overhead of
processing a model. To compute the pgens of many sequences, it is suggested to
read the sequences in from a file and use either mode 2 or 3.

It is also possible to condition the Pgen computation on V and J identity by
specifying the V or J usages as a mask. However, note that these V/J masks will
be applied to ALL of the sequences provided as arguments. Read Options on
v_mask and j_mask for more info.

--------------------------------------------------------------------------------
Example calls:

$ olga-compute_pgen --humanTRB CASSTGQANYGYTF
Pgen of the amino acid sequence CASSTGQANYGYTF: 5.26507446955e-08

$ olga-compute_pgen --humanTRB TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT
Pgen of the nucleotide sequence TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT: 1.31873701121e-17
Pgen of the amino acid sequence nt2aa(TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT) = CASSDAQGRNRGTEAFF: 4.70599549953e-13

#For a regular expression sequence, backslashes may be needed to specify the
#characters {} if seq read in as args
$ olga-compute_pgen --humanTRB CASSTGX\{1,5\}QAN[YA]GYTF
Pgen of the regular expression sequence CASSTGX{1,5}QAN[YA]GYTF: 7.588241802e-08

$ olga-compute_pgen --humanTRB CASSTGQANYGYTF CASSTGX\{1,5\}QAN[YA]GYTF TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT
Pgen of the amino acid sequence CASSTGQANYGYTF: 5.26507446955e-08
Pgen of the regular expression sequence CASSTGX{1,5}QAN[YA]GYTF: 7.588241802e-08
Pgen of the nucleotide sequence TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT: 1.31873701121e-17
Pgen of the amino acid sequence nt2aa(TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT) = CASSDAQGRNRGTEAFF: 4.70599549953e-13

$ olga-compute_pgen --humanTRB CASSTGQANYGYTF --v_mask TRBV14 --j_mask TRBJ1-2
Pgen of the amino acid sequence CASSTGQANYGYTF: 5.5513032863e-10

--------------------------------------------------------------------------------

Modes 2/3):

These read in sequences from a file. The script has only minimal file parsing
built in, so reading in sequences from a file requires the file to be structured
with delimiter spaced values (i.e. the  data is organized in columns separated
by delimiter like a .tsv or .csv file). Read Options on delimiter for more info.

It is not recommended to read in regular expression sequences from a file. These
sequences require enumerating out the amino acid sequences which correspond to
them and computing pgen for each of them individually -- this can require a
large time cost. Instead consider defining a custom 'amino acid' alphabet to
define the symbols used in the regular expressions if possible. Furthermore,
BE CAREFUL if reading in from a .csv file -- if commas are used in a regex
sequence and comma is used as the delimiter of the .csv file, the sequence will
not be read in properly.

If nucleotide sequences are to be read in it is possible to specify if the
output should be the nucleotide sequence Pgen and/or the translated amino acid
sequence Pgen (the default is to compute and output both). See Options.

To read in sequences, the index of column of CDR3 sequences is needed. The
default is to assume that the sequences to be read in are in the first column
(index 0), meaning that a text file with only a sequence on each line will be
read in okay by default. Read Options on seq_in for more info.

It is also possible to condition the Pgen computation on V and J identity by
specifying what index the column that V and J masks are stored for each line.

Mode 2 does not have a specified output file and so will print the sequences
and their pgens to stdout.

Mode 3 does have a specified output file. By default in this mode there is a
running display of the last few sequences/pgens written to the output file as
well as time elapsed, current rate of computation, and estimated time remaining.
This display can be disabled (see Options).

As it is rare for datasets to be >> 1e4, parallelization is not built in.
However, there are options to skip N lines of the file and to load at most M
sequences so, if wanted, one could build a parallelized wrapper around this
script (though it would be recommended to instead just import the modules and
build from there).

--------------------------------------------------------------------------------
Example calls (assumes a file example_seqs.tsv with the line structure
ntseq   aaseq   V_mask  J_mask):

#Reads in the ntseqs, prints the ntseq, aaseq and their pgens to stdout
$ olga-compute_pgen -i example_seqs.tsv --humanTRB

#Reads in the ntseqs, writes the ntseq, aaseq and their pgens to example_pgens.tsv
$ olga-compute_pgen -i example_seqs.tsv --humanTRB -o example_pgens.tsv

#Specifies the V/J mask indices
$ olga-compute_pgen -i example_seqs.tsv --humanTRB -o example_pgens.tsv --v_in 2 --j_in 3

#Reads in the aaseq column
$ olga-compute_pgen -i example_seqs.tsv --humanTRB -o example_pgens.tsv --seq_in 1

#Only runs the first 100 sequences
$ olga-compute_pgen -i example_seqs.tsv --humanTRB -o example_pgens.tsv --seq_in 1 -m 100
--------------------------------------------------------------------------------

To specify a custom model folder use:
--set_custom_model_VJ (generative model of VJ recombination, e.g. T alpha chain)
--set_custom_model_VDJ (generative model of VDJ recombination, e.g. T beta chain)

Note, if specifying a custom model folder for either a VJ recombination model
(e.g. an alpha or light chain model) or a VDJ recombination model
(e.g. a beta or heavy chain model), the folder must contain the following files
with the exact naming convention:

model_params.txt (IGoR inference param file)
model_marginals.txt (IGoR inference marginal file)
V_gene_CDR3_anchors.csv (V residue anchor and functionality file)
J_gene_CDR3_anchors.csv (J residue anchor and functionality file)

--------------------------------------------------------------------------------
Options:
  -h, --help            show this help message and exit
  --humanTRA, --human_T_alpha
                        use default human TRA model (T cell alpha chain)
  --humanTRB, --human_T_beta
                        use default human TRB model (T cell beta chain)
  --mouseTRB, --mouse_T_beta
                        use default mouse TRB model (T cell beta chain)
  --humanIGH, --human_B_heavy
                        use default human IGH model (B cell heavy chain)
  --humanIGK
                        use default human IGK model
  --set_custom_model_VDJ=PATH/TO/FOLDER/
                        specify PATH/TO/FOLDER/ for a custom VDJ generative
                        model
  --set_custom_model_VJ=PATH/TO/FOLDER/
                        specify PATH/TO/FOLDER/ for a custom VJ generative
                        model
  -i PATH/TO/FILE, --infile=PATH/TO/FILE
                        read in CDR3 sequences (and optionally V/J masks) from
                        PATH/TO/FILE
  -o PATH/TO/FILE, --outfile=PATH/TO/FILE
                        write CDR3 sequences and pgens to PATH/TO/FILE
  --seq_in=INDEX, --seq_index=INDEX
                        specifies sequences to be read in are in column INDEX.
                        Default is index 0 (the first column).
  --v_in=INDEX, --v_mask_index=INDEX
                        specifies V_masks are found in column INDEX in the
                        input file. Default is no V mask.
  --j_in=INDEX, --j_mask_index=INDEX
                        specifies J_masks are found in column INDEX in the
                        input file. Default is no J mask.
  --v_mask=V_MASK       specify V usage to condition Pgen on for seqs read in
                        as arguments.
  --j_mask=J_MASK       specify J usage to condition Pgen on for seqs read in
                        as arguments.
  -m N, --max_number_of_seqs=N
                        compute Pgens for at most N sequences.
  --lines_to_skip=N     skip the first N lines of the file. Default is 0.
  -a PATH/TO/FILE, --alphabet_filename=PATH/TO/FILE
                        specify PATH/TO/FILE defining a custom 'amino acid'
                        alphabet. Default is no custom alphabet.
  --seq_type_out=SEQ_TYPE
                        if read in sequences are ntseqs, declare what type of
                        sequence to compute pgen for. Default is all. Choices:
                        'all', 'ntseq', 'nucleotide', 'aaseq', 'amino_acid'
  --skip_off, --skip_empty_off
                        stop skipping empty or blank sequences/lines (if for
                        example you want to keep line index fidelity between
                        the infile and outfile).
  --display_off         turn the sequence display off (only applies in safe
                        mode). Default is on.
  --num_lines_for_display=N
                        N lines of the output file are displayed when sequence
                        display is on. Also used to determine the number of
                        sequences to average over for speed and time
                        estimates.
  --time_updates_off    turn time updates off (only applies when sequence
                        display is disabled).
  --seqs_per_time_update=N
                        specify the number of sequences between time updates.
                        Default is 1e5.
  -d DELIMITER, --delimiter=DELIMITER
                        declare infile delimiter. Default is tab for .tsv
                        input files, comma for .csv files, and any whitespace
                        for all others. Choices: 'tab', 'space', ',', ';', ':'
  --raw_delimiter=DELIMITER
                        declare infile delimiter as a raw string.
  --delimiter_out=DELIMITER_OUT
                        declare outfile delimiter. Default is tab for .tsv
                        output files, comma for .csv files, and the infile
                        delimiter for all others. Choices: 'tab', 'space',
                        ',', ';', ':'
  --raw_delimiter_out=DELIMITER_OUT
                        declare for the delimiter outfile as a raw string.
  --gene_mask_delimiter=GENE_MASK_DELIMITER
                        declare gene mask delimiter. Default comma unless
                        infile delimiter is comma, then default is a
                        semicolon. Choices: 'tab', 'space', ',', ';', ':'
  --raw_gene_mask_delimiter=GENE_MASK_DELIMITER
                        declare delimiter of gene masks as a raw string.
  --comment_delimiter=COMMENT_DELIMITER
                        character or string to indicate comment or header
                        lines to skip.
--------------------------------------------------------------------------------


@author: zacharysethna
"""

#Function assumes that it is in the same directory that the folder app/ is
#in (which should contain all the modules imported).
from __future__ import print_function, division
import os
import sys
import time
import subprocess
reqs = subprocess.check_output([sys.executable, '-m', 'pip', 'freeze'])
installed_packages = [r.decode().split('==')[0] for r in reqs.split()]
if 'olga' not in installed_packages: sys.path.insert(0, os.path.split(os.path.dirname(os.path.abspath(__file__)))[0])
import olga.load_model as load_model
import olga.generation_probability as generation_probability
from olga.utils import nt2aa, determine_seq_type, gene_to_num_str
#
#import load_model
#import generation_probability
#from utils import nt2aa, determine_seq_type, gene_to_num_str

from optparse import OptionParser

#Set input = raw_input for python 2
try:
    import __builtin__
    input = getattr(__builtin__, 'raw_input')
except (ImportError, AttributeError):
    pass


#Need to determine what mode to run in.
#1) Sequences as args
#2) 'Safe mode', i.e. reading from file, writing to file
#3) not 'safe mode'. i.e. reading from file, printing to stdout


def main():
    """Compute Pgens from a file and output to another file."""

    parser = OptionParser(conflict_handler="resolve")
    parser.add_option('--humanTRA', '--human_T_alpha', action='store_true', dest='humanTRA', default=False, help='use default human TRA model (T cell alpha chain)')
    parser.add_option('--humanTRB', '--human_T_beta', action='store_true', dest='humanTRB', default=False, help='use default human TRB model (T cell beta chain)')
    parser.add_option('--mouseTRB', '--mouse_T_beta', action='store_true', dest='mouseTRB', default=False, help='use default mouse TRB model (T cell beta chain)')
    parser.add_option('--humanIGH', '--human_B_heavy', action='store_true', dest='humanIGH', default=False, help='use default human IGH model (B cell heavy chain)')
    parser.add_option('--mouseTRA', '--mouse_T_alpha', action='store_true', dest='mouseTRA', default=False, help='use default mouse TRA model (T cell alpha chain)')
    parser.add_option('--humanIGL', '--human_B_lambda', action='store_true', dest='humanIGL', default=False, help='use default human IGL model (B cell light lambda chain)')
    parser.add_option('--humanIGK', '--human_B_kappa', action='store_true', dest='humanIGK', default=False, help='use default human IGK model (B cell light kappa chain)')
    parser.add_option('--set_custom_model_VDJ', dest='vdj_model_folder', metavar='PATH/TO/FOLDER/', help='specify PATH/TO/FOLDER/ for a custom VDJ generative model')
    parser.add_option('--set_custom_model_VJ', dest='vj_model_folder', metavar='PATH/TO/FOLDER/', help='specify PATH/TO/FOLDER/ for a custom VJ generative model')

    parser.add_option('-i', '--infile', dest = 'infile_name',metavar='PATH/TO/FILE', help='read in CDR3 sequences (and optionally V/J masks) from PATH/TO/FILE')
    parser.add_option('-o', '--outfile', dest = 'outfile_name', metavar='PATH/TO/FILE', help='write CDR3 sequences and pgens to PATH/TO/FILE')
    parser.add_option('--seq_in', '--seq_index', type='int', metavar='INDEX', dest='seq_in_index', default = 0, help='specifies sequences to be read in are in column INDEX. Default is index 0 (the first column).')

    parser.add_option('--v_in', '--v_mask_index', type='int', metavar='INDEX', dest='V_mask_index', help='specifies V_masks are found in column INDEX in the input file. Default is no V mask.')
    parser.add_option('--j_in', '--j_mask_index', type='int', metavar='INDEX', dest='J_mask_index', help='specifies J_masks are found in column INDEX in the input file. Default is no J mask.')

    parser.add_option('--v_mask', type='string', dest='V_mask', help='specify V usage to condition Pgen on for seqs read in as arguments.')
    parser.add_option('--j_mask', type='string', dest='J_mask', help='specify J usage to condition Pgen on for seqs read in as arguments.')

    parser.add_option('-m', '--max_number_of_seqs', type='int',metavar='N', dest='max_number_of_seqs', help='compute Pgens for at most N sequences.')
    parser.add_option('--lines_to_skip', type='int',metavar='N', dest='lines_to_skip', default = 0, help='skip the first N lines of the file. Default is 0.')
    parser.add_option('-a', '--alphabet_filename', dest='alphabet_filename', metavar='PATH/TO/FILE', help="specify PATH/TO/FILE defining a custom 'amino acid' alphabet. Default is no custom alphabet.")
    parser.add_option('--seq_type_out', type='choice',metavar='SEQ_TYPE', dest='seq_type_out',  choices=['all', 'ntseq', 'nucleotide', 'aaseq', 'amino_acid'], help="if read in sequences are ntseqs, declare what type of sequence to compute pgen for. Default is all. Choices: 'all', 'ntseq', 'nucleotide', 'aaseq', 'amino_acid'")
    parser.add_option('--skip_off','--skip_empty_off', action='store_true', dest = 'skip_empty', default=True, help='stop skipping empty or blank sequences/lines (if for example you want to keep line index fidelity between the infile and outfile).')

    parser.add_option('--display_off', action='store_false', dest='display_seqs', default=True, help='turn the sequence display off (only applies in write-to-file mode). Default is on.')
    parser.add_option('--num_lines_for_display', type='int', metavar='N', default = 50, dest='num_lines_for_display', help='N lines of the output file are displayed when sequence display is on. Also used to determine the number of sequences to average over for speed and time estimates.')
    parser.add_option('--time_updates_off', action='store_false', dest='time_updates', default=True, help='turn time updates off (only applies when sequence display is disabled).')
    parser.add_option('--seqs_per_time_update', type='float', metavar='N', default = 100, dest='seqs_per_time_update', help='specify the number of sequences between time updates. Default is 1e5.')

    parser.add_option('-d', '--delimiter', type='choice', dest='delimiter',  choices=['tab', 'space', ',', ';', ':'], help="declare infile delimiter. Default is tab for .tsv input files, comma for .csv files, and any whitespace for all others. Choices: 'tab', 'space', ',', ';', ':'")
    parser.add_option('--raw_delimiter', type='str', dest='delimiter', help="declare infile delimiter as a raw string.")
    parser.add_option('--delimiter_out', type='choice', dest='delimiter_out',  choices=['tab', 'space', ',', ';', ':'], help="declare outfile delimiter. Default is tab for .tsv output files, comma for .csv files, and the infile delimiter for all others. Choices: 'tab', 'space', ',', ';', ':'")
    parser.add_option('--raw_delimiter_out', type='str', dest='delimiter_out', help="declare for the delimiter outfile as a raw string.")
    parser.add_option('--gene_mask_delimiter', type='choice', dest='gene_mask_delimiter',  choices=['tab', 'space', ',', ';', ':'], help="declare gene mask delimiter. Default comma unless infile delimiter is comma, then default is a semicolon. Choices: 'tab', 'space', ',', ';', ':'")
    parser.add_option('--raw_gene_mask_delimiter', type='str', dest='gene_mask_delimiter', help="declare delimiter of gene masks as a raw string.")
    parser.add_option('--comment_delimiter', type='str', dest='comment_delimiter', help="character or string to indicate comment or header lines to skip.")


    (options, args) = parser.parse_args()

    #Check that the model is specified properly
    main_folder = os.path.dirname(__file__)

    default_models = {}
    default_models['humanTRA'] = [os.path.join(main_folder, 'default_models', 'human_T_alpha'),  'VJ']
    default_models['humanTRB'] = [os.path.join(main_folder, 'default_models', 'human_T_beta'), 'VDJ']
    default_models['mouseTRB'] = [os.path.join(main_folder, 'default_models', 'mouse_T_beta'), 'VDJ']
    default_models['humanIGH'] = [os.path.join(main_folder, 'default_models', 'human_B_heavy'), 'VDJ']
    default_models['humanIGK'] = [os.path.join(main_folder, 'default_models', 'human_B_kappa'), 'VJ']
    default_models['humanIGL'] = [os.path.join(main_folder, 'default_models', 'human_B_lambda'),  'VJ']
    default_models['mouseTRA'] = [os.path.join(main_folder, 'default_models', 'mouse_T_alpha'), 'VJ']


    num_models_specified = sum([1 for x in list(default_models.keys()) + ['vj_model_folder', 'vdj_model_folder'] if getattr(options, x)])

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
        print('Need to indicate generative model.')
        print('Exiting...')
        return -1
    elif num_models_specified > 1:
        print('Only specify one model')
        print('Exiting...')
        return -1

    #Check that all model and genomic files exist in the indicated model folder
    if not os.path.isdir(model_folder):
        print('Check pathing... cannot find the model folder: ' + model_folder)
        print('Exiting...')
        return -1

    params_file_name = os.path.join(model_folder,'model_params.txt')
    marginals_file_name = os.path.join(model_folder,'model_marginals.txt')
    V_anchor_pos_file = os.path.join(model_folder,'V_gene_CDR3_anchors.csv')
    J_anchor_pos_file = os.path.join(model_folder,'J_gene_CDR3_anchors.csv')

    for x in [params_file_name, marginals_file_name, V_anchor_pos_file, J_anchor_pos_file]:
        if not os.path.isfile(x):
            print('Cannot find: ' + x)
            print('Please check the files (and naming conventions) in the model folder ' + model_folder)
            print('Exiting...')
            return -1

    alphabet_filename = options.alphabet_filename #used if a custom alphabet is to be specified
    if alphabet_filename is not None:
        if not os.path.isfile(alphabet_filename):
            print('Cannot find custom alphabet file: ' + alphabet_filename)
            print('Exiting...')
            return -1

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

    aa_alphabet = ''.join(pgen_model.codons_dict.keys())

    if options.infile_name is not None:
        infile_name = options.infile_name

        if not os.path.isfile(infile_name):
            print('Cannot find input file: ' + infile_name)
            print('Exiting...')
            return -1

    if options.outfile_name is not None:
        outfile_name = options.outfile_name
        if os.path.isfile(outfile_name):
            if not input(outfile_name + ' already exists. Overwrite (y/n)? ').strip().lower() in ['y', 'yes']:
                print('Exiting...')
                return -1

    #Parse delimiter
    delimiter = options.delimiter
    if delimiter is None: #Default case
        if options.infile_name is None:
            delimiter = '\t'
        elif infile_name.endswith('.tsv'): #parse TAB separated value file
            delimiter = '\t'
        elif infile_name.endswith('.csv'): #parse COMMA separated value file
            delimiter = ','
    else:
        try:
            delimiter = {'tab': '\t', 'space': ' ', ',': ',', ';': ';', ':': ':'}[delimiter]
        except KeyError:
            pass #Other string passed as the delimiter.

    #Parse delimiter_out
    delimiter_out = options.delimiter_out
    if delimiter_out is None: #Default case
        if delimiter is None:
            delimiter_out = '\t'
        else:
            delimiter_out = delimiter
        if options.outfile_name is None:
            pass
        elif outfile_name.endswith('.tsv'): #output TAB separated value file
            delimiter_out = '\t'
        elif outfile_name.endswith('.csv'): #output COMMA separated value file
            delimiter_out = ','
    else:
        try:
            delimiter_out = {'tab': '\t', 'space': ' ', ',': ',', ';': ';', ':': ':'}[delimiter_out]
        except KeyError:
            pass #Other string passed as the delimiter.

    #Parse gene_delimiter
    gene_mask_delimiter = options.gene_mask_delimiter
    if gene_mask_delimiter is None: #Default case
        gene_mask_delimiter = ','
        if delimiter == ',':
            gene_mask_delimiter = ';'
    else:
        try:
            gene_mask_delimiter = {'tab': '\t', 'space': ' ', ',': ',', ';': ';', ':': ':'}[gene_mask_delimiter]
        except KeyError:
            pass #Other string passed as the delimiter.


    #More options
    time_updates = options.time_updates
    display_seqs = options.display_seqs
    num_lines_for_display = options.num_lines_for_display
    seq_in_index = options.seq_in_index #where in the line the sequence is after line.split(delimiter)
    lines_to_skip = options.lines_to_skip #one method of skipping header
    comment_delimiter = options.comment_delimiter #another method of skipping header
    seqs_per_time_update = options.seqs_per_time_update
    max_number_of_seqs = options.max_number_of_seqs
    V_mask_index = options.V_mask_index #Default is not conditioning on V identity
    J_mask_index = options.J_mask_index #Default is not conditioning on J identity
    skip_empty = options.skip_empty

    seq_type_out = options.seq_type_out #type of pgens to be computed. Can be ntseq, aaseq, or both
    if seq_type_out is not None:
        seq_type_out = {'all': None, 'ntseq': 'ntseq', 'nucleotide': 'ntseq', 'aaseq': 'aaseq', 'amino_acid': 'aaseq'}[seq_type_out]

    if options.infile_name is None: #No infile specified -- args should be the input seqs
        print_warnings = True
        seqs = args
        seq_types = [determine_seq_type(seq, aa_alphabet) for seq in seqs]
        unrecognized_seqs = [seq for i, seq in enumerate(seqs) if seq_types[i] is None]
        if len(unrecognized_seqs) > 0 and print_warnings:
            print('The following sequences/arguments were not recognized: ' + ', '.join(unrecognized_seqs))
        seqs = [seq for i, seq in enumerate(seqs) if seq_types[i] is not None]
        seq_types = [seq_type for seq_type in seq_types if seq_type is not None]


        #Format V and J masks -- uniform for all argument input sequences
        try:
            V_mask = options.V_mask.split(',')
            unrecognized_v_genes = [v for v in V_mask if gene_to_num_str(v, 'V') not in pgen_model.V_mask_mapping.keys()]
            V_mask = [v for v in V_mask if gene_to_num_str(v, 'V') in pgen_model.V_mask_mapping.keys()]
            if len(unrecognized_v_genes) > 0:
                print('These V genes/alleles are not recognized: ' + ', '.join(unrecognized_v_genes))
            if len(V_mask) == 0:
                print('No recognized V genes/alleles in the provided V_mask. Continuing without conditioning on V usage.')
                V_mask = None
        except AttributeError:
            V_mask = options.V_mask #Default is None, i.e. not conditioning on V identity

        try:
            J_mask = options.J_mask.split(',')
            unrecognized_j_genes = [j for j in J_mask if gene_to_num_str(j, 'J') not in pgen_model.J_mask_mapping.keys()]
            J_mask = [j for j in J_mask if gene_to_num_str(j, 'J') in pgen_model.J_mask_mapping.keys()]
            if len(unrecognized_j_genes) > 0:
                print('These J genes/alleles are not recognized: ' + ', '.join(unrecognized_j_genes))
            if len(J_mask) == 0:
                print('No recognized J genes/alleles in the provided J_mask. Continuing without conditioning on J usage.')
                J_mask = None
        except AttributeError:
            J_mask = options.J_mask #Default is None, i.e. not conditioning on J identity

        print('')
        start_time = time.time()
        for seq, seq_type in zip(seqs, seq_types):
            if seq_type == 'aaseq':
                c_pgen = pgen_model.compute_aa_CDR3_pgen(seq, V_mask, J_mask, print_warnings)
                print('Pgen of the amino acid sequence ' + seq + ': ' + str(c_pgen))
                print('')
            elif seq_type == 'regex':
                c_pgen = pgen_model.compute_regex_CDR3_template_pgen(seq, V_mask, J_mask, print_warnings)
                print('Pgen of the regular expression sequence ' + seq + ': ' + str(c_pgen))
                print('')
            elif seq_type == 'ntseq':
                if seq_type_out is None or seq_type_out == 'ntseq':
                    c_pgen_nt = pgen_model.compute_nt_CDR3_pgen(seq, V_mask, J_mask, print_warnings)
                    print('Pgen of the nucleotide sequence ' + seq + ': ' + str(c_pgen_nt))
                if seq_type_out is None or seq_type_out == 'aaseq':
                    c_pgen_aa = pgen_model.compute_aa_CDR3_pgen(nt2aa(seq), V_mask, J_mask, print_warnings)
                    print('Pgen of the amino acid sequence nt2aa(' + seq + ') = ' + nt2aa(seq) + ': ' + str(c_pgen_aa))
                print('')

        c_time = time.time() - start_time
        if c_time > 86400: #more than a day
            c_time_str = '%d days, %d hours, %d minutes, and %.2f seconds.'%(int(c_time)//86400, (int(c_time)//3600)%24, (int(c_time)//60)%60, c_time%60)
        elif c_time > 3600: #more than an hr
            c_time_str = '%d hours, %d minutes, and %.2f seconds.'%((int(c_time)//3600)%24, (int(c_time)//60)%60, c_time%60)
        elif c_time > 60: #more than a min
            c_time_str = '%d minutes and %.2f seconds.'%((int(c_time)//60)%60, c_time%60)
        else:
            c_time_str = '%.2f seconds.'%(c_time)

        print('Completed pgen computation in: ' + c_time_str)

    else: #Read sequences in from file
        print_warnings = False #Most cases of reading in from file should have warnings disabled
        seqs = []
        seq_types = []
        V_usage_masks = []
        J_usage_masks = []

        infile = open(infile_name, 'r')

        for i, line in enumerate(infile):
            if comment_delimiter is not None: #Default case -- no comments/header delimiter
                if line.startswith(comment_delimiter): #allow comments
                    continue
            if i < lines_to_skip:
                continue

            if delimiter is None: #Default delimiter is any whitespace
                split_line = line.split()
            else:
                split_line = line.split(delimiter)

            #Find the seq
            try:
                seq = split_line[seq_in_index].strip()
                if len(seq.strip()) == 0:
                    if skip_empty:
                        continue
                    else:
                        seqs.append(seq) #keep the blank seq as a placeholder
                        seq_types.append('aaseq')
                else:
                    seqs.append(seq)
                    seq_types.append(determine_seq_type(seq, aa_alphabet))
            except IndexError: #no index match for seq
                if skip_empty and len(line.strip()) == 0:
                    continue
                print('seq_in_index is out of range')
                print('Exiting...')
                infile.close()
                return -1

            #Find and format V_usage_mask
            if V_mask_index is None:
                V_usage_masks.append(None) #default mask
            else:
                try:
                    V_usage_mask = split_line[V_mask_index].strip().split(gene_mask_delimiter)
                    #check that all V gene/allele names are recognized
                    if all([gene_to_num_str(v, 'V') in pgen_model.V_mask_mapping for v in V_usage_mask]):
                        V_usage_masks.append(V_usage_mask)
                    else:
                        print(str(V_usage_mask) + " is not a usable V_usage_mask composed exclusively of recognized V gene/allele names")
                        print('Unrecognized V gene/allele names: ' + ', '.join([v for v in V_usage_mask if gene_to_num_str(v, 'V') not in pgen_model.V_mask_mapping.keys()]))
                        print('Exiting...')
                        infile.close()
                        return -1
                except IndexError: #no index match for V_mask_index
                    print('V_mask_index is out of range')
                    print('Exiting...')
                    infile.close()
                    return -1

            #Find and format J_usage_mask
            if J_mask_index is None:
                J_usage_masks.append(None) #default mask
            else:
                try:
                    J_usage_mask = split_line[J_mask_index].strip().split(gene_mask_delimiter)
                    #check that all V gene/allele names are recognized
                    if all([gene_to_num_str(j, 'J') in pgen_model.J_mask_mapping for j in J_usage_mask]):
                        J_usage_masks.append(J_usage_mask)
                    else:
                        print(str(J_usage_mask) + " is not a usable J_usage_mask composed exclusively of recognized J gene/allele names")
                        print('Unrecognized J gene/allele names: ' + ', '.join([j for j in J_usage_mask if gene_to_num_str(j, 'J') not in pgen_model.J_mask_mapping.keys()]))
                        print('Exiting...')
                        infile.close()
                        return -1
                except IndexError: #no index match for J_mask_index
                    print('J_mask_index is out of range')
                    print('Exiting...')
                    infile.close()
                    return -1

            if max_number_of_seqs is not None:
                if len(seqs) >= max_number_of_seqs:
                    break


        unrecognized_seqs = [seq for i, seq in enumerate(seqs) if seq_types[i] is None]
        if len(unrecognized_seqs) > 0 and len(unrecognized_seqs) < len(seqs):
            if print_warnings or options.outfile_name is not None:
                print('Some strings read in were not parsed as sequences -- they will be omitted.')
                print('Examples of improperly read strings: ')
                for unrecognized_seq in unrecognized_seqs[:10]:
                    print(unrecognized_seq)
            seqs = [seq for i, seq in enumerate(seqs) if seq_types[i] is not None]
            V_usage_masks = [V_usage_mask for i, V_usage_mask in enumerate(V_usage_masks) if seq_types[i] is not None]
            J_usage_masks = [J_usage_mask for i, J_usage_mask in enumerate(J_usage_masks) if seq_types[i] is not None]
            seq_types = [seq_type for seq_type in seq_types if seq_type is not None]
        elif len(unrecognized_seqs) > 0 and len(unrecognized_seqs) == len(seqs):
            print('None of the read in strings were parsed as sequences. Check input file.')
            print('Examples of improperly read strings:')
            for unrecognized_seq in unrecognized_seqs[:10]:
                print(unrecognized_seq)
            print('Exiting...')
            return -1

        infile.close()


        if options.outfile_name is not None: #OUTFILE SPECIFIED, allow printed info/display

            print('Successfully read in and formatted ' + str(len(seqs)) + ' sequences and any V or J usages.')
            if display_seqs:
                sys.stdout.write('\r'+'Continuing to Pgen computation in 3... ')
                sys.stdout.flush()
                time.sleep(0.4)
                sys.stdout.write('\r'+'Continuing to Pgen computation in 2... ')
                sys.stdout.flush()
                time.sleep(0.4)
                sys.stdout.write('\r'+'Continuing to Pgen computation in 1... ')
                sys.stdout.flush()
                time.sleep(0.4)
            else:
                print('Continuing to Pgen computation.')
                print_warnings = True #Display is off, can print warnings

            if display_seqs:
                lines_for_display = []
                times_for_speed_calc = [time.time()]

            outfile = open(outfile_name, 'w')
            start_time = time.time()
            for i, seq in enumerate(seqs):
                if seq_types[i] == 'aaseq':
                    #Compute Pgen and print out
                    c_pgen_line = seq + delimiter_out + str(pgen_model.compute_aa_CDR3_pgen(seq, V_usage_masks[i], J_usage_masks[i], print_warnings))
                if seq_types[i] == 'regex':
                    #Compute Pgen and print out
                    c_pgen_line = seq + delimiter_out + str(pgen_model.compute_regex_CDR3_template_pgen(seq, V_usage_masks[i], J_usage_masks[i], print_warnings))
                elif seq_types[i] == 'ntseq':
                    ntseq = seq
                    if len(ntseq) % 3 == 0: #inframe sequence
                        aaseq = nt2aa(ntseq)
                        #Compute Pgen and print out based on recomb_type and seq_type_out
                        if seq_type_out is None:
                            c_pgen_line = ntseq + delimiter_out + str(pgen_model.compute_nt_CDR3_pgen(ntseq, V_usage_masks[i], J_usage_masks[i], print_warnings)) + delimiter_out + aaseq + delimiter_out +  str(pgen_model.compute_aa_CDR3_pgen(aaseq, V_usage_masks[i], J_usage_masks[i], print_warnings))
                        elif seq_type_out == 'ntseq':
                            c_pgen_line = ntseq + delimiter_out + str(pgen_model.compute_nt_CDR3_pgen(ntseq, V_usage_masks[i], J_usage_masks[i], print_warnings))
                        elif seq_type_out == 'aaseq':
                            c_pgen_line = aaseq + delimiter_out + str(pgen_model.compute_aa_CDR3_pgen(aaseq, V_usage_masks[i], J_usage_masks[i], print_warnings))
                    else: #out of frame sequence -- Pgens are 0 and use 'out_of_frame' for aaseq
                        if seq_type_out is None:
                            c_pgen_line = ntseq + delimiter_out + '0' + delimiter_out + 'out_of_frame' + delimiter_out + '0'
                        elif seq_type_out == 'ntseq':
                            c_pgen_line = ntseq + delimiter_out + '0'
                        elif seq_type_out == 'aaseq':
                            c_pgen_line = 'out_of_frame' + delimiter_out + '0'

                outfile.write(c_pgen_line + '\n')

                #Print time update
                if display_seqs:
                    cc_time = time.time()
                    c_time = cc_time - start_time
                    times_for_speed_calc = [cc_time] + times_for_speed_calc[:num_lines_for_display]
                    c_avg_speed = (len(times_for_speed_calc)-1)/(times_for_speed_calc[0] - times_for_speed_calc[-1])

                    #eta = ((len(seqs) - (i+1))/float(i+1))*c_time

                    eta = (len(seqs) - (i+1))/c_avg_speed

                    lines_for_display = [c_pgen_line] + lines_for_display[:num_lines_for_display]


                    c_time_str = '%s hours, %s minutes, and %s seconds.'%(repr(int(c_time)//3600).rjust(3), repr((int(c_time)//60)%60).rjust(2), repr(int(c_time)%60).rjust(2))
                    eta_str = '%s hours, %s minutes, and %s seconds.'%(repr(int(eta)//3600).rjust(3), repr((int(eta)//60)%60).rjust(2), repr(int(eta)%60).rjust(2))
                    time_str = 'Time to compute Pgen on %s seqs: %s \nEst. time for remaining %s seqs: %s'%(repr(i+1).rjust(9), c_time_str, repr(len(seqs) - (i + 1)).rjust(9), eta_str)
                    speed_str = 'Current Pgen computation speed: %s seqs/min'%(repr(round((len(times_for_speed_calc)-1)*60/(times_for_speed_calc[0] - times_for_speed_calc[-1]), 2)).rjust(8))
                    display_str = '\n'.join(lines_for_display[::-1]) + '\n' + '-'*80 + '\n' + time_str + '\n' + speed_str + '\n' + '-'*80
                    print('\033[2J' + display_str)
                elif (i+1)%seqs_per_time_update == 0 and time_updates:
                    c_time = time.time() - start_time
                    eta = ((len(seqs) - (i+1))/(i+1))*c_time
                    if c_time > 86400: #more than a day
                        c_time_str = '%d days, %d hours, %d minutes, and %.2f seconds.'%(int(c_time)//86400, (int(c_time)//3600)%24, (int(c_time)//60)%60, c_time%60)
                    elif c_time > 3600: #more than an hr
                        c_time_str = '%d hours, %d minutes, and %.2f seconds.'%((int(c_time)//3600)%24, (int(c_time)//60)%60, c_time%60)
                    elif c_time > 60: #more than a min
                        c_time_str = '%d minutes and %.2f seconds.'%((int(c_time)//60)%60, c_time%60)
                    else:
                        c_time_str = '%.2f seconds.'%(c_time)

                    if eta > 86400: #more than a day
                        eta_str = '%d days, %d hours, %d minutes, and %.2f seconds.'%(int(eta)//86400, (int(eta)//3600)%24, (int(eta)//60)%60, eta%60)
                    elif eta > 3600: #more than an hr
                        eta_str = '%d hours, %d minutes, and %.2f seconds.'%((int(eta)//3600)%24, (int(eta)//60)%60, eta%60)
                    elif eta > 60: #more than a min
                        eta_str = '%d minutes and %.2f seconds.'%((int(eta)//60)%60, eta%60)
                    else:
                        eta_str = '%.2f seconds.'%(eta)

                    print('Pgen computed for %d sequences in: %s Estimated time remaining: %s'%(i+1, c_time_str, eta_str))

            c_time = time.time() - start_time
            if c_time > 86400: #more than a day
                c_time_str = '%d days, %d hours, %d minutes, and %.2f seconds.'%(int(c_time)//86400, (int(c_time)//3600)%24, (int(c_time)//60)%60, c_time%60)
            elif c_time > 3600: #more than an hr
                c_time_str = '%d hours, %d minutes, and %.2f seconds.'%((int(c_time)//3600)%24, (int(c_time)//60)%60, c_time%60)
            elif c_time > 60: #more than a min
                c_time_str = '%d minutes and %.2f seconds.'%((int(c_time)//60)%60, c_time%60)
            else:
                c_time_str = '%.2f seconds.'%(c_time)
            print('Completed Pgen computation for %d sequences: in %s'%(len(seqs), c_time_str))

            outfile.close()

        else: #NO OUTFILE -- print directly to stdout
            start_time = time.time()
            for i, seq in enumerate(seqs):
                if seq_types[i] == 'aaseq':
                    #Compute Pgen and print out
                    c_pgen_line = seq + delimiter_out + str(pgen_model.compute_aa_CDR3_pgen(seq, V_usage_masks[i], J_usage_masks[i], print_warnings))
                if seq_types[i] == 'regex':
                    #Compute Pgen and print out
                    c_pgen_line = seq + delimiter_out + str(pgen_model.compute_regex_CDR3_template_pgen(seq, V_usage_masks[i], J_usage_masks[i], print_warnings))
                elif seq_types[i] == 'ntseq':
                    ntseq = seq
                    if len(ntseq) % 3 == 0: #inframe sequence
                        aaseq = nt2aa(ntseq)
                        #Compute Pgen and print out based on recomb_type and seq_type_out
                        if seq_type_out is None:
                            c_pgen_line = ntseq + delimiter_out + str(pgen_model.compute_nt_CDR3_pgen(ntseq, V_usage_masks[i], J_usage_masks[i], print_warnings)) + delimiter_out + aaseq + delimiter_out +  str(pgen_model.compute_aa_CDR3_pgen(aaseq, V_usage_masks[i], J_usage_masks[i], print_warnings))
                        elif seq_type_out == 'ntseq':
                            c_pgen_line = ntseq + delimiter_out + str(pgen_model.compute_nt_CDR3_pgen(ntseq, V_usage_masks[i], J_usage_masks[i], print_warnings))
                        elif seq_type_out == 'aaseq':
                            c_pgen_line = aaseq + delimiter_out + str(pgen_model.compute_aa_CDR3_pgen(aaseq, V_usage_masks[i], J_usage_masks[i], print_warnings))
                    else: #out of frame sequence -- Pgens are 0 and use 'out_of_frame' for aaseq
                        if seq_type_out is None:
                            c_pgen_line = ntseq + delimiter_out + '0' + delimiter_out + 'out_of_frame' + delimiter_out + '0'
                        elif seq_type_out == 'ntseq':
                            c_pgen_line = ntseq + delimiter_out + '0'
                        elif seq_type_out == 'aaseq':
                            c_pgen_line = 'out_of_frame' + delimiter_out + '0'

                print(c_pgen_line)

if __name__ == '__main__': main()
