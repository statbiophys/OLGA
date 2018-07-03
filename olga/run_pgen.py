#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""Command line script to compute Pgens of sequences from file.

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

This program will compute  the generation probabilities (Pgens) of sequences
read in from a file as defined by a specified generative model.
This program has only minimal parsing built in, so the file read in must be
structured with delimiter spaced values (i.e. the data is organized in columns
separated by delimiter like a .tsv or .csv file). Read the Options on delimiter
for more info.

The default is to assume that the sequences to be read in are in the first
column (index 0), meaning that a text file with only a sequence on each line
will be read in okay by default.

The sequences which are read in must be TRIMMED TO ONLY THE CDR3 region as
defined by the V and J anchor files (default is to INCLUDE the conserved
residues of the C in the V region and the F/W in the J region).

It is also possible to condition the Pgen computation on V and J identity by
specifying what index column the V and J gene/allele information is stored.
Look at the example files/calls to see some examples of this syntax.

The default is for the program to run in 'safe mode' where all of the sequences
and V/J masks are loaded up and checked before running. This mode also has a
display (default on) that will show some number of lines of the outfile as a
visual check that the program is running properly (it also includes time, speed
and progress). These can be disabled.

This program runs on only one type of sequence (either nucleotide or amino
acid), which will be inferred by default in safe mode. Regular expression
sequences are not accepted in this program as they can incur large time costs
(instead consider defining a custom 'amino acid' alphabet to define the symbols
used in the regular expressions if possible, or use the single sequence
function). If nucleotide sequences are read in it is possible to specify if the
output should be the nucleotide sequence Pgen and/or the translated amino acid
sequence Pgen (the default is to compute and output both).

If not running in safe mode, sequences will be read in one by one and Pgen
computed/written to the outfile as they are read in (so only one sequence is
held in memory at a time). In this mode, it is required to specify whether the
sequences to be read in are nucleotide sequences or 'amino acid' sequences.

As it is rare for datasets to be >> 1e4, parallelization is not built in.
However, there are options to skip N lines of the file and to load at most M
sequences so, if wanted, one could build a parallelized wrapper around this
script (though it would be recommended to instead just import the modules and
build from there).

Required inputs:

1) File name for input file: PATH/TO/INFILE
THIS IS ASSUMED TO BE THE FIRST ARGUMENT (see example calls)

2) File name for output file: PATH/TO/OUTFILE
THIS IS ASSUMED TO BE THE SECOND ARGUMENT (see example calls)

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
Example calls (example calls are formatted as executed functions instead of the
console script entry point. The arguments are identical in either case.). These
calls assume the existance of example_seqs.tsv.

./run_pgen.py example_seqs.tsv example_pgens.tsv --humanTCRB
./run_pgen.py example_seqs.tsv example_pgens.tsv --VDJ_model_folder models/human_T_beta/

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
  --seq_in_index=INDEX  specifies sequences to be read in are in column INDEX.
                        Default is index 0 (the first column).
  -v INDEX, --v_mask_index=INDEX
                        specifies V_masks are found in column INDEX. Default
                        is no V mask.
  -j INDEX, --j_mask_index=INDEX
                        specifies J_masks are found in column INDEX. Default
                        is no J mask.
  -m N, --max_number_of_seqs=N
                        compute Pgens for at most N sequences.
  --lines_to_skip=N     skip the first N lines of the file. Default is 0.
  -a PATH/TO/FILE, --alphabet_filename=PATH/TO/FILE
                        specify PATH/TO/FILE defining a custom 'amino acid'
                        alphabet. Default is no custom alphabet.
  --seq_type=SEQ_TYPE   declare sequence type. Infers seq_type by default in
                        safe_mode. Need to declare when safe_mode is off.
                        Choices: 'ntseq', 'nucleotide', 'aaseq', 'amino_acid'
  --seq_type_out=SEQ_TYPE
                        if read in sequences are ntseqs, declare what type of
                        sequence to compute pgen for. Default is all. Choices:
                        'all', 'ntseq', 'nucleotide', 'aaseq', 'amino_acid'
  --skip_off, --skip_empty_sequences_off
                        stop skipping empty or blank sequences/lines (if for
                        example you want to keep line index fidelity between
                        the infile and outfile).
  --safe_mode_off       turn safe_mode off. Default is on.
  --display_seqs_off    turn the sequence display off (only applies in safe
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
  --print_warnings_off  turn Pgen print warnings off
  --overwrite_on        overwrites outfile without prompt when safe mode is
                        turned off.
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
-------------------------------------------------------------------------------
Example calls with options

#Only run on the first 250 seqs (all the other examples will also be restricted to 250 seqs)
./run_pgen.py example_seqs.tsv example_pgens.tsv --humanTCRB -m 250

#Read in amino acid sequences (which happen to be in the column of index 1 (second column) in example_seqs.txt)
./run_pgen.py example_seqs.tsv example_pgens.tsv --humanTCRB --seq_in_index 1 -m 250

#Specify the custom alphabet file expanded_alphabet_files/example_expanded_amino_acid_alphabet.txt
#(Note, to take advantage of the expanded 'amino acid' alphabet need to read in amino acid sequences
#as the nucleotide sequences always translate to the standard amino acids.)
./run_pgen.py example_seqs.tsv example_pgens.tsv --humanTCRB -a expanded_alphabet_files/example_expanded_amino_acid_alphabet.txt -m 250

#Read in nucleotide sequences but only run Pgen on the translated amino acid sequences
./run_pgen.py example_seqs.tsv example_pgens.tsv --humanTCRB --seq_type ntseq --seq_type_out aaseq -m 250

#Compute Pgen conditioned on the V and J gene/alleles provided in columns of index 2 and 3.
./run_pgen.py example_seqs.tsv example_pgens.tsv --humanTCRB -v 2 -j 3 -m 250

#Turn the sequence display off
./run_pgen.py example_seqs.tsv example_pgens.tsv --humanTCRB --display_seqs_off -m 250

#Change sequence display to 20 sequences (default is 50)
./run_pgen.py example_seqs.tsv example_pgens.tsv --humanTCRB --num_lines_for_display 20 -m 250

#Skip the first 5 lines
./run_pgen.py example_seqs.tsv example_pgens.tsv --humanTCRB --lines_to_skip 5 -m 250

#Safe mode off
./run_pgen.py example_seqs.tsv example_pgens.tsv --humanTCRB --seq_type ntseq --safe_mode_off -m 250


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
    """Compute Pgens from a file and output to another file."""

    parser = OptionParser(conflict_handler="resolve")

    parser.add_option('--humanTCRA', '--human_T_alpha', action='store_true', dest='humanTCRA', default=False, help='use default human TCRA model (T cell alpha chain)')
    parser.add_option('--humanTCRB', '--human_T_beta', action='store_true', dest='humanTCRB', default=False, help='use default human TCRB model (T cell beta chain)')
    parser.add_option('--mouseTCRB', '--mouse_T_beta', action='store_true', dest='mouseTCRB', default=False, help='use default mouse TCRB model (T cell beta chain)')
    parser.add_option('--humanIGH', '--human_B_heavy', action='store_true', dest='humanIGH', default=False, help='use default human IGH model (B cell heavy chain)')
    parser.add_option('--VDJ_model_folder', dest='vdj_model_folder', metavar='PATH/TO/FOLDER/', help='specify PATH/TO/FOLDER/ for a custom VDJ generative model')
    parser.add_option('--VJ_model_folder', dest='vj_model_folder', metavar='PATH/TO/FOLDER/', help='specify PATH/TO/FOLDER/ for a custom VJ generative model')
    parser.add_option('--seq_in_index', type='int', metavar='INDEX', dest='seq_in_index', default = 0, help='specifies sequences to be read in are in column INDEX. Default is index 0 (the first column).')
    parser.add_option('-v', '--v_mask_index', type='int', metavar='INDEX', dest='V_mask_index', help='specifies V_masks are found in column INDEX. Default is no V mask.')
    parser.add_option('-j', '--j_mask_index', type='int', metavar='INDEX', dest='J_mask_index', help='specifies J_masks are found in column INDEX. Default is no J mask.')
    parser.add_option('-m', '--max_number_of_seqs', type='int',metavar='N', dest='max_number_of_seqs', help='compute Pgens for at most N sequences.')
    parser.add_option('--lines_to_skip', type='int',metavar='N', dest='lines_to_skip', default = 0, help='skip the first N lines of the file. Default is 0.')
    parser.add_option('-a', '--alphabet_filename', dest='alphabet_filename', metavar='PATH/TO/FILE', help="specify PATH/TO/FILE defining a custom 'amino acid' alphabet. Default is no custom alphabet.")
    parser.add_option('--seq_type', type='choice', dest='seq_type',  choices=['ntseq', 'nucleotide', 'aaseq', 'amino_acid'], help="declare sequence type. Infers seq_type by default in safe_mode. Need to declare when safe_mode is off. Choices: 'ntseq', 'nucleotide', 'aaseq', 'amino_acid'")
    parser.add_option('--seq_type_out', type='choice',metavar='SEQ_TYPE', dest='seq_type_out',  choices=['all', 'ntseq', 'nucleotide', 'aaseq', 'amino_acid'], help="if read in sequences are ntseqs, declare what type of sequence to compute pgen for. Default is all. Choices: 'all', 'ntseq', 'nucleotide', 'aaseq', 'amino_acid'")
    parser.add_option('--skip_off','--skip_empty_sequences_off', action='store_true', dest = 'skip_empty_sequences', default=True, help='stop skipping empty or blank sequences/lines (if for example you want to keep line index fidelity between the infile and outfile).')

    parser.add_option('--safe_mode_off', action='store_false', dest='safe_mode', default=True, help='turn safe_mode off. Default is on.')
    parser.add_option('--display_seqs_off', action='store_false', dest='display_seqs', default=True, help='turn the sequence display off (only applies in safe mode). Default is on.')
    parser.add_option('--num_lines_for_display', type='int', metavar='N', default = 50, dest='num_lines_for_display', help='N lines of the output file are displayed when sequence display is on. Also used to determine the number of sequences to average over for speed and time estimates.')
    parser.add_option('--time_updates_off', action='store_false', dest='time_updates', default=True, help='turn time updates off (only applies when sequence display is disabled).')
    parser.add_option('--seqs_per_time_update', type='float', metavar='N', default = 100, dest='seqs_per_time_update', help='specify the number of sequences between time updates. Default is 1e5.')
    parser.add_option('--print_warnings_off', action='store_false', dest="print_warnings", default=True, help='turn Pgen print warnings off')
    parser.add_option('--overwrite_on', action='store_true', dest="overwrite", default=False, help='overwrites outfile without prompt when safe mode is turned off.')

    parser.add_option('-d', '--delimiter', type='choice', dest='delimiter',  choices=['tab', 'space', ',', ';', ':'], help="declare infile delimiter. Default is tab for .tsv input files, comma for .csv files, and any whitespace for all others. Choices: 'tab', 'space', ',', ';', ':'")
    parser.add_option('--raw_delimiter', type='str', dest='delimiter', help="declare infile delimiter as a raw string.")
    parser.add_option('--delimiter_out', type='choice', dest='delimiter_out',  choices=['tab', 'space', ',', ';', ':'], help="declare outfile delimiter. Default is tab for .tsv output files, comma for .csv files, and the infile delimiter for all others. Choices: 'tab', 'space', ',', ';', ':'")
    parser.add_option('--raw_delimiter_out', type='str', dest='delimiter_out', help="declare for the delimiter outfile as a raw string.")
    parser.add_option('--gene_mask_delimiter', type='choice', dest='gene_mask_delimiter',  choices=['tab', 'space', ',', ';', ':'], help="declare gene mask delimiter. Default comma unless infile delimiter is comma, then default is a semicolon. Choices: 'tab', 'space', ',', ';', ':'")
    parser.add_option('--raw_gene_mask_delimiter', type='str', dest='gene_mask_delimiter', help="declare delimiter of gene masks as a raw string.")
    parser.add_option('--comment_delimiter', type='str', dest='comment_delimiter', help="character or string to indicate comment or header lines to skip.")


    (options, args) = parser.parse_args()

    #INFILE IS THE FIRST ARGUMENT
    try:
        infile_name = args[0]
    except IndexError:
        print 'Need to specify PATH/TO/INFILE as the first argument!'
        print 'Exiting...'
        return -1

    if len(infile_name.strip()) == 0:
        print 'Need to specify PATH/TO/INFILE as the first argument!'
        print 'Exiting...'
        return -1

    if not os.path.isfile(infile_name):
        print 'Cannot find file: ' + infile_name
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


    safe_mode = options.safe_mode
    overwrite = options.overwrite
    #OUTFILE IS THE SECOND ARGUMENT
    try:
        outfile_name = args[1]
    except IndexError:
        print 'Need to specify PATH/TO/OUTFILE as the second argument!'
        print 'Exiting...'
        return -1

    if len(outfile_name.strip()) == 0:
        print 'Need to specify PATH/TO/OUTFILE as the second argument!'
        print 'Exiting...'
        return -1

    if os.path.isfile(outfile_name):
        if safe_mode:
            if not raw_input(outfile_name + ' already exists. Overwrite (y/n)? ').strip().lower() in ['y', 'yes']:
                print 'Exiting...'
                return -1
        else:
            if not overwrite:
                print outfile_name + ' already exists! To overwrite use --overwrite_on'
                print 'Exiting...'
                return -1




    #Parse delimiter
    delimiter = options.delimiter
    if delimiter is None: #Default case
        if infile_name.endswith('.tsv'): #parse TAB separated value file
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
        if outfile_name.endswith('.tsv'): #output TAB separated value file
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
    alphabet_filename = options.alphabet_filename #used if a custom alphabet is to be specified
    seq_type_in = options.seq_type
    print_warnings = options.print_warnings
    time_updates = options.time_updates
    display_seqs = options.display_seqs
    num_lines_for_display = options.num_lines_for_display
    seq_type_out = options.seq_type_out #type of pgens to be computed. Can be ntseq, aaseq, or both
    seq_in_index = options.seq_in_index #where in the line the sequence is after line.split(delimiter)
    lines_to_skip = options.lines_to_skip #one method of skipping header
    comment_delimiter = options.comment_delimiter #another method of skipping header
    seqs_per_time_update = options.seqs_per_time_update
    max_number_of_seqs = options.max_number_of_seqs
    V_mask_index = options.V_mask_index #Default is not conditioning on V identity
    J_mask_index = options.J_mask_index #Default is not conditioning on J identity
    skip_empty_sequences = options.skip_empty_sequences


    #Set seq_types
    if seq_type_in is not None:
        seq_type_in = {'ntseq': 'ntseq', 'nucleotide': 'ntseq', 'aaseq': 'aaseq', 'amino_acid': 'aaseq'}[seq_type_in]
    if seq_type_out is not None:
        seq_type_out = {'all': None, 'ntseq': 'ntseq', 'nucleotide': 'ntseq', 'aaseq': 'aaseq', 'amino_acid': 'aaseq'}[seq_type_out]


    if seq_type_in == 'aaseq' and seq_type_out == 'ntseq':
        print "Can't compute nucleotide Pgen from an amino acid sequence (seq_type_in = aaseq, seq_type_out = ntseq)"
        print 'Exiting...'
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

    #SAFE MODE
    if safe_mode:
        seqs = []
        V_usage_masks = []
        J_usage_masks = []

        infile = open(infile_name, 'r')

        for i, line in enumerate(infile):
            if comment_delimiter is not None: #Default case -- no comments/header delimiter
                if line.startswith(comment_delimiter): #allow comments
                    continue
            if i < lines_to_skip:
                continue

            if delimiter is None: #Default delimiter is just any whitespace
                split_line = line.split()
            else:
                split_line = line.split(delimiter)

            #Find seq
            try:
                seq = split_line[seq_in_index].strip()
                if skip_empty_sequences and len(seq.strip()) == 0:
                    continue
                if seq_type_in == 'aaseq':
                    #check seq is composed only of symbols that are in the codons_dict
                    if all([x in pgen_model.codons_dict.keys() for x in seq]):
                        seqs.append(seq)
                    else:
                        print seq + " is not a CDR3 'amino acid' sequence composed exclusively of recognized symbols"
                        print 'Unrecognized symbols: ' + ', '.join([x for x in seq if not x in pgen_model.codons_dict.keys()])
                        print 'Exiting...'
                        infile.close()
                        return -1
                elif seq_type_in == 'ntseq':
                    #check that seq is only composed of nucleotides
                    if all([x in 'ACGTacgt' for x in seq]):
                        if len(seq)%3 == 0:
                            seqs.append(seq)
                        else:
                            #Out of frame sequence -- print warning.
                            print 'WARNING: ' + seq + ' is out of frame --- cannot run Pgen for it (will output 0).'
                            seqs.append(seq)
                    else:
                        print seq + " is not a CDR3 nucleotide sequence composed exclusively of A, C, G, or T."
                        print 'Unrecognized symbols: ' + ', '.join([x for x in seq if not x in 'ACGTacgt'])
                        print 'Exiting...'
                        infile.close()
                        return -1
                else: #Will infer seq_type_in and vet sequences later
                    seqs.append(seq)

            except IndexError: #no index match for seq
                if skip_empty_sequences:
                    if len(line.strip()) == 0:
                        continue
                print 'seq_in_index is out of range'
                print 'Exiting...'
                infile.close()
                return -1

            #Find and format V_usage_mask
            if V_mask_index is None:
                V_usage_masks.append(None) #default mask
            else:
                try:
                    V_usage_mask = split_line[V_mask_index].strip().split(gene_mask_delimiter)
                    #check that all V gene/allele names are recognized
                    if all([v in pgen_model.V_mask_mapping for v in V_usage_mask]):
                        V_usage_masks.append(V_usage_mask)
                    else:
                        print str(V_usage_mask) + " is not a usable V_usage_mask composed exclusively of recognized V gene/allele names"
                        print 'Unrecognized V gene/allele names: ' + ', '.join([v for v in V_usage_mask if not v in pgen_model.V_mask_mapping.keys()])
                        print 'Exiting...'
                        infile.close()
                        return -1
                except IndexError: #no index match for V_mask_index
                    print 'V_mask_index is out of range'
                    print 'Exiting...'
                    infile.close()
                    return -1

            #Find and format J_usage_mask
            if J_mask_index is None:
                J_usage_masks.append(None) #default mask
            else:
                try:
                    J_usage_mask = split_line[J_mask_index].strip().split(gene_mask_delimiter)
                    #check that all V gene/allele names are recognized
                    if all([j in pgen_model.J_mask_mapping for j in J_usage_mask]):
                        J_usage_masks.append(J_usage_mask)
                    else:
                        print str(J_usage_mask) + " is not a usable J_usage_mask composed exclusively of recognized J gene/allele names"
                        print 'Unrecognized J gene/allele names: ' + ', '.join([j for j in J_usage_mask if not j in pgen_model.J_mask_mapping.keys()])
                        print 'Exiting...'
                        infile.close()
                        return -1
                except IndexError: #no index match for J_mask_index
                    print 'J_mask_index is out of range'
                    print 'Exiting...'
                    infile.close()
                    return -1

            if max_number_of_seqs is not None:
                if len(seqs) >= max_number_of_seqs:
                    break

        infile.close()

        if seq_type_in == 'aaseq': #Check if the 'aaseqs' read in are actually 'ntseqs'
            all_symbols_nts = True
            for seq in seqs:
                if not all([x in 'ACGTacgt' for x in seq]):
                    all_symbols_nts = False
                    break
            if all_symbols_nts:
                print 'It looks like the sequences read in are all nucleotide sequences composed of A, C, G, or T.'
                if raw_input('Change seq_type_in to nucleotide? (y/n)? ') in ['y', 'yes']:
                    seq_type_in = 'ntseq'
                else:
                    print 'Okay... you are the boss... chances are all the Pgens will be 0.'

        elif seq_type_in is None: #no specified seq_type_in --- infer if possible
            ntseq_poss = True
            aaseq_poss = True
            ntseq_found = False
            for seq in seqs:
                if len(seq) == 0: #skip empty entry
                    continue
                if all([x in 'ACGTacgt' for x in seq]): #ntseq
                    ntseq_found = True
                elif all([x in pgen_model.codons_dict.keys() for x in seq]): #aaseq that CAN'T be ntseq
                    ntseq_poss = False
                else: #Can't be either a ntseq or aaseq
                    ntseq_poss = False
                    aaseq_poss = False
                    break

            if ntseq_poss: #If nucleotide sequence is possible set seq_type_in to nucleotide
                seq_type_in = 'ntseq'
            elif aaseq_poss:
                if not ntseq_found: #No sequences that appear to be nucleotide sequences found
                    seq_type_in = 'aaseq'
                else:
                    print "It looks like the some sequences are 'amino acid' sequences, while others are nucleotide sequences."
                    if raw_input('Compute amino acid Pgen for all sequences (including nucleotide sequences)? y/n (if no will break)? ') in ['y', 'yes']:
                        seq_type_in = 'aaseq'
                    else:
                       print 'Exiting...'
                       return -1
            else:
                print 'Sequence ' + seq + " is unrecognized as either an 'amino acid' sequence or a nucleotide sequence."
                print 'Exiting...'
                return -1

            if seq_type_in == 'aaseq' and seq_type_out == 'ntseq':
                print "Can't compute nucleotide Pgen from an amino acid sequence (seq_type_in = aaseq, seq_type_out = ntseq)"
                print 'Exiting...'
                return -1

        print 'Successfully read in and formatted ' + str(len(seqs)) + ' sequences and any V or J usages.'
        if display_seqs:
            print_warnings = False #Can't print warnings with display on
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
            print 'Continuing to Pgen computation.'

        if display_seqs:
            lines_for_display = []
            times_for_speed_calc = [time.time()]

        outfile = open(outfile_name, 'w')
        start_time = time.time()
        for i, seq in enumerate(seqs):
            if seq_type_in == 'aaseq':
                aaseq = seq
                #Compute Pgen and print out
                c_pgen_line = aaseq + delimiter_out + str(pgen_model.compute_aa_CDR3_pgen(aaseq, V_usage_masks[i], J_usage_masks[i], print_warnings))
            elif seq_type_in == 'ntseq':
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
                c_avg_speed = (len(times_for_speed_calc)-1)/float(times_for_speed_calc[0] - times_for_speed_calc[-1])

                #eta = ((len(seqs) - (i+1))/float(i+1))*c_time

                eta = (len(seqs) - (i+1))/c_avg_speed

                lines_for_display = [c_pgen_line] + lines_for_display[:num_lines_for_display]


                c_time_str = '%s hours, %s minutes, and %s seconds.'%(repr(int(c_time)/3600).rjust(3), repr((int(c_time)/60)%60).rjust(2), repr(int(c_time)%60).rjust(2))
                eta_str = '%s hours, %s minutes, and %s seconds.'%(repr(int(eta)/3600).rjust(3), repr((int(eta)/60)%60).rjust(2), repr(int(eta)%60).rjust(2))
                time_str = 'Time to compute Pgen on %s seqs: %s \nEst. time for remaining %s seqs: %s'%(repr(i+1).rjust(9), c_time_str, repr(len(seqs) - (i + 1)).rjust(9), eta_str)
                speed_str = 'Current Pgen computation speed: %s seqs/min'%(repr(round((len(times_for_speed_calc)-1)*60/float(times_for_speed_calc[0] - times_for_speed_calc[-1]), 2)).rjust(8))
                display_str = '\n'.join(lines_for_display[::-1]) + '\n' + '-'*80 + '\n' + time_str + '\n' + speed_str + '\n' + '-'*80
                print '\033[2J' + display_str
            elif (i+1)%seqs_per_time_update == 0 and time_updates:
                c_time = time.time() - start_time
                eta = ((len(seqs) - (i+1))/float(i+1))*c_time
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

                print 'Pgen computed for %d sequences in: %s Estimated time remaining: %s'%(i+1, c_time_str, eta_str)

        c_time = time.time() - start_time
        if c_time > 86400: #more than a day
            c_time_str = '%d days, %d hours, %d minutes, and %.2f seconds.'%(int(c_time)/86400, (int(c_time)/3600)%24, (int(c_time)/60)%60, c_time%60)
        elif c_time > 3600: #more than an hr
            c_time_str = '%d hours, %d minutes, and %.2f seconds.'%((int(c_time)/3600)%24, (int(c_time)/60)%60, c_time%60)
        elif c_time > 60: #more than a min
            c_time_str = '%d minutes and %.2f seconds.'%((int(c_time)/60)%60, c_time%60)
        else:
            c_time_str = '%.2f seconds.'%(c_time)
        print 'Completed Pgen computation for %d sequences: in %s'%(len(seqs), c_time_str)
        outfile.close()


    #Non-safe mode
    else:
        print 'Starting Pgen computation with safe mode disabled.'
        pgens_computed = 0

        infile = open(infile_name, 'r')
        outfile = open(outfile_name, 'w')

        start_time = time.time()
        for i, line in enumerate(infile):
            if comment_delimiter is not None: #Default case -- no comments/header delimiter
                if line.startswith(comment_delimiter): #allow comments
                    continue
            if i < lines_to_skip:
                continue

            if delimiter is None: #Default delimiter is just any whitespace
                split_line = line.split()
            else:
                split_line = line.split(delimiter)

            #Find seq
            try:
                seq = split_line[seq_in_index].strip()
                if skip_empty_sequences and len(seq.strip()) == 0:
                    continue
            except IndexError: #no index match for seq
                if skip_empty_sequences:
                    if len(line.strip()) == 0:
                        continue
                print 'seq_in_index is out of range'
                print 'Exiting...'
                infile.close()
                outfile.close()
                return -1
            #Find and format V_usage_mask
            if V_mask_index is None:
                V_usage_mask = None #default mask
            else:
                try:
                    V_usage_mask = split_line[V_mask_index].strip().split(gene_mask_delimiter)
                except IndexError: #no index match for V_mask_index
                    print 'V_mask_index is out of range'
                    print 'Exiting...'
                    infile.close()
                    outfile.close()
                    return -1

            #Find and format J_usage_mask
            if J_mask_index is None:
                J_usage_mask = None #default mask
            else:
                try:
                    J_usage_mask = split_line[J_mask_index].strip().split(gene_mask_delimiter)
                except IndexError: #no index match for J_mask_index
                    print 'J_mask_index is out of range'
                    print 'Exiting...'
                    infile.close()
                    outfile.close()
                    return -1

            if seq_type_in == 'aaseq':
                aaseq = seq
                #Compute Pgen and print out
                outfile.write(aaseq + delimiter_out + str(pgen_model.compute_aa_CDR3_pgen(aaseq, V_usage_mask, J_usage_mask, print_warnings)) + '\n')
            elif seq_type_in == 'ntseq':
                ntseq = seq
                aaseq = nt2aa(ntseq)
                #Compute Pgen and print out based on recomb_type and seq_type_out
                if seq_type_out == None:
                    outfile.write(ntseq + delimiter_out + str(pgen_model.compute_nt_CDR3_pgen(ntseq, V_usage_mask, J_usage_mask, print_warnings)) + delimiter_out + aaseq + delimiter_out +  str(pgen_model.compute_aa_CDR3_pgen(aaseq, V_usage_mask, J_usage_mask, print_warnings)) + '\n')
                elif seq_type_out == 'ntseq':
                    outfile.write(ntseq + delimiter_out + str(pgen_model.compute_nt_CDR3_pgen(ntseq, V_usage_mask, J_usage_mask, print_warnings)) + '\n')
                elif seq_type_out == 'aaseq':
                    outfile.write(aaseq + delimiter_out + str(pgen_model.compute_aa_CDR3_pgen(aaseq, V_usage_mask, J_usage_mask, print_warnings)) + '\n')

            pgens_computed += 1
            #Print time update
            if (pgens_computed)%seqs_per_time_update == 0 and time_updates:
                c_time = time.time() - start_time
                if c_time > 86400: #more than a day
                    c_time_str = '%d days, %d hours, %d minutes, and %.2f seconds.'%(int(c_time)/86400, (int(c_time)/3600)%24, (int(c_time)/60)%60, c_time%60)
                elif c_time > 3600: #more than an hr
                    c_time_str = '%d hours, %d minutes, and %.2f seconds.'%((int(c_time)/3600)%24, (int(c_time)/60)%60, c_time%60)
                elif c_time > 60: #more than a min
                    c_time_str = '%d minutes and %.2f seconds.'%((int(c_time)/60)%60, c_time%60)
                else:
                    c_time_str = '%.2f seconds.'%(c_time)


                print 'Pgen computed for %d sequences in %s'%(pgens_computed, c_time_str)

            if max_number_of_seqs is not None:
                if pgens_computed >= max_number_of_seqs:
                    break

        c_time = time.time() - start_time
        if c_time > 86400: #more than a day
            c_time_str = '%d days, %d hours, %d minutes, and %.2f seconds.'%(int(c_time)/86400, (int(c_time)/3600)%24, (int(c_time)/60)%60, c_time%60)
        elif c_time > 3600: #more than an hr
            c_time_str = '%d hours, %d minutes, and %.2f seconds.'%((int(c_time)/3600)%24, (int(c_time)/60)%60, c_time%60)
        elif c_time > 60: #more than a min
            c_time_str = '%d minutes and %.2f seconds.'%((int(c_time)/60)%60, c_time%60)
        else:
            c_time_str = '%.2f seconds.'%(c_time)
        print 'Completed Pgen computation for %d sequences in %s'%(pgens_computed, c_time_str)

        infile.close()
        outfile.close()

if __name__ == '__main__': main()
