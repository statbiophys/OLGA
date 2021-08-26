## Synopsis

OLGA (Optimized Likelihood estimate of immunoGlobulin Amino-acid sequences) is a python 2.7/3.6 software developed to compute the generation probability of amino acid and in-frame nucleotide CDR3 sequences from a generative model of V(D)J recombination.

## Motivation

The recent ubiquity in adaptive immune system repertoire sequencing has led to the need for a variety of probabilistic and bioinformatic tools to analyze CDR3 sequence. Of particular interest are the generative models of V(D)J recombination and an inference procedure first introduced in [Murugan 2012](http://www.pnas.org/content/early/2012/09/10/1212755109.short). Recently, an efficient software package [IGoR](https://github.com/qmarcou/IGoR) (Inference and Generation Of Repertoires) was published [(Marcou 2018)](http://www.nature.com/articles/s41467-018-02832-w) and is available on [GitHub](https://github.com/qmarcou/IGoR) and should be used for any generative model inference.

This software implements the dynamic programming algorithm discussed in [Sethna 2018](https://arxiv.org/abs/1807.04425) to compute the generation probability of amino acid and in-frame nucleotide CDR3 sequences.

As adaptive repertoire researchers span the full gamut of technical coding skills, the goal of this README and software is to be as usable and useful as possible for researchers with any programming background. To this end, the core algorithm provided as Python modules with a few command line console scripts provided. The technical user can thus import the Python modules and incorporate them into their own analysis pipeline while the researcher who only requires a cursory analysis can use the provided console scripts. (Furthermore, there is a functional/procedural implementation of the algorithm -- please contact via email if you prefer the modules structured in this manner).

Documentation and examples provided below. Ideally, even the interested party who is uncomfortable with coding should be able to manage from the examples.

## Version
Latest released version: 1.2.4

## Installation
OLGA is a python 2.7 software which only uses standard python libraries and requires no additional dependencies. It is available on PyPI and can be downloaded and installed through pip: ```pip install olga```.

OLGA is also available on [GitHub](https://github.com/zsethna/OLGA). The command line entry points can be installed by using the setup.py script: ```$ python setup.py install```. If the command line console scripts are not wanted, no installation is necessary and the scripts ```compute_pgen.py``` and ```generate_sequences.py``` can be run as executables.



Directory architecture:
```
olga/
│   README.md
│   LICENSE
│   setup.py
│   example_expanded_amino_acid_alphabet.txt  
│
└───olga/
    │   __init__.py
    │   compute_pgen.py
    │   generate_sequences.py
    │   generation_probability.py
    │   preprocess_generative_model_and_data.py
    │   load_model.py
    │   sequence_generation.py
    │
    └───default_models/
        └───human_T_alpha/
        │       model_marginals.txt
        │       model_params.txt
        │       J_gene_CDR3_anchors.txt
        │       V_gene_CDR3_anchors.txt
        │
        └───human_T_beta/
        │       model_marginals.txt
        │       model_params.txt
        │       J_gene_CDR3_anchors.txt
        │       V_gene_CDR3_anchors.txt
        │
        └───human_B_heavy/
        │       model_marginals.txt
        │       model_params.txt
        │       J_gene_CDR3_anchors.txt
        │       V_gene_CDR3_anchors.txt
        │
        └───mouse_T_beta/
                model_marginals.txt
                model_params.txt
                J_gene_CDR3_anchors.txt
                V_gene_CDR3_anchors.txt
```

## Command line console scripts and Examples

There are two command line console scripts (the scripts can still be called as executables if OLGA is not installed):
1. olga-compute_pgen
  * computes the generation probability CDR3 sequences according to a generative V(D)J model
2. olga-generate_sequences
  * generates CDR3 sequences from a generative V(D)J model

For any of them you can execute with the -h or --help flags to get the options.

### Quick Demo
After installing OLGA, we offer a quick demonstration of the console scripts. This will demonstrate generating sequences and computing pgens from the default model for human TCR beta chains that ships with OLGA. This demo will generate two files, example_seqs.tsv and example_pgens.tsv, in the directory that these command line calls are made.

1. ```$ olga-compute_pgen --humanTRB CASSLGRDGGHEQYF```
  * This computes the pgen of the amino acid sequence CASSLGRDGGHEQYF (you should get ~7.25342176315e-10)

2. ```$ olga-generate_sequences --humanTRB -n 5```
  * Generate 5 human TRB CDR3 sequences (both nucleotide and amino acid sequences) and print to stdout along with the V and J genes used to generate them.

3. ```$ olga-generate_sequences --humanTRB -o example_seqs.tsv -n 1e2```
  * This generates a file example_seqs.tsv and writes 100 generated human TRB CDR3 sequences.

4. ```$ olga-compute_pgen --humanTRB -i example_seqs.tsv -m 5```
  * This reads in the first 5 sequences from the file we just generated, example_seqs.tsv, and and computes both the nucleotide sequence and the amino acid sequence pgens for the 5 sequences and prints it to stdout.

5. ```$ olga-compute_pgen --humanTRB -i example_seqs.tsv -o example_pgens.tsv```
  * This reads in the full file example_seqs.tsv, and computes both the nucleotide sequence and the amino acid sequence pgen for each of the 100 sequences and writes them to the file example_pgens.tsv. A running display will be printed to stdout with the last few lines computed along with time/rate updates.

### Specifying a default V(D)J model (or a custom model folder)
All of the console scripts require specifying a V(D)J generative model and genomic data. OLGA ships with 4 default models that can be indicated by flags, or a custom model folder can be indicated.

| Options                                        | Description                                      |
|------------------------------------------------|--------------------------------------------------|
| **--humanTRA**                                 | Default human T cell alpha chain model (VJ)      |
| **--humanTRB**                                 | Default human T cell beta chain model (VDJ)      |
| **--mouseTRB**                                 | Default mouse T cell beta chain model (VDJ)      |
| **--humanIGH**                                 | Default human B cell heavy chain model (VDJ)     |
| **--set_custom_model_VJ** PATH/TO/MODEL_FOLDER/ | Specifies the directory PATH/TO/MODEL_FOLDER/ of a custom VJ generative model      |
| **--set_custom_model_VDJ** PATH/TO/MODEL_FOLDER/| Specifies the directory PATH/TO/MODEL_FOLDER/ of a custom VDJ generative model     |

Note, if specifying a folder for a custom VJ recombination model
(e.g. an alpha or light chain model) or a custom VDJ recombination model
(e.g. a beta or heavy chain model), the folder must contain the following files
with the exact naming convention:

* model_params.txt ([IGoR](https://github.com/qmarcou/IGoR) inference param file)
* model_marginals.txt ([IGoR](https://github.com/qmarcou/IGoR) inference marginal file)
* V_gene_CDR3_anchors.csv (V anchor residue position and functionality file)
* J_gene_CDR3_anchors.csv (J anchor residue position and functionality file)

The console scripts can only read files of the assumed anchor.csv/[IGoR](https://github.com/qmarcou/IGoR) syntaxes. If the advanced user wants to read in models from files of other formats please read the discussion in the Python module section and the documentation of load_model.py.

### compute_pgen.py (olga-compute_pgen)

This script will compute  the generation probabilities (pgens) of sequences, as defined by a specified generative model. The sequences must be TRIMMED TO ONLY THE CDR3 region as defined by the V and J anchor files (default is to INCLUDE the conserved residues of the C in the V region and the F/W in the J region).

Each sequence will be determined to be one of:
* 'Amino acid' sequence including any ambiguous amino acid symbols (see
Options for alphabet_file to specify custom symbols)
* Nucleotide sequence (in-frame)
* Regular expression template for 'amino acid sequences'

This script can read in sequences and output pgens in one of three modes (determined automatically based on the input arguments):

1. Pass CDR3 sequences as arguments, output is printed to stdout.
2. Read in CDR3 sequences from a file, output is printed to stdout.
3. Read in CDR3 sequences from a file, output to a file, dynamic display with time updates printed to stdout.

Full options can be printed by: ```$ olga-compute_pgen -h```

**Mode 1):**

This mode is provided as a way to quickly compute pgen of just a single or couple of sequences. It is not advisable to use this method to set up a loop over a lot of sequences as each call of the script demands the overhead of processing a model. To compute the pgens of many sequences, it is suggested to read the sequences in from a file and use either mode 2 or 3.

It is also possible to condition the Pgen computation on V and J identity by specifying the V or J usages as a mask. However, note that these V/J masks will be applied to ALL of the sequences provided as arguments. Read Options on v_mask and j_mask for more info.


Example calls:
* Compute the pgen of an amino acid CDR3 sequence
```
$ olga-compute_pgen --humanTRB CASSTGQANYGYTF
Pgen of the amino acid sequence CASSTGQANYGYTF: 5.26507446955e-08
```
* Compute the pgen of an in-frame nucleotide sequence *and* the amino acid sequence it translates to.
```
$ olga-compute_pgen --humanTRB TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT
Pgen of the nucleotide sequence TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT: 1.31873701121e-17
Pgen of the amino acid sequence nt2aa(TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT) = CASSDAQGRNRGTEAFF: 4.70599549953e-13
```
* Compute the pgen of a regular expression template of CDR3 amino acid sequences. Note, for a regular expression sequence, provided as an argument, backslashes may be needed to specify the characters {} for the sequence to be read in properly.
```
$ olga-compute_pgen --humanTRB CASSTGX\{1,5\}QAN[YA]GYTF
Pgen of the regular expression sequence CASSTGX{1,5}QAN[YA]GYTF: 7.588241802e-08
```
* Compute the pgens of all three sequences.
```
$ olga-compute_pgen --humanTRB CASSTGQANYGYTF CASSTGX\{1,5\}QAN[YA]GYTF TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT
Pgen of the amino acid sequence CASSTGQANYGYTF: 5.26507446955e-08
Pgen of the regular expression sequence CASSTGX{1,5}QAN[YA]GYTF: 7.588241802e-08
Pgen of the nucleotide sequence TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT: 1.31873701121e-17
Pgen of the amino acid sequence nt2aa(TGTGCCAGCAGTGACGCACAGGGGCGTAATCGTGGGACTGAAGCTTTCTTT) = CASSDAQGRNRGTEAFF: 4.70599549953e-13
```
* Specify a comma delimited V or J mask to condition the pgen computation on V and/or J gene usage.
```
$ olga-compute_pgen --humanTRB CASSTGQANYGYTF --v_mask TRBV2,TRBV14 --j_mask TRBJ1-2
Pgen of the amino acid sequence CASSTGQANYGYTF: 1.39165562898e-09
```

It is also possible to restrict the Pgen computation to specified V and/or J genes or alleles (to reflect any alignment outside of the CDR3 region) by using the options -v or -j (see example below). You can specify multiple V or J genes/alleles by using a comma as a delimiter.

The only required inputs are the sequence and specifying the generative V(D)J model. Additional options can be found by using -h.

Modes 2/3):

These read in sequences from a file. The script has only minimal file parsing built in, so reading in sequences from a file requires the file to be structured with delimiter spaced values (i.e. the  data is organized in columns separated by delimiter like a .tsv or .csv file). Read Options on delimiter for more info.

To read in sequences, the index of column of CDR3 sequences is needed. The default is to assume that the sequences to be read in are in the first column (index 0), meaning that a text file with only a sequence on each line will be read in okay by default. Read Options on seq_in for more info.

It is not recommended to read in regular expression sequences from a file. These sequences require enumerating out the amino acid sequences which correspond to them and computing pgen for each of them individually -- this can require a large time cost. Instead consider defining a custom 'amino acid' alphabet to define the symbols used in the regular expressions if possible. Furthermore, BE CAREFUL if reading in from a .csv file -- if commas are used in a regex sequence and comma is used as the delimiter of the .csv file, the sequence will not be read in properly.

If nucleotide sequences are to be read in it is possible to specify if the
output should be the nucleotide sequence Pgen and/or the translated amino acid sequence Pgen (the default is to compute and output both). See Options.

It is also possible to condition the Pgen computation on V and J identity by specifying what index the column that V and J masks are stored for each line.

Mode 2 does not have a specified output file and so will print the sequences and their pgens to stdout.

Mode 3 does have a specified output file. By default in this mode there is a running display of the last few sequences/pgens written to the output file as well as time elapsed, current rate of computation, and estimated time remaining. This display can be disabled (see Options).

As it is rare for datasets to have >> 1e4 unique sequences, parallelization is not built in to compute_pgen. However, there are options to skip N lines of the file and to load at most M sequences so, if wanted, one could build a parallelized wrapper around this script (though it would be recommended to instead just import the modules and build from there).

Example calls (assuming a file example_seqs.tsv with the line structure
ntseq   aaseq   V_mask  J_mask):


* Read in the ntseqs and print the ntseq, aaseq = nt2aa(ntseq) and their pgens to stdout
```
$ olga-compute_pgen -i example_seqs.tsv --humanTRB
```

* Only read in the first 10 sequences
```
$ olga-compute_pgen -i example_seqs.tsv --humanTRB -m 10
```

* Read in the ntseqs, write the ntseq, aaseq = nt2aa(ntseq), and their pgens to example_pgens.tsv
```
$ olga-compute_pgen -i example_seqs.tsv --humanTRB -o example_pgens.tsv
```

* Specify the V/J mask indices
```
$ olga-compute_pgen -i example_seqs.tsv --humanTRB -o example_pgens.tsv --v_in 2 --j_in 3
```

* Read in the aaseq column instead of the ntseq column
```
$ olga-compute_pgen -i example_seqs.tsv --humanTRB -o example_pgens.tsv --seq_in 1
```

| Selected Options                               | Description                                      |
|------------------------------------------------|--------------------------------------------------|
|   **-h**, **--help**                           |   show full Options list and exit                          |
|   **-i** PATH/TO/FILE                          |   read in CDR3 sequences (and optionally V/J masks) from PATH/TO/FILE |
|   **-o** PATH/TO/FILE                          |   write CDR3 sequences and pgens to PATH/TO/FILE |
|   **--seq_in** INDEX                           |   specifies sequences to be read in are in column INDEX. Default is index 0 (the first column).    |
|   **--v_in** INDEX                             |   specifies V_masks are found in column INDEX in the input file. Default is no V mask.   |
|   **--j_in** INDEX                             |   specifies J_masks are found in column INDEX in the input file. Default is no J mask.   |
|   **--v_mask** V_MASK                          |   specify V usage to condition Pgen on for seqs read in as arguments.   |
|   **--j_mask** J_MASK                          |   specify J usage to condition Pgen on for seqs read in as arguments.   |
|   **-m** N                                     |   compute Pgens for at most N sequences.         |
|   **--lines_to_skip** N                        |   skip the first N lines of the file. Default is 0. |
|   **-a** PATH/TO/FILE                          |  specify PATH/TO/FILE defining a custom 'amino acid' alphabet. Default is no custom alphabet.     |
|   **--seq_type_out** SEQ_TYPE                  |   if read in sequences are ntseqs, declare what type of sequence to compute pgen for. Default is all. Choices: 'all', 'ntseq', 'aaseq'    |
|  **--display_off**                             |   turn the sequence display off (only applies in write-to-file mode). Default is on.    |
|  **-d** DELIMITER                              |   declare infile delimiter. Default is tab for .tsv input files, comma for .csv files, and any whitespace for all others. Choices: 'tab', 'space', ',', ';', ':'    |
|  **--raw_delimiter** DELIMITER                 |   declare infile delimiter as a raw string.    |
|  **--delimiter_out**,  **--raw_delimiter_out**,  **--gene_mask_delimiter**, **--raw_gene_mask_delimiter**      |   declare delimiters for the outfile and for gene masks (read in from the columns of v_mask_index and j_mask_index). Same syntax as the infile delimiter.  |
|   **--comment_delimiter** COMMENT_DELIMITER         |  character or string to indicate comment or header lines to skip.     |

### generate_sequences.py (olga-generate_sequences)

This program will generate a file of Monte Carlo sampling from a specified generative V(D)J model. The sequences generated will have NO ERRORS.

It is required to specify the number of sequences to be generated. This is done with -n (see Options).

If a file is specified to write to (using -o, see Options), the generated sequences are written to the file, otherwise they are printed to stdout.

The default is to record both the nucleotide CDR3 sequence and the amino acid CDR3 sequence. This can be specified (see Options).

The V/J genes used to generate each sequence can be recorded or not. Default is to record them, but this can be toggled off with --record_genes_off (see Options).

Full options can be printed by: ```$ olga-generate_sequences -h```

Example calls:

* Print 20 generated sequences to stdout
```
$ olga-generate_sequences --humanTRB -n 20
```

* Write 200 generated sequences to example_seqs.tsv
```
$ olga-generate_sequences --humanTRB -o example_seqs.tsv -n 200
```

* Write 20,000 generated sequences to example_seqs.tsv
```
$ olga-generate_sequences --humanTRB -o example_seqs.tsv -n 2e4
```

* Write only the amino acid sequences
```
$ olga-generate_sequences --humanTRB -o example_seqs.tsv -n 200 --seq_type amino_acid
```

| Selected Options                               | Description                                      |
|------------------------------------------------|--------------------------------------------------|
|   **-h**, **--help**                           |   show full Options list and exit                |
|   **-o** PATH/TO/FILE                          |   write CDR3 sequences to PATH/TO/FILE           |
|   **-n** N                                     |   specify the number of sequences to generate.   |
|   **--seed** SEED                              |   set seed for pseudorandom number generator. Default is to not set a seed. |
|   **--seq_type** SEQ_TYPE                      |   declare sequence type for output sequences. Choices: 'all' [default], 'ntseq', 'aaseq' |
|   **--time_updates_off**                       |   turn time updates off.
|   **--record_genes_off**                       |   turn off recording V and J gene info.
|   **-d** DELIMITER                             |   declare delimiter choice. Default is tab for .tsv output files, comma for .csv files, and tab for all others. Choices: 'tab', 'space', ',', ';', ':' |
|   **--raw_delimiter** DELIMITER                |   declare delimiter choice as a raw string.      |

## Using the OLGA modules in a Python script (advanced users)
In order to incorporate the core algorithm into an analysis pipeline (or to write your own script wrappers) all that is needed is to import the modules and load up a generative model. Each module defines some classes that only a few methods get called on.

As the generative model structure is different between VDJ recombination and VJ recombination the algorithms to compute Pgen for the two are different. For this reason, different objects are defined for VDJ recombination models and VJ recombination models, however the *methods* that get called are the same.

The modules are:

| Module name                                    | Classes                                          |    
|------------------------------------------------|--------------------------------------------------|
| load_model.py                                  | GenomicDataVDJ, GenomicDataVJ, GenerativeModelVDJ, GenerativeModelVJ|
| preprocess_generative_model_and_data.py        | PreprocessedParametersVDJ, PreProcessedParametersVJ|
| generation_probability.py                      | GenerationProbabilityVDJ, GenerationProbabilityVJ|
| sequence_generation.py                         | SequenceGenerationVDJ, SequenceGenerationVJ      |
| utils.py                                       | N/A (contains util functions)                    |

The classes with methods that are of interest will be GenerationProbabilityV(D)J (to compute Pgens) and SequenceGenerationV(D)J (to generate sequences).

There is a fair amount of parameter processing that must go on to call these methods, however this is generally all done by instantiating a particular class. An exception to this rule are the classes GenerativeModelV(D)J and GenomicDataV(D)J. Normally the genomic data and model parameters are read in from [IGoR](https://github.com/qmarcou/IGoR) inference files (and prepared V and J anchor files that have been prepared), however this is not mandated in order to make it easier for people to adapt the code to read in models/genomic data from other sources.

Instantiating GenerativeModelV(D)J and GenomicDataV(D)J leaves the attributes as dummies, and calling the methods load_and_process_igor_model and load_igor_genomic_data will load up [IGoR](https://github.com/qmarcou/IGoR) files.

If you want to load models/data from other sources, you will need to write your own methods to set the attributes in GenerativeModelV(D)J and GenomicDataV(D)J. Please see the documentation of load_model.py for more details.

Here is an example of loading the default human TRB model to compute some sequence Pgens and to generate some random CDR3 sequences:

```
>>> import olga.load_model as load_model
>>> import olga.generation_probability as pgen
>>> import olga.sequence_generation as seq_gen
>>>
>>> #Define the files for loading in generative model/data
... params_file_name = 'default_models/human_T_beta/model_params.txt'
>>> marginals_file_name = 'default_models/human_T_beta/model_marginals.txt'
>>> V_anchor_pos_file ='default_models/human_T_beta/V_gene_CDR3_anchors.csv'
>>> J_anchor_pos_file = 'default_models/human_T_beta/J_gene_CDR3_anchors.csv'
>>>
>>> #Load data
... genomic_data = load_model.GenomicDataVDJ()
>>> genomic_data.load_igor_genomic_data(params_file_name, V_anchor_pos_file, J_anchor_pos_file)
>>> #Load model
... generative_model = load_model.GenerativeModelVDJ()
>>> generative_model.load_and_process_igor_model(marginals_file_name)
>>>
>>> #Process model/data for pgen computation by instantiating GenerationProbabilityVDJ
... pgen_model = pgen.GenerationProbabilityVDJ(generative_model, genomic_data)
>>>
>>> #Compute some sequence pgens
... pgen_model.compute_regex_CDR3_template_pgen('CASSAX{0,5}SARPEQFF')
6.846877804096558e-10
>>> pgen_model.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF', 'TRBV30*01', 'TRBJ1-2*01')
1.203646865765782e-10
>>> pgen_model.compute_nt_CDR3_pgen('TGTGCCAGTAGTATAACAACCCAGGGCTTGTACGAGCAGTACTTC')
3.9945642868171824e-14
>>>
>>>
>>> #Process model/data for sequence generation by instantiating SequenceGenerationVDJ
... seq_gen_model = seq_gen.SequenceGenerationVDJ(generative_model, genomic_data)
>>>
>>> #Generate some random sequences
... seq_gen_model.gen_rnd_prod_CDR3()
('TGTGCCAGCAGTGAAAAAAGGCAATGGGAAAGCGGGGAGCTGTTTTTT', 'CASSEKRQWESGELFF', 27, 8)
>>> seq_gen_model.gen_rnd_prod_CDR3()
('TGTGCCAGCAGTTTAGTGGGAAGGGCGGGGCCCTATGGCTACACCTTC', 'CASSLVGRAGPYGYTF', 14, 1)
>>> seq_gen_model.gen_rnd_prod_CDR3()
('TGTGCCAGCTGGACAGGGGGCAACTACGAGCAGTACTTC', 'CASWTGGNYEQYF', 55, 13)
```
Additional documentation of the modules is found in their docstrings (accessible either through pydoc or using help() within the python interpreter).

## Notes about CDR3 sequence definition

This code is quite flexible, however it does demand a very consistent definition of CDR3 sequences.

**CHECK THE DEFINITION OF THE CDR3 REGION OF THE SEQUENCES YOU INPUT.** This will likely be the most often problem that occurs.

The default models/genomic data are set up to define the CDR3 region from the conserved cysteine C (INCLUSIVE) in the V region to the conserved F or W (INCLUSIVE) in the J. This corresponds to positions X and X according to IMGT. This can be changed by altering the anchor position files, however the user is strongly recommended against this. If a sequence which does NOT match the definition of this definition of CDR3 region it will likely output a generation probability of 0.

One of the most powerful aspects of the algorithm is the ability to define a custom 'amino acid' alphabet. By this we mean an 'amino acid' symbol is a character that corresponds to some collection of codons. A CDR3 sequence that is fed into the algorithm is a string of such characters. Indeed, the way the algorithm runs on nucleotide sequences is by defining symbols for all 64 codons. The default 'amino acid' alphabet includes all standard amino acids, symbols for the 64 codons, and the standard ambiguous amino acids B, J, X, and Z.

It is possible to add custom symbols to this alphabet. For instance it may be useful to have symbols for collections of amino acids with certain electrostatic/steric properties or (if interested in mutations/errors) the collection of codons hamming distance 1 (in nucleotides) from a particular amino acid.

To facilitate adding to the base alphabet files of a particular syntax can be specified and read in. Generally the syntax of the file is:

```
symbolA: list,of,amino,acids,or,codons,separated,by,commas
symbolB: another,comma,delimited,list
```
etc, where symbolA and symbolB are single characters and aren't one of the protected characters. Please see utils.py for more documentation, and the options of the console scripts to see how to call such alphabet files.

In addition to allowing customization of the 'amino acid' alphabet, we include functions that can parse regular expressions with a limited vocabulary. In particular only [] and {} are supported (the symbol for a Kleene star, \*, must be reserved for stop codons). The sequences corresponding to the regular expression can then be enumerated and their Pgens summed. Note, this can be slow as Pgen must be computed for each of the enumerated sequences independently. If a particular combination of amino acids is being used in a [] frequently you may consider defining a symbol for that combination and adding it to the alphabet. For example the regular expression 'CASS**[ACDEFGHIKLMNPQRSTVWY]**SARPEQFF' will list out 20 sequences, but its Pgen could be computed by considering the single CDR3 sequence 'CASS**X**SARPEQFF' Please see documentation in pgen.py for more info and examples. To prevent overflow, the arbitrary length repetition {x,} is cutoff at a maximum of 40 (if you want a higher number {x,y} will go higher than 40). Furthermore, if more than 10,000 sequences correspond to the template, the regular expression will be rejected.

## Contact

Any issues or questions should be addressed to [us](mailto:zachary.sethna@gmail.com).

## License

Free use of OLGA is granted under the terms of the GNU General Public License version 3 (GPLv3).
