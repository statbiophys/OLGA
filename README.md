## Synopsis

OLGA is a python 2.7 software developed to compute the generation probability of amino acid and in-frame nucleotide CDR3 sequences from a generative model of V(D)J recombination.

## Motivation

The recent ubiquity in adaptive immune system repertoire sequencing has led to the need for a variety of probabilistic and bioinformatic tools to analyze CDR3 sequence. Of particular interest are the generative models of V(D)J recombination and an inference procedure first introduced in [Murugan 2012](http://www.pnas.org/content/early/2012/09/10/1212755109.short). Recently, a more complete and efficient software package IGoR (Inference and Generation Of Repertoires) was published [(Marcou 2018)](http://www.nature.com/articles/s41467-018-02832-w) and is available on GitHub [(IGoR)](https://github.com/qmarcou/IGoR) and should be used for any generative model inference.

This software implements the dynamic programming algorithm discussed *CITE PAPER* to compute the generation probability of amino acid and in-frame nucleotide CDR3 sequences.

As adaptive repertoire researchers span the full gamut of technical coding skills, the goal of this README and software is to be as usable and useful as possible for researchers with any programming background. To this end, the core algorithm provided as Python modules with a few command line console scripts provided. The technical user can thus import the Python modules and incorporate them into their own analysis pipeline while the researcher who only requires a cursory analysis can use the provided console scripts. (Furthermore, there is a functional/procedural implementation of the algorithm -- please contact via email if you prefer the modules structured in this manner).

Documentation and examples provided below, and assumes only minimal familiarity with Python/Bash. Ideally, even the interested party who is very uncomfortable with coding should be able to manage from the examples.

## Version
Latest released version: 0.1.0

## Installation
OLGA is a python 2.7 software which only uses standard python libraries and requires no additional dependencies. It is available on PyPI and can be downloaded and installed through pip: ```pip -install olga```.

OLGA is also available on [GitHub](https://github.com/zsethna/OLGA). The command line entry points can be installed by using the setup.py script: ```$ python setup.py install```. If the command line console scripts are not wanted, no installation is necessary and the scripts ```run_pgen.py```, ```compute_single_sequence_pgen.py```, and ```generate_synthetic_sequences.py``` can all be run as executables.



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
    │   run_pgen.py
    │   compute_single_sequence_pgen.py
    │   generate_synthetic_sequences.py
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

There are three command line console scripts (the scripts can still be called as executables if OLGA is not installed):
1. olga-compute_single_sequence_pgen
  * computes the generation probability of a single sequence.
2. olga-run_pgen
  * computes the generation probability of sequences read from a file and outputs the results to another file.
3. olga-generate_synthetic_sequences
  * generates sequences and writes them to a file.

For any of them you can execute with the -h or --help flags to get the options.

### Quick Demo
After installing OLGA, we offer a quick demonstration of the console scripts. This will demonstrate generating sequences and computing pgens from the default model for human TCR beta chains that ships with OLGA. This demo will generate two files, example_seqs.tsv and example_pgens.tsv, in the directory that these command line calls are made.

1. ```$ olga-compute_single_sequence_pgen CASSLGRDGGHEQYF --humanTCRB```
  * This computes the pgen of the amino acid sequence CASSLGRDGGHEQYF (you should get 7.25342176315e-10)

2. ```$ olga-generate_synthetic_sequences.py example_seqs.tsv --humanTCRB -n 1e3```
  * This generates a file example_seqs.tsv with 1000 human TCRB CDR3 sequences (both the nucleotide sequence and the amino acid sequence) along with the V and J genes used to generate them.

3. ```$ olga-run_pgen example_seqs.tsv example_pgens.tsv --humanTCRB```
  * This reads in the file we just generated, example_seqs.tsv, and computes both the nucleotide sequence and the amino acid sequence pgen for each of the 1000 sequences we generated and writes it to the file example_pgens.tsv.

Did it work? If not, check these issues:
* Do you have Python 2.7? (Should be built into the OS for Unix/Mac)
* Is OLGA installed?
  * If OLGA is not 'installed,' the scripts can still be called as executables (or through python):
    1. ```$ ./compute_single_sequence_pgen.py CASSLGRDGGHEQYF --humanTCRB```
    2. ```$ ./generate_synthetic_sequences.py example_seqs.tsv --humanTCRB -n 1e3```
    3. ```$ ./run_pgen.py example_seqs.tsv example_pgens.tsv --humanTCRB```

### Specifying a default V(D)J model (or a custom model folder)
All of the console scripts require specifying a V(D)J generative model and genomic data. OLGA ships with 4 default models that can be indicated by flags, or a custom model folder can be indicated.

| Options                                        | Description                                      |
|------------------------------------------------|--------------------------------------------------|
| **--humanTCRA**, **--human_T_alpha**           | Default human T cell alpha chain model (VJ)      |
| **--humanTCRB**, **--human_T_beta**            | Default human T cell beta chain model (VDJ)      |
| **--mouseTCRB**, **--mouse_T_beta**            | Default mouse T cell beta chain model (VDJ)      |
| **--humanIGH**, **--human_B_heavy**            | Default human B cell heavy chain model (VDJ)     |
| **-VDJ_model_folder** PATH/TO/MODEL_FOLDER/    | Specifies a custom VDJ model directory           |
| **-VJ_model_folder** PATH/TO/MODEL_FOLDER/     | Specifies a custom VDJ model directory           |

Note, if specifying a custom model folder for either a VJ recombination model
(e.g. an alpha or light chain model) or a VDJ recombination model
(e.g. a beta or heavy chain model), the folder must contain the following files
with the exact naming convention:

* model_params.txt (IGoR inference param file)
* model_marginals.txt (IGoR inference marginal file)
* V_gene_CDR3_anchors.csv (V residue anchor and functionality file)
* J_gene_CDR3_anchors.csv (J residue anchor and functionality file)

The console scripts can only read files of the assumed IGoR/anchor.csv syntaxes. In order to read in models from files of other formats, please read the discussion in the Python module section and the documentation of load_model.py.

### compute_single_sequence_pgen.py

syntax: ```olga-compute_single_sequence_pgen CDR3_SEQ **kwarg ```

This console script is used to compute the generation probability (Pgen) of a single CDR3 sequence as defined by a provided generative V(D)J model. The sequence itself can be an 'amino acid' sequence, an in-frame nucleotide sequence, or an 'amino acid' regular expression.

The CDR3 sequence is the first argument. The program has a minimal sequence parser that will guess if the inputted sequence is an 'amino acid' sequence, a nucleotide sequence, or a regular expression and the printed output will indicate what the parser interpreted the sequence to be.

It is also possible to restrict the Pgen computation to specified V and/or J genes or alleles (to reflect any alignment outside of the CDR3 region) by using the options -v or -j (see example below). You can specify multiple V or J genes/alleles by using a comma as a delimiter.

The only required inputs are the sequence and specifying the generative V(D)J model. Additional options can be found by using -h.

Examples:
```
$ olga-compute_single_sequence_pgen CASSLGRDGGHEQYF --humanTCRB

------------------------------------------------------------------------------------------
Pgen of the amino acid sequence CASSLGRDGGHEQYF: 7.25342176315e-10
------------------------------------------------------------------------------------------


$ olga-compute_single_sequence_pgen TGTGCCAGCAGCTTAGGTAGGGATGGAGGTCACGAGCAGTACTTC --humanTCRB

------------------------------------------------------------------------------------------
Pgen of the nucleotide sequence TGTGCCAGCAGCTTAGGTAGGGATGGAGGTCACGAGCAGTACTTC: 8.50807296092e-14
Pgen of the amino acid sequence nt2aa(TGTGCCAGCAGCTTAGGTAGGGATGGAGGTCACGAGCAGTACTTC) = CASSLGRDGGHEQYF: 7.25342176315e-10
------------------------------------------------------------------------------------------


$ olga-compute_single_sequence_pgen CASSLX\{0,5\}DG[GAR]HEQYF --humanTCRB

------------------------------------------------------------------------------------------
Pgen of the regular expression sequence CASSLX{0,5}DG[GAR]HEQYF: 2.42813803819e-07
------------------------------------------------------------------------------------------


$ olga-compute_single_sequence_pgen.py CASSLGRDGGHEQYF --humanTCRB -v TRBV11-1 -j TRBJ2-7

------------------------------------------------------------------------------------------
Pgen of the amino acid sequence CASSLGRDGGHEQYF: 6.66613706172e-12

(Conditioned on the V and J gene/allele usages: ['TRBV11-1'], ['TRBJ2-7'])
------------------------------------------------------------------------------------------

$ olga-compute_single_sequence_pgen CASSLGRDGGHEQYF --humanTCRB -v TRBV2,TRBV11-1 -j TRBJ2-7,not_a_J_gene

Unfamiliar J gene/allele: not_a_J_gene
------------------------------------------------------------------------------------------
Pgen of the amino acid sequence CASSLGRDGGHEQYF: 7.93830612446e-12

(Conditioned on the V and J gene/allele usages: ['TRBV2', 'TRBV11-1'], ['TRBJ2-7'])
------------------------------------------------------------------------------------------

```


### run_pgen.py

syntax: ```olga-run_pgen PATH/TO/INFILE PATH/TO/OUTFILE **kwarg ```

This console script is used to compute the generation probabilities (Pgens), as defined by a specified generative V(D)J model, of CDR3 sequences from a file and writes the output to another file. The infile is assumed to be a DELIMITER SPACED VALUE file (e.g. a tab separated file .tsv or a comma separated file .csv, though other delimiters can be specified. For more info read the options/docstring).

The CDR3 sequences can be either 'amino acid' sequences or in-frame nucleotide sequences and are read in as a specific column of each line (as defined by the delimiter). The default index for the sequence column is 0 (the first column) to ensure that a file composed of only a single sequence on each line will be correctly read in by default.

It is also possible to restrict the Pgen computation to specified V and/or J genes or alleles to reflect any alignment outside of the CDR3 region by specifying the index of columns that contain the V/J mask to be read in (syntax of the mask string is similar to what was used for compute_single_sequence_pgen.py - additional information about this syntax in the options or docstring).

Required inputs are the PATH/TO/INFILE and PATH/TO/OUTFILE (which are the first two arguments), and specifying the generative V(D)J model.


Further documentation of additional options along with more examples can be found in the options/docstring.


Examples:
```
$ olga-run_pgen example_seqs.tsv example_pgen.tsv --humanTCRB

$ olga-run_pgen example_seqs.tsv example_pgen.tsv --humanTCRB -seq_in_index 1

```

### generate_synthetic_sequences.py

syntax: ```olga-generate_synthetic_sequences PATH/TO/OUTFILE **kwarg```

This console script generates a file of CDR3 sequences generated from Monte Carlo sampling of a specified generative V(D)J model.

Required inputs are the PATH/TO/OUTFILE (assumed to be the first argument), the number of sequences to generate, and specifying the generative V(D)J model.

Further documentation of all options along with more examples can be found in the options/docstring.

Example:
```
$ olga-generate_synthetic_sequences data/example_seqs.tsv --humanTCRB -n 1e3

```

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

There is a fair amount of parameter processing that must go on to call these methods, however this is generally all done by instantiating a particular class. An exception to this rule are the classes GenerativeModelV(D)J and GenomicDataV(D)J. Normally the genomic data and model parameters are read in from IGoR inference files (and prepared V and J anchor files that have been prepared), however this is not mandated in order to make it easier for people to adapt the code to read in models/genomic data from other sources.

Instantiating GenerativeModelV(D)J and GenomicDataV(D)J leaves the attributes as dummies, and calling the methods load_and_process_igor_model and load_igor_genomic_data will load up IGoR files.

If you want to load models/data from other sources, you will need to write your own methods to set the attributes in GenerativeModelV(D)J and GenomicDataV(D)J. Please see the documentation of load_model.py for more details.

Here is an example of loading the default human TCRB model to compute some sequence Pgens and to generate some random CDR3 sequences:

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

In addition to allowing customization of the 'amino acid' alphabet, we include functions that can parse regular expressions with a limited vocabulary. In particular only [] and {} are supported (the symbol for a Kleene star, \*, must be reserved for stop codons). The sequences corresponding to the regular expression can then be enumerated and their Pgens summed. Note, this can be slow as Pgen must be computed for each of the enumerated sequences independently. If a particular combination of amino acids is being used in a [] frequently you may consider defining a symbol for that combination and adding it to the alphabet. For example the regular expression 'CASS**[ACDEFGHIKLMNPQRSTVWY]**SARPEQFF' will list out 20 sequences, but its Pgen could be computed by considering the single CDR3 sequence 'CASS**X**SARPEQFF' Please see documentation in pgen.py for more info and examples.

## Contact

Any issues or questions should be addressed to [us](mailto:sethna@princeton.edu).

## License

Free use of OLGA is granted under the terms of the GNU General Public License version 3 (GPLv3).
