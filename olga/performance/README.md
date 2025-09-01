# Performance submodule

This submodule provides an re-implementation of `compute_CDR3_pgen` with a MUCH faster
Pgen calculation, thanks to `numba`. 

Comes with an updated CLI flag for `olga-compute_pgen` called `--fast_pgen`,
as well as a class called `FastPgen` for developers. This is also useful for
[SONIA](https://github.com/statbiophys/SONIA) and [soNNia](https://github.com/statbiophys/soNNia) users.

> [!WARNING]
> This module is new and was only tested using the `compute_aa_CDR3_pgen` function with V and J genes. 
It is highly recommended you conduct your own tests before using in production.

## Benchmarks

The new method was tested on human AIRR-seq datasets for all available models, each with 500,000 sequences, 
providing the `junction_aa`, `v_call`, and `j_call` as input, with and without the new `--fast_pgen` flag.

| Model (# of seqs) | OLGA | OLGA `--fast_pgen`  | Speedup |
| ----- | ---- | ------------------- | ------- | 
| humanTRA (500k)  | 15m26s | 1m00s | 15x |
| humanTRB (500k)  | 163m43s | 2m06s | **78x** |
| humanIGK (500k)  | 29m45s | 1m25s | 21x |
| humanIGL (500k)  | 24m47s | 1m04s | 23x |
| humanIGH (500k)* |  | 22m31s |  |
| humanIGH (10k)   | 50m34s | 27s | **112x** |

* For IGH, the original method was taking forever, so I had to switch to 10k sequences.

> Timings were recorded on a linux server with an Intel Xeon Gold 5220 CPU @ 2.20GHz
in a micromamba environment using python3.6 and numba 0.53.1.

## Note on floating-point precision

The new `numba` implementation produces Pgen values that are functionally equivalent to the original,
but since the number of float operations has changed, it is not possible to achieve an exact bit match. 
However, during testing, 100% of all sequences were <= 9 [ULPs]https://en.wikipedia.org/wiki/Unit_in_the_last_place)
of the old value, with 99.9% <= 4. Here's what that looks like in practice.

```bash
# 1 ULP off
1.5244400593599107e-16
1.5244400593599110e-16

# 4 ULP off (99.9% of data)
1.8518084083473990e-12
1.8518084083473973e-12

# 9 ULP off (highest observed)
3.4568224954178370e-12
3.4568224954178407e-12
```

Here is the histogram of ULP distances between the old and new Pgen scores for all models tested.

|ULP distance|IGH |IGK   |IGL   |TRA   |TRB   |
|------------|----|------|------|------|------|
|0           |2379|262490|321109|180113|135962|
|1           |3696|162691|130382|212673|210079|
|2           |2296|58032 |40049 |83953 |107167|
|3           |1042|13676 |7453  |19972 |36337 |
|4           |421 |2644  |900   |2934  |8650  |
|5           |118 |430   |99    |319   |1578  |
|6           |36  |34    |8     |36    |194   |
|7           |12  |4     |1     |1     |28    |
|8           |    |      |      |      |5     |
|9           |    |      |      |      |1     |

## Installation

The installation process is the same, just using this fork. You can directly `pip install` from github

```bash
# pip uninstall olga (if already installed)
pip install git+https://github.com/dweb0/olga
```

You will also need `numba` installed in your environment

```bash
# For pip
pip install numba

# For conda/mamba/micromamba 
conda install -c conda-forge numba
```

## Usage

To use from the command line, simply add the flag `--fast_pgen` to your regular `olga-compute_pgen`
command. For example,

```bash
olga-compute_pgen \
    --seq_in 0 \
    --v_in 2 \
    --j_in 3 \
    --lines_to_skip 1 \
    -d, \
    -i demo.csv \
    -o results.csv \
    --humanIGH \
    --fast_pgen  # <-- new flag
```

Library usage. Just 1 line to wrap the base model.

```python
from olga import GenerationProbabilityVJ
from olga.performance.fast_pgen import FastPgen

base_model = GenerationProbabilityVJ(...)

# Wrap the base model here
fast_model = FastPgen(base_model)

# Call compute_aa_CDR3_pgen as normal
p = fast_model.compute_aa_CDR3_pgen('CAWSVAPDRGGYTF', 'TRBV30', 'TRBJ1-2')

# Or use any of the other compute functions
p = fast_model.compute_nt_CDR3_pgen('TGTGCCAGTAGTATAACAACCCAGGGCTTGTACGAGCAGTACTTC')
```

> [!NOTE]
> The first invocation will have a bit of a delay since numba needs to compile (~ 45s), but
that only happens once.

## Limitations

- This submodule is only for computing Pgen, not generating sequences. It seems
there is already a tool for that though called [OLHA](https://github.com/statbiophys/olha).

## Tips

For even more throughput, you can use [SONIA](https://github.com/statbiophys/SONIA),
which provides multiprocessing. Just wrap `self.pgen_model` with `FastPgen`.

```python
# In sonia.evaluate_model.py.compute_all_pgens
pgen_model = FastPgen(self.pgen_model)
```
