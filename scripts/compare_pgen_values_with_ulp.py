#!/usr/bin/env python3
import sys
import math
import struct
from pathlib import Path
from functools import reduce

import click
import polars as pl #v0.18.0
import numpy as np


def ulp_distance_1(a: float, b: float) -> int:
    a, b = sorted((a, b))
    cur = a
    ulps = 0
    while cur != b:
        nex = np.nextafter(cur, b)
        ulps += 1
        cur = nex
    
    return ulps


def ulp_distance_2(a: float, b: float) -> int:
    if math.isnan(a) or math.isnan(b):
        return math.inf

    ai = struct.unpack('>q', struct.pack('>d', a))[0]
    bi = struct.unpack('>q', struct.pack('>d', b))[0]

    if ai < 0:
        ai = (0x8000000000000000 - ai) & 0xFFFFFFFFFFFFFFFF
    if bi < 0:
        bi = (0x8000000000000000 - bi) & 0xFFFFFFFFFFFFFFFF

    return abs(ai - bi)


@click.command()
@click.argument('files', type=click.Path(exists=True), nargs=-1, required=True)
def main(files):
    """ Compare old pgen results vs new by calculating the ULP distance.

    Input is a 2-col csv containing the old Pgen and new Pgen scores.
    (can be built from from `xsv cat columns old.csv new.csv`)
    """

    hists = []
    for file in files:
        label = Path(file).stem

        df = pl.read_csv(
            file,
            has_header=False,
            new_columns=['old', 'new']
        )

        # We'll use 2 implementations of ULP distance in case one is faulty
        df = (
            df.with_columns([
                pl.struct(['old', 'new'])
                    .apply(lambda x: ulp_distance_1(x['old'], x['new']))
                    .alias('ulp_distance'),
                pl.struct(['old', 'new'])
                    .apply(lambda x: ulp_distance_2(x['old'], x['new']))
                    .alias('ulp_distance_2'),
            ])
        )

        if not df.filter(pl.col('ulp_distance') != pl.col('ulp_distance_2')).is_empty():
            raise RuntimeError('1 of the ULP distance functions is off.')

        df_hist = (
            df
            .groupby('ulp_distance')
            .count()
            .sort('ulp_distance')
            .with_columns(
                pl.col('count').cumsum().alias('cumulative')
            )
            .with_columns(
                (pl.col('cumulative') / pl.col('cumulative').max()) * 100.0
            )
            .select(['ulp_distance', 'count'])
            .rename({'count': label})
        )
        hists.append(df_hist)

    #print (df_agg)
    
    df_hist = reduce(lambda l, r: l.join(r, on='ulp_distance', how='outer'), hists)
    df_hist.write_csv(sys.stdout.buffer)


if __name__ == "__main__":
    main()
