import csv
import os
os.environ["POLARS_MAX_THREADS"] = "12"
import polars as pl
import argparse
from Bio import SeqIO
from itertools import groupby
from itertools import filterfalse
import numpy as np
from datetime import datetime as dt

# generate unique edge list with aggregate protein IDs

def map_kmer(kmer, dct):
    if kmer in dct.keys():
        return(dct[kmer])
    return(np.nan)

def unique(input_file, output_file, return_ids=False):
    
    # read input to dataframe & create new file name
    if not output_file:
        output_file = input_file.replace(".csv","")
        output_file = "{}_unique.csv".format(output_file)

    print(f"[{dt.now()}] Reading in de Bruijn graph from {input_file} ...")
    query = (
        pl.scan_csv(input_file,with_column_names=lambda cols: [col.lower() for col in cols])
        .select(pl.col("node1"),pl.col("node2").str.tail(1).alias("node2"))
        .select(pl.concat_str([pl.col("node1"),pl.col("node2")]).alias("node_edge"))
        .group_by("node_edge")
        .len("count")
        .select(
            pl.col("node_edge").str.head(-1).alias("node1"),
            pl.col("node_edge").str.tail(-1).alias("node2"),
            pl.col("count")
        )
        .sort("count",descending=True)
    )
    print(f"[{dt.now()}] Counting unique de Bruijn kmers ...")
    df = query.collect()
    
    # write the output to a new file
    print(f"[{dt.now()}] Writing results to {output_file} ...")
    df.write_csv(output_file)


def main():

    unique_kmer_file = unique(args.input_file, args.outfile, args.return_ids)

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Collapses an annotated deBruijn kmer file into uniques while aggregating IDs")
    parser.add_argument("--input_file", action="store", required=True,
                                        help="Filename for kmer edge list (.csv)")
    parser.add_argument("--sep", action="store", required=False, default=',',
                                        help="Column separator for input file, default=,")
    parser.add_argument("--outfile", action="store", required=False,
                                        help="(Optional) Filename for unique weighted kmer edge list (.csv)")

    parser.add_argument("--return_ids", action="store_true", default=False, help="(Optional) Return IDs associated with each unique edge; not recommended for huge data sets.")
    
    args = parser.parse_args()
    main()
