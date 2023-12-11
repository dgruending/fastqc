from argparse import ArgumentParser
from os import path

import numpy as np
import pyfastx


def read_fastq_file(filepath):
    if path.exists(filepath):
        return pyfastx.Fastq(args.file_path)
    else:
        raise FileNotFoundError(f"FastQ file \"{args.file_path}\" does not exist.")


def aggregate(fastq, quality_scores, start_index=0, end_index=0):
    """

    :return:
    """
    for ind in range(start_index, end_index + 1):
        seq = fastq[ind].seq
        qualities = fastq[ind].quali
        quality_scores[ind, 0:len(qualities)] = qualities


def qual_per_base(quali_arr):
    pass


if __name__ == '__main__':
    parser = ArgumentParser()
    # TODO add help text
    parser.add_argument("-f", "--file", dest="file_path", help="", required=True)
    parser.add_argument("-o", "--output", dest="output", help="")

    args = parser.parse_args()

    fq_file = read_fastq_file(args.file_path)

    num_seqs = len(fq_file)

    qual_scores = np.full((num_seqs, fq_file.maxlen), np.nan)
    aggregate(fq_file, qual_scores, start_index=0, end_index=num_seqs - 1)

    qual_per_base(qual_scores)
