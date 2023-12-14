from argparse import ArgumentParser
from os import path
import numpy as np
import pyfastx
import pandas as pd


def read_fastq_file(filepath):
    if path.exists(filepath):
        return pyfastx.Fastq(args.file_path)
    else:
        raise FileNotFoundError(f"FastQ file \"{args.file_path}\" does not exist.")


def aggregate(fastq, quality_scores, start_index=0, end_index=0):
    """

    :return:
    """
    lengths = np.zeros(fastq.maxlen)

    for ind in range(start_index, end_index + 1):
        read = fastq[ind]
        # seq = read.seq
        qualities = read.quali

        # get quality scores for per base sequence quality statistic
        quality_scores[ind, 0:len(qualities)] = qualities

        # get length data for length distribution
        lengths[len(read) - 1] += 1

    return lengths


def qual_per_base(quali_arr, plot=False):
    base_count = quali_arr.shape[1]
    if plot:
        pass
    else:
        df = pd.DataFrame(index=range(1, base_count + 1))
        df.index.name = "Base"
        df["Mean"] = np.mean(quali_arr, axis=0)
        df["Median"] = np.median(quali_arr, axis=0)
        df["Lower Quartile"] = np.percentile(quali_arr, 25, axis=0)
        df["Upper Quartile"] = np.percentile(quali_arr, 75, axis=0)
        df["10th Percentile"] = np.percentile(quali_arr, 10, axis=0)
        df["90th Percentile"] = np.percentile(quali_arr, 90, axis=0)
        return df


def length_distribution(lengths, plot=False):
    if plot:
        pass
    else:
        non_zero_indices = np.nonzero(lengths)[0]
        df = pd.DataFrame(index=non_zero_indices + 1, data=lengths[non_zero_indices], columns=["Count"], dtype=int)
        df.index.name = "Length"
        return df


if __name__ == '__main__':
    parser = ArgumentParser()
    # TODO add help text
    parser.add_argument("-f", "--file", dest="file_path", help="", required=True)
    parser.add_argument("-o", "--output", dest="output", help="")

    args = parser.parse_args()

    fq_file = read_fastq_file(args.file_path)

    num_seqs = len(fq_file)

    qual_scores = np.full((num_seqs, fq_file.maxlen), np.nan)
    length_data = aggregate(fq_file, qual_scores, start_index=0, end_index=num_seqs - 1)

    qual_per_base(qual_scores)
    length_distribution(length_data)
