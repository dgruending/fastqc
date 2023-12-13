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
    for ind in range(start_index, end_index + 1):
        # seq = fastq[ind].seq
        qualities = fastq[ind].quali
        quality_scores[ind, 0:len(qualities)] = qualities


def qual_per_base(quali_arr, plot=False):
    base_count = quali_arr.shape[1]
    if plot:
        pass
    else:
        df = pd.DataFrame(index=range(1, base_count+1))
        df.index.name = "Base"
        df["Mean"] = np.mean(quali_arr, axis=0)
        df["Median"] = np.median(quali_arr, axis=0)
        df["Lower Quartile"] = np.percentile(quali_arr, 25, axis=0)
        df["Upper Quartile"] = np.percentile(quali_arr, 75, axis=0)
        df["10th Percentile"] = np.percentile(quali_arr, 10, axis=0)
        df["90th Percentile"] = np.percentile(quali_arr, 90, axis=0)
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
    aggregate(fq_file, qual_scores, start_index=0, end_index=num_seqs - 1)

    qual_per_base(qual_scores)
