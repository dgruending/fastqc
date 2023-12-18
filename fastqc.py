import numbers
import re
from argparse import ArgumentParser
from enum import Enum
from os import path
import numpy as np
import pyfastx
import pandas as pd


class Base(Enum):
    G = 0
    A = 1
    T = 2
    C = 3
    N = 4

    def __str__(self):
        if self == Base.A:
            return "A"
        elif self == Base.C:
            return "C"
        elif self == Base.G:
            return "G"
        elif self == Base.T:
            return "T"
        else:
            return "N"

    def __eq__(self, other):
        if self.__class__ is other.__class__:
            return self.value == other.value
        elif isinstance(other, numbers.Number):
            return self.value == other
        elif isinstance(other, str):
            return self.__str__() == other
        else:
            raise ValueError(f"{type(other)} can not be compared: {other}")


def read_fastq_file(filepath):
    if path.exists(filepath):
        return pyfastx.Fastq(args.file_path)
    else:
        raise FileNotFoundError(f"FastQ file \"{args.file_path}\" does not exist.")


def aggregate(fastq, quality_scores, motifs=None, start_index=0, end_index=0):
    """

    :return:
    """
    lengths = np.zeros(fastq.maxlen)
    base_content = np.zeros((5, fastq.maxlen))

    for ind in range(start_index, end_index + 1):
        read = fastq[ind]
        seq = read.seq
        qualities = read.quali

        # get quality scores for per base sequence quality statistic
        quality_scores[ind, 0:len(qualities)] = qualities

        # get length data for length distribution
        lengths[len(read) - 1] += 1

        # get base content information
        # TODO optimize?
        for base_ind in range(len(seq)):
            base_content[Base[seq[base_ind]].value, base_ind] += 1

        # motif matching
        for mot_ind, motif in enumerate(motifs):
            for match in re.finditer(motif, seq):
                pass

    return lengths, base_content


def qual_per_base(quali_arr, output_file, plot=False):
    base_count = quali_arr.shape[1]
    if plot:
        pass
    else:
        df = pd.DataFrame(index=range(1, base_count + 1))
        df.index.name = "Base"
        # use nan methods to ignore nan values. There shouldn't be empty slices to worry about.
        df["Mean"] = np.nanmean(quali_arr, axis=0)
        df["Median"] = np.nanmedian(quali_arr, axis=0)
        df["Lower Quartile"] = np.nanpercentile(quali_arr, 25, axis=0)
        df["Upper Quartile"] = np.nanpercentile(quali_arr, 75, axis=0)
        df["10th Percentile"] = np.nanpercentile(quali_arr, 10, axis=0)
        df["90th Percentile"] = np.nanpercentile(quali_arr, 90, axis=0)
        df["Minimum"] = np.nanmin(quali_arr, axis=0)
        df["Maximum"] = np.nanmax(quali_arr, axis=0)
        df.to_csv(output_file, sep="\t")


def length_distribution(lengths, output_file, plot=False):
    if plot:
        pass
    else:
        non_zero_indices = np.nonzero(lengths)[0]
        df = pd.DataFrame(index=non_zero_indices + 1, data=lengths[non_zero_indices], columns=["Count"], dtype=int)
        df.index.name = "Length"
        df.to_csv(output_file, sep="\t")


def per_base_seq_content(base_counts, output_file, plot=False):
    base_count = base_counts.shape[1]
    sums = np.sum(base_counts, axis=0)
    df = pd.DataFrame(index=range(1, base_count + 1))
    df.index.name = "Base"
    for base in Base:
        df[str(base)] = base_counts[base.value, :] / sums * 100
    df.to_csv(output_file, sep="\t")


if __name__ == '__main__':
    parser = ArgumentParser()
    # TODO add help text
    parser.add_argument("-f", "--file", dest="file_path", help="", required=True)
    parser.add_argument("-o", "--output", dest="output_prefix", help="", default="./")
    parser.add_argument("-m", "--motifs", dest="motifs", nargs='+', help="", default=[])
    parser.add_argument("-p", "--plot", dest="plotting", action="store_true", help="")
    parser.add_argument("-a", "--polya", dest="poly_a", action="store_true", help="")
    parser.add_argument("-n", "--njobs", dest="n_jobs", help="", default=1)

    args = parser.parse_args()

    fq_file = read_fastq_file(args.file_path)

    num_seqs = len(fq_file)

    # prepare storage for quality scores
    qual_scores = np.full((num_seqs, fq_file.maxlen), np.nan)
    # compile motifs as regex pattern
    pattern_list = [re.compile(pattern) for pattern in args.motifs]
    if args.poly_a:
        pattern_list.append(re.compile(r"A{4,}"))
    length_data, bases_data = aggregate(fq_file, qual_scores, motifs=pattern_list, start_index=0, end_index=num_seqs - 1)

    # TODO use os.path for path data
    qual_per_base(qual_scores, args.output_prefix + "seq_qual_per_base.csv")
    length_distribution(length_data, args.output_prefix + "seq_length.csv")
    per_base_seq_content(bases_data, args.output_prefix + "per_base_seq_content.csv")
