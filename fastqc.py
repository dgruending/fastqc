from argparse import ArgumentParser
from os import path

import numpy as np
import pyfastx


def aggregate(fastq, quality_scores, start_index=0, end_index=0):
    """

    :return:
    """
    for ind in range(start_index, end_index+1):
        seq = fastq[ind].seq
        qualities = fastq[ind].quali
        quality_scores[ind, 0:len(qualities)] = qualities

if __name__ == '__main__':
    parser = ArgumentParser()
    # TODO add help text
    parser.add_argument("-f", "--file", dest="file_path", help="", required=True)
    parser.add_argument("-o", "--output", dest="output", help="")

    args = parser.parse_args()

    if path.exists(args.file_path):
        fq_file = pyfastx.Fastq(args.file_path)
        qual_scores = np.zeros((len(fq_file), fq_file.maxlen))
        aggregate(fq_file, qual_scores, start_index=0, end_index=len(fq_file)-1)
    else:
        raise FileNotFoundError(f"FastQ file \"{args.file_path}\" does not exist.")
