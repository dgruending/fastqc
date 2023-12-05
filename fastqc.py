from argparse import ArgumentParser
from os import path
import pyfastx

if __name__ == '__main__':
    parser = ArgumentParser()
    # TODO add help text
    parser.add_argument("-f", "--file", dest="file_path", help="", required=True)
    parser.add_argument("-o", "--output", dest="output", help="")

    args = parser.parse_args()

    if path.exists(args.file_path):
        fq_file = pyfastx.Fastq(args.file_path)
