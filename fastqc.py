import numbers
import os.path
import sys
from argparse import ArgumentParser
from enum import Enum
from os import path
import multiprocessing as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyfastx
import suffix_tree.tree
from suffix_tree import Tree
import plotly.express as px


class Base(Enum):
    G = 0
    A = 1
    T = 2
    C = 3
    N = 4

    def __str__(self):
        if self == Base.A:
            if self == Base.C:
                return "N"
            return "A"
        elif self == Base.C:
            return "C"
        elif self == Base.G:
            return "G"
        elif self == Base.T:
            return "T"

    def __eq__(self, other):
        if self.value == 4:
            return True
        if self.__class__ is other.__class__:
            if other.value == 4:
                return True
            else:
                return self.value == other.value
        elif isinstance(other, numbers.Number):
            return self.value == other
        elif isinstance(other, str):
            return self.__str__() == other
        else:
            raise ValueError(f"{type(other)} can not be compared: {other}")


def find_index(tree, motif):
    node = tree.root
    return _find_index(node, motif, 0)


def _find_index(node, motif, ind):
    if isinstance(node, suffix_tree.node.Leaf):
        match_ind = node.start
        for (char_motif, char_seq) in zip(motif[ind:], node.S[match_ind + ind:]):
            if not isinstance(char_seq, suffix_tree.util.UniqueEndChar) and Base[char_motif] != Base[char_seq]:
                return node.end
        return node.start

    else:
        if ind == len(motif):
            return node.start
        keys = [key for key in node.children.keys() if isinstance(key, suffix_tree.util.UniqueEndChar)
                or Base[key] == Base[motif[ind]]]
        indices = [_find_index(node.children[key], motif, ind + 1) for key in keys]
        if indices:
            return min(indices)
        else:
            return sys.maxsize


def get_barcode_ind(barcode):
    if barcode:
        return [ind for ind, char in enumerate(barcode, start=1) if char == 'X']
    else:
        return []


def read_fastq_file(filepath):
    if path.exists(filepath):
        return pyfastx.Fastq(args.file_path)
    else:
        raise FileNotFoundError(f"FastQ file \"{args.file_path}\" does not exist.")


def bar_plot(data, name):
    fig, ax = plt.subplots()
    ax.boxplot(data)
    fig.savefig(name)


def plot_qual_per_base(data, output_path, zoom=None):
    px.box(data)


def qual_per_base(quali_arr, zoom_indexes, output_path, plot=False):
    base_count = quali_arr.shape[1]
    df = pd.DataFrame(index=range(1, base_count + 1))
    df.index.name = "Base"
    df["Mean"] = np.mean(quali_arr, axis=0)
    df["Median"] = np.median(quali_arr, axis=0)
    df["Lower Quartile"] = np.percentile(quali_arr, 25, axis=0)
    df["Upper Quartile"] = np.percentile(quali_arr, 75, axis=0)
    df["10th Percentile"] = np.percentile(quali_arr, 10, axis=0)
    df["90th Percentile"] = np.percentile(quali_arr, 90, axis=0)
    df["Minimum"] = np.min(quali_arr, axis=0)
    df["Maximum"] = np.max(quali_arr, axis=0)
    df.to_csv(os.path.join(output_path, "seq_qual_per_base.csv"), sep="\t")
    if plot:
        bar_plot(quali_arr, "qual_per_base.png")
        plot_qual_per_base(pd.DataFrame(quali_arr), output_path, zoom=zoom_indexes)
    if zoom_indexes:
        zoom_df = df.loc[zoom_indexes]
        zoom_df.to_csv(os.path.join(output_path, "zoom.csv"), sep="\t")
        if plot:
            bar_plot(quali_arr[zoom_indexes], "barcode_error.png")


def length_plot(data, non_zero, output_path):
    # start, stop = non_zero
    # if start != 0:
    #     start -= 1
    # if stop != len(data) - 1:
    #     stop += 1
    # x = np.array(range(start, stop + 1))
    # fig, ax = plt.subplots()
    # ax.plot(x + 1, data[x], 'r')
    # fig.savefig(name)
    fig = px.line(data, y="Count", x=[data.index], title="Length Counts")
    fig.update_layout(xaxis_title="Sequence length", showlegend=False)
    fig.write_html(os.path.join(output_path, "length_distribution.html"))


def length_distribution(lengths, output_path, plot=False):
    non_zero_indices = np.nonzero(lengths)[0]
    df = pd.DataFrame(index=range(1, len(lengths) + 1), data=lengths, columns=["Count"], dtype=int)
    df.index.name = "Length"
    df.iloc[non_zero_indices].to_csv(os.path.join(output_path, "seq_length.csv"), sep="\t")
    if plot:
        length_plot(df, (non_zero_indices[0], non_zero_indices[-1]), output_path)


def seq_content_plot(data, name):
    data = data.reset_index()
    fig = px.line(data, x="Base", y=data.columns[1:], title="Sequence content")
    fig.update_layout(yaxis_title="Base content %")
    fig.write_html(name)


def per_base_seq_content(base_counts, output_path, plot=False):
    base_count = base_counts.shape[1]
    sums = np.sum(base_counts, axis=0)
    df = pd.DataFrame(index=range(1, base_count + 1))
    df.index.name = "Base"
    for base in Base:
        df[str(base)] = base_counts[base.value, :] / sums * 100
    df.to_csv(os.path.join(output_path, "per_base_seq_content.csv"), sep="\t")
    if plot:
        seq_content_plot(df, os.path.join(output_path, "base_seq_content.html"))


def plot_motif_counts(data, name):
    data = data.reset_index()
    fig = px.line(data, x="Base", y=data.columns[1:], title="Motif start positions")
    fig.update_layout(yaxis_title="Count")
    fig.write_html(name)


def motif_statistic(motif_data, motifs, output_path, plot=False):
    non_zero_columns_ind = np.unique(np.nonzero(motif_data)[1])
    df = pd.DataFrame(index=range(1, motif_data.shape[1] + 1))
    df.index.name = "Base"
    for ind, motif in enumerate(motifs):
        df[motif] = motif_data[ind]
    df.iloc[non_zero_columns_ind].to_csv(os.path.join(output_path, "motif_counts.csv"), sep="\t")
    if plot:
        plot_motif_counts(df, os.path.join(output_path, "motif_counts.html"))


def main(file_path, num_seqs, max_len, pattern_list, barcode_ind, plotting, n_jobs, output_path):
    block_size = num_seqs / n_jobs
    arg_list = [(file_path, max_len, pattern_list, (int(i * block_size), int((i + 1) * block_size))) for i in
                range(n_jobs)]
    with mp.Pool(processes=n_jobs) as pool:
        result = pool.starmap(aggregate_async, arg_list)
        length_data = 0
        base_content_data = 0
        motif_data = 0
        qual_data = np.zeros((num_seqs, max_len), dtype=np.uint8)

        ind = 0
        for qual_scores, lengths, bases, motif_count in result:
            block_start, block_end = arg_list[ind][3]
            qual_data[block_start:block_end, :] = qual_scores
            length_data += lengths
            base_content_data += bases
            motif_data += motif_count
            ind += 1

        processes = []
        processes.append(pool.apply_async(qual_per_base, (qual_data, barcode_ind, output_path)))
        processes.append(pool.apply_async(length_distribution, (length_data, output_path, plotting)))
        processes.append(pool.apply_async(per_base_seq_content, (base_content_data, output_path, plotting)))
        processes.append(pool.apply_async(motif_statistic, (motif_data, pattern_list, output_path, plotting)))
        for process in processes:
            process.wait()


def aggregate_async(file_path, max_len, motifs, block):
    fq = read_fastq_file(file_path)
    start, end = block
    quality_scores = np.zeros((end - start, fq_file.maxlen), dtype=np.uint8)
    lengths = np.zeros(max_len)
    base_content = np.zeros((5, max_len))
    motifs_occ = np.zeros((len(motifs), max_len))
    for qual_ind, ind in enumerate(range(start, end)):
        read = fq[ind]
        seq = read.seq
        qualities = read.quali

        # get quality scores for per base sequence quality statistic
        quality_scores[qual_ind, 0:len(qualities)] = qualities

        # get length data for length distribution
        lengths[len(read) - 1] += 1

        # get base content information
        # TODO optimize?
        for base_ind in range(len(seq)):
            base_content[Base[seq[base_ind]].value, base_ind] += 1

        # motif matching
        tree = Tree({0: seq})
        for motif_ind, motif in enumerate(motifs):
            motif_start = find_index(tree, motif)
            if motif_start < len(seq):
                motifs_occ[motif_ind, motif_start] += 1

    return quality_scores, lengths, base_content, motifs_occ


if __name__ == '__main__':
    parser = ArgumentParser()
    # TODO add help text
    parser.add_argument("-f", "--file", dest="file_path", help="", required=True)
    parser.add_argument("-o", "--output", dest="output", help="", default="./")
    parser.add_argument("-m", "--motifs", dest="motifs", nargs='+', help="", default=[])
    parser.add_argument("-p", "--plot", dest="plotting", action="store_true", help="")
    parser.add_argument("-a", "--polya", dest="poly_a", type=int, help="")
    parser.add_argument("-n", "--njobs", dest="n_jobs", help="", default=1, type=int)
    parser.add_argument("-b", "--bar", dest="barcode", help="")

    args = parser.parse_args()

    fq_file = read_fastq_file(args.file_path)

    num_seqs = len(fq_file)
    pattern_list = args.motifs
    if args.poly_a:
        pattern_list.append("A" * args.poly_a)
    # get relevant barcode data (indices)
    barcode_ind = get_barcode_ind(args.barcode)
    if args.n_jobs <= 0:
        args.n_jobs = 1
    main(args.file_path, num_seqs, fq_file.maxlen, pattern_list, barcode_ind, args.plotting, args.n_jobs, args.output)
