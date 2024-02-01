import multiprocessing as mp
import numbers
import os.path
import sys
from argparse import ArgumentParser
from enum import Enum
from os import path
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import pyfastx
import suffix_tree.tree
from suffix_tree import Tree
import time
import logging

from viztracer import VizTracer

TIMING = True


class Base(Enum):
    G = 0
    A = 1
    T = 2
    C = 3
    N = 4

    def __str__(self):
        """
        Single letter string representation of the Base.

        :return: string
        """
        if self == Base.A:
            # needed because Base.N matches everything
            if self == Base.C:
                return "N"
            return "A"
        elif self == Base.C:
            return "C"
        elif self == Base.G:
            return "G"
        elif self == Base.T:
            return "T"

    # Can compare(equal) to Base, number and string
    def __eq__(self, other):
        """
        Check for equality with Base, number or string.

        :param other:  Base, number or string.
        :return: bool
        """
        # if self = Base.N
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
            return str(self) == other
        else:
            raise ValueError(f"{type(other)} can not be compared: {other}")


def find_index(tree, pattern):
    """
    Finds the first index, where the pattern is fully represented or reaches the text end.

    :param tree: suffix_tree to be searched.
    :param pattern: string-like object.
    :return: number of pattern start position or a number >= len(text)
    """
    node = tree.root
    return _find_index(node, pattern, 0)


def _find_index(node, pattern, pattern_ind):
    """
    Recursive pattern search in the subtree.

    :param node: currently searched node
    :param pattern: string-like object.
    :param pattern_ind: match position in pattern
    :return: number of pattern start position or a number >= len(text)
    """
    if isinstance(node, suffix_tree.node.Leaf):
        match_ind = node.start
        # match rest of pattern against rest of text
        for (char_motif, char_seq) in zip(pattern[pattern_ind:], node.S[match_ind + pattern_ind:]):
            # Text ended or mismatch
            if not isinstance(char_seq, suffix_tree.util.UniqueEndChar) and Base[char_motif] != Base[char_seq]:
                return node.end
        return node.start

    else:
        # check for pattern ending
        if pattern_ind == len(pattern):
            return node.start
        keys = [key for key in node.children.keys() if isinstance(key, suffix_tree.util.UniqueEndChar)
                or Base[key] == Base[pattern[pattern_ind]]]
        indices = [_find_index(node.children[key], pattern, pattern_ind + 1) for key in keys]
        if indices:
            return min(indices)
        else:
            return sys.maxsize


def get_barcode_ind(barcode):
    """
    Extract barcode indices from string representation.

    :param barcode: string-like object; representing a barcode.
    :return: list of indices
    """
    if barcode:
        return [ind for ind, char in enumerate(barcode, start=1) if char == 'X']
    else:
        return []


def read_fastq_file(file_path):
    """
    Reads a Fastq-file.

    :param file_path: path-like object to file
    :return: Fastq-file object
    """
    if path.exists(file_path):
        return pyfastx.Fastq(args.file_path)
    else:
        raise FileNotFoundError(f"FastQ file \"{file_path}\" does not exist.")


def box_plot(df, arr, name, save_loc, indexes=None):
    """
    Makes a plotly boxplot for a given array of quality values.

    :param df: dataframe with Median, Lower_Quartile and Upper_Quartile
    :param arr: array of quality values.
    :param name: string-like object for plot title
    :param save_loc: path-like object for plot save location
    :param indexes: iterable of columns to use; if None all are used; default=None
    """
    # TODO add y-axis label
    layout = go.Layout(showlegend=False, title=name)
    if indexes is None:
        indexes = range(arr.shape[1])
    inter_quantile_range = df["Upper_Quartile"].values - df["Lower_Quartile"].values
    traces = []
    for base_ind in indexes:
        box_trace = go.Box(
            x=[base_ind + 1],
            q1=[df["Lower_Quartile"].values[base_ind]],
            median=[df["Median"].values[base_ind]],
            q3=[df["Upper_Quartile"].values[base_ind]],
            upperfence=[df["Upper_Quartile"].values[base_ind] + 1.5 * inter_quantile_range[base_ind]],
            lowerfence=[df["Lower_Quartile"].values[base_ind] - 1.5 * inter_quantile_range[base_ind]],
            name=f"Base {base_ind + 1}",
            marker=dict(color="blue")
        )
        traces.append(box_trace)
    traces.append(go.Scatter(x=df.index.values + 1, y=df["Mean"].values, line=dict(color="red"), mode="lines",
                             name="Mean"))
    # TODO add outliers
    fig = go.Figure(data=traces, layout=layout)
    fig.write_html(save_loc)


def quality_per_base(quality_arr, zoom_indexes, output_path, plot=False):
    """
    Makes base position quality statistic.
    Calculates mean, median, minimum, maximum, 10-, 25-, 75- and 90-quartile.

    :param quality_arr: array(sequence x base_position) of quality values
    :param zoom_indexes: list or array of positions to save/plot additionally
    :param output_path: path-like object to save location
    :param plot: boolean value, if plotting will be done; defaults to False
    """
    base_count = quality_arr.shape[1]
    df = pd.DataFrame(index=range(1, base_count + 1))
    df.index.name = "Base"
    # quality array could have values of 0 representing no quality entry,
    # because for example the sequence has ended early
    # TODO handle missing values separately
    df["Mean"] = np.mean(quality_arr, axis=0)
    df["Median"] = np.median(quality_arr, axis=0)
    df["Lower_Quartile"] = np.percentile(quality_arr, 25, axis=0)
    df["Upper_Quartile"] = np.percentile(quality_arr, 75, axis=0)
    df["10th_Percentile"] = np.percentile(quality_arr, 10, axis=0)
    df["90th_Percentile"] = np.percentile(quality_arr, 90, axis=0)
    df["Minimum"] = np.min(quality_arr, axis=0)
    df["Maximum"] = np.max(quality_arr, axis=0)
    df.to_csv(os.path.join(output_path, "seq_qual_per_base.csv"), sep="\t")
    if plot:
        box_plot(df, quality_arr, "Sequence Quality per Base",
                 os.path.join(output_path, "seq_qual_per_base.html"))
    if zoom_indexes:
        zoom_df = df.loc[zoom_indexes]
        zoom_df.to_csv(os.path.join(output_path, "zoom.csv"), sep="\t")
        if plot:
            box_plot(df, quality_arr, "Barcode quality", os.path.join(output_path, "barcode_quality.html"),
                     indexes=zoom_indexes)


def length_plot(data, output_path):
    """
    Makes a plotly line plot for the length distribution.

    :param data: dataframe with "Count" column
    :param output_path: path-like object for save location
    """
    fig = px.line(data, y="Count", x=[data.index], title="Length Counts")
    fig.update_layout(xaxis_title="Sequence length", showlegend=False)
    fig.write_html(os.path.join(output_path, "length_distribution.html"))


def length_distribution(lengths, output_path, plot=False):
    """
    Saves only non-zero entries of lengths as csv-file.

    :param lengths: 1d-array of length counts
    :param output_path: path-like object for save location
    :param plot: boolean value, if plotting will be done; defaults to False
    """
    non_zero_indices = np.nonzero(lengths)[0]
    df = pd.DataFrame(index=range(1, len(lengths) + 1), data=lengths, columns=["Count"], dtype=int)
    df.index.name = "Length"
    # save only non-zero entries
    df.iloc[non_zero_indices].to_csv(os.path.join(output_path, "seq_length.csv"), sep="\t")
    if plot:
        length_plot(df, output_path)


def seq_content_plot(data, save_file):
    """
    Makes plotly line plot with base content for each base.

    :param data: dataframe with "Base" index and Bases as columns
    :param save_file: path-like object for save location
    """
    data = data.reset_index()
    fig = px.line(data, x="Base", y=data.columns[1:], title="Sequence content",
                  color_discrete_sequence=["grey", "green", "red", "blue", "black"])
    fig.update_layout(yaxis_title="Base content %")
    fig.write_html(save_file)


def per_base_seq_content(base_counts, output_path, plot=False):
    """
    Makes per position base content statistic.

    :param base_counts: array with counts for each base;
    row for each base(are assumed to be in Base-enum order), column is position
    :param output_path: path-like object for save location
    :param plot: boolean value, if plotting will be done; defaults to False
    """
    base_count = base_counts.shape[1]
    sums = np.sum(base_counts, axis=0)
    df = pd.DataFrame(index=range(1, base_count + 1))
    df.index.name = "Base"
    for base in Base:
        # TODO catch error if base count is to low
        df[str(base)] = base_counts[base.value, :] / sums * 100
    df.to_csv(os.path.join(output_path, "per_base_seq_content.csv"), sep="\t")
    if plot:
        seq_content_plot(df, os.path.join(output_path, "base_seq_content.html"))


def plot_motif_counts(data, save_loc):
    """
    Makes plotly line plot of motif-start count at a given position.

    :param data: dataframe with "Base" index and motifs as column names
    :param save_loc: path-like object for save location
    """
    data = data.reset_index()
    fig = px.line(data, x="Base", y=data.columns[1:], title="Motif start positions")
    fig.update_layout(yaxis_title="Count")
    fig.write_html(save_loc)


def motif_statistic(motif_data, motifs, output_path, plot=False):
    """
    Saves motif counts as csv-file and optionally plots them.

    :param motif_data: array of motif counts
    :param motifs: iterable of motifs
    :param output_path: path-like object for save location
    :param plot: boolean value, if plotting will be done; defaults to False
    """
    non_zero_columns_ind = np.unique(np.nonzero(motif_data)[1])
    df = pd.DataFrame(index=range(1, motif_data.shape[1] + 1))
    df.index.name = "Base"
    for ind, motif in enumerate(motifs):
        df[motif] = motif_data[ind]
    # only save non-zero columns
    df.iloc[non_zero_columns_ind].to_csv(os.path.join(output_path, "motif_counts.csv"), sep="\t")
    if plot:
        plot_motif_counts(df, os.path.join(output_path, "motif_counts.html"))


def main(file_path, numb_seqs, max_len, patterns_list, barcode_indices, plotting, n_jobs, output_path):
    block_size = numb_seqs / n_jobs
    # due to int conversion and potential float block_size even the last sequence is considered
    arg_list = [(file_path, max_len, patterns_list, (int(i * block_size), int((i + 1) * block_size)), i) for i in
                range(n_jobs)]
    with mp.Pool(processes=n_jobs) as pool:
        result = pool.starmap(aggregate, arg_list)
        length_data = 0
        base_content_data = 0
        motif_data = 0
        # quality data is doubled in memory, could be a place for memory optimization
        qual_data = np.zeros((numb_seqs, max_len), dtype=np.uint8)

        # an index is needed to get the start and end position for saving quality data
        ind = 0
        for qual_scores, lengths, bases, motif_count in result:
            block_start, block_end = arg_list[ind][3]
            qual_data[block_start:block_end, :] = qual_scores
            length_data += lengths
            base_content_data += bases
            motif_data += motif_count
            ind += 1

        # use multiple process, if available, to make statistics, they are independent of each other
        processes = [pool.apply_async(quality_per_base, (qual_data, barcode_indices, output_path, plotting)),
                     pool.apply_async(length_distribution, (length_data, output_path, plotting)),
                     pool.apply_async(per_base_seq_content, (base_content_data, output_path, plotting)),
                     pool.apply_async(motif_statistic, (motif_data, patterns_list, output_path, plotting))]
        # waiting for started process to avoid confusion and terminate main-program if everything is done
        for process in processes:
            process.wait()


def aggregate(file_path, max_len, motifs, block, process_id):
    """
    Aggregate quality, length, base count and pattern count data.

    :param process_id: number used to distinguish processes for debugging and profiling
    :param file_path: path-like object to fastq file
    :param max_len: maximal sequence length
    :param motifs: iterable of motifs/pattern
    :param block: tuple with inclusive start index and exclusive end index of sequences to process
    :return: arrays for quality, length, base count and pattern count data
    """
    with VizTracer(output_file=f"{process_id}.json", tracer_entries=block[1]*10) as tracer:
        fq = read_fastq_file(file_path)
        start_time = time.time()
        start, end = block
        logging.debug(f"Process {process_id}: Starting aggregation: {end - start} sequences to do.")
        # init arrays
        quality_scores = np.zeros((end - start, fq_file.maxlen), dtype=np.uint8)
        logging.debug(f"Process {process_id}: {sys.getsizeof(quality_scores) / 1024 ** 3} GB allocated for quality scores.")
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
            if motifs:
                tree = Tree({0: seq})
                for motif_ind, motif in enumerate(motifs):
                    motif_start = find_index(tree, motif)
                    if motif_start < len(seq):
                        motifs_occ[motif_ind, motif_start] += 1

            # debugging log
            if (start == 0 and qual_ind % 1000 == 0) or qual_ind % 1000000 == 0:
                logging.debug(f"Process {process_id}: {qual_ind + 1} sequences done elapsed time "
                              f"{(time.time() - start_time)} s")

    return quality_scores, lengths, base_content, motifs_occ


if __name__ == '__main__':
    logging.basicConfig(filename='debug.log', filemode='w', encoding='utf-8', level=logging.DEBUG)
    parser = ArgumentParser()
    # could do settings as file, for easier calls
    parser.add_argument("-f", "--file", dest="file_path", help="un-/zipped fastq file", required=True)
    parser.add_argument("-o", "--output", dest="output", help="directory for output data", default="./")
    parser.add_argument("-m", "--motifs", dest="motifs", nargs='+', help="search motifs", default=[])
    parser.add_argument("-p", "--plot", dest="plotting", action="store_true", help="set to make plots")
    parser.add_argument("-a", "--polya", dest="poly_a", type=int,
                        help="minimal length of Poly-A tail to consider")
    parser.add_argument("-n", "--njobs", dest="n_jobs",
                        help="Number of processes to create for data aggregation, evaluation and plotting.",
                        default=1, type=int)
    parser.add_argument("-b", "--bar", dest="barcode", help="Sequence of Barcode format to extract indices")

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
    # check if output directory exists, else make it
    Path(args.output).mkdir(parents=True, exist_ok=True)

    if TIMING:
        time_start = time.time()
    main(args.file_path, num_seqs, fq_file.maxlen, pattern_list, barcode_ind, args.plotting, args.n_jobs,
         args.output)
    if TIMING:
        time_end = time.time()
        # noinspection PyUnboundLocalVariable
        logging.info(f"The time of execution of above program is : {(time_end - time_start)} s")
