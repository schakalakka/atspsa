#################
# reads all the .(c)act files in every directory and puts them in a dictionary
# creates two line graph plots (currently not used because of tikz)
# writes a CSV file alignmenttime_stats.csv containing the raw data in one file
# additionally writes a .tex file alignmenttime_stats.tex containing a table
#


import os
from collections import OrderedDict

import matplotlib.pyplot as plt
import matplotlib.pyplot as ticker

from config import *


def plot(filename):
    fig = plt.figure()
    fig.suptitle('Alignment Computation Benchmark', fontsize=14)
    ax1 = fig.add_axes([0.2, 0.1, 0.7, 0.77])

    ax1.set_ylabel("Time in seconds")
    ax1.set_xlabel("Average read length in hundred base pairs")
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_yscale('log')
    ax1.plot([i + 1 for i, x in enumerate(seq_align_dict.values()) if x > 0],
             [x for x in seq_align_dict.values() if x > 0], 'ro-', label='SeqAn')
    ax1.plot([i + 1 for i, x in enumerate(seq_m20_20_align_dict.values()) if x > 0],
             list([x for x in seq_m20_20_align_dict.values() if x > 0]), 'go-', label='SeqAn20')
    ax1.plot([i + 1 for i, x in enumerate(c_align_dict.values()) if x > 0],
             list([x for x in c_align_dict.values() if x > 0]), 'bo-', label='Calign')
    ax1.grid(True)
    ax1.set_xticks(range(27))
    # ax1.set_xticklabels([x[2] // 100 1for x in list(seq_align_dict.keys())])
    ax1.xaxis.set_major_formatter(ticker.NullFormatter())
    ax1.xaxis.set_major_locator(ticker.FixedLocator([i + 1 for i in range(27)]))
    ax1.xaxis.set_major_formatter(ticker.FixedFormatter([x[2] // 100 for x in list(seq_align_dict.keys())]))
    # ax1.xaxis.set_minor_locator(ticker.FixedLocator([i+1 for i in range(27)]))
    # ax1.xaxis.set_minor_formatter(ticker.FixedFormatter([x[2] // 100 for x in list(seq_align_dict.keys())]))
    # ax1.xaxis.set_ticks_position('bottom')


    legend = ax1.legend(['SeqAn', 'SeqAn_-20_20', 'Calign'], loc='lower right')

    ax2 = fig.add_axes([0.2, 0.87, 0.7, 0])
    ax2.yaxis.set_visible(False)
    ax2.set_xticks(range(10))
    ax2.xaxis.set_major_formatter(ticker.NullFormatter())
    ax2.xaxis.set_major_locator(ticker.FixedLocator(range(10)))
    ax2.xaxis.set_minor_locator(ticker.FixedLocator([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]))
    ax2.xaxis.set_minor_formatter(ticker.FixedFormatter(['C5', 'C20', 'C40', 'C5', 'C20', 'C40', 'C5', 'C20', 'C40']))
    ax2.xaxis.set_ticks_position('bottom')

    ax3 = fig.add_axes([0.2, 0.871, 0.7, 0])
    ax3.yaxis.set_visible(False)
    ax3.set_xticks(range(4))
    ax3.xaxis.set_major_formatter(ticker.NullFormatter())
    ax3.xaxis.set_major_locator(ticker.FixedLocator([0, 1, 2, 3]))
    ax3.xaxis.set_minor_locator(ticker.FixedLocator([0.5, 1.5, 2.5]))
    ax3.xaxis.set_minor_formatter(ticker.FixedFormatter([ref1_name, ref2_name, ref3_name]))
    ax3.xaxis.set_ticks_position('top')

    plt.savefig(filename + "_plot.pdf")
    plt.close()


def plotbanded(filename):
    fig = plt.figure()
    fig.suptitle('Alignment Computation Benchmark', fontsize=14)
    ax1 = fig.add_axes([0.2, 0.1, 0.7, 0.77])

    ax1.set_ylabel("Time in seconds")
    ax1.set_xlabel("Average read length in hundred base pairs")
    ax1.xaxis.set_ticks_position('bottom')
    ax1.set_yscale('log')
    ax1.plot([i for i, x in enumerate(seq_m5_5_align_dict.values()) if x > 0],
             [x for x in seq_m5_5_align_dict.values() if x > 0], 'ro-', label='SeqAn5')
    ax1.plot([i for i, x in enumerate(seq_m20_20_align_dict.values()) if x > 0],
             list([x for x in seq_m20_20_align_dict.values() if x > 0]), 'go-', label='SeqAn20')

    # ax1.grid(True)
    ax1.set_xticks(range(len(seq_align_dict)))
    ax1.set_xticklabels([x[2] // 100 for x in list(seq_align_dict.keys())])

    legend = ax1.legend(['SeqAn_-5_5', 'SeqAn_-20_20'], loc='lower right')

    ax2 = fig.add_axes([0.2, 0.87, 0.7, 0])
    ax2.yaxis.set_visible(False)
    ax2.set_xticks(range(10))
    ax2.xaxis.set_major_formatter(ticker.NullFormatter())
    ax2.xaxis.set_major_locator(ticker.FixedLocator(range(10)))
    ax2.xaxis.set_minor_locator(ticker.FixedLocator([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5]))
    ax2.xaxis.set_minor_formatter(ticker.FixedFormatter(['C5', 'C20', 'C40', 'C5', 'C20', 'C40', 'C5', 'C20', 'C40']))
    ax2.xaxis.set_ticks_position('bottom')

    ax3 = fig.add_axes([0.2, 0.871, 0.7, 0])
    ax3.yaxis.set_visible(False)
    ax3.set_xticks(range(4))
    ax3.xaxis.set_major_formatter(ticker.NullFormatter())
    ax3.xaxis.set_major_locator(ticker.FixedLocator([0, 1, 2, 3]))
    ax3.xaxis.set_minor_locator(ticker.FixedLocator([0.5, 1.5, 2.5]))
    ax3.xaxis.set_minor_formatter(ticker.FixedFormatter([ref1_name, ref2_name, ref3_name]))
    ax3.xaxis.set_ticks_position('top')

    plt.savefig(filename + "_plot.pdf")
    plt.close()


def write_csv(filename):
    csv = 'Nr\tRef\tCoverage\tLength\tSeqAn\tSeqAn(0,max)\tSeqAn(-5,5)\tSeqAn(-20,20)\tSeqAn(0,5)\tSeqAn(0,20)\tCalign\n'
    counter = 1
    for ref in [1, 2, 3]:
        for coverage in coverages:
            for length in average_length_list:
                csv += '{0}\t{1}\t{2}\t{3}\t{4:.0f}\t{5:.0f}\t{6:.0f}\t{7:.0f}\t{8:.0f}\t{9:.0f}\t{10:.0f}\n'.format(
                    counter, ref,
                    coverage,
                    length,
                    seq_align_dict[
                        (ref,
                         coverage,
                         length)],
                    seq_0_max_align_dict[
                        (
                            ref,
                            coverage,
                            length)],
                    seq_m5_5_align_dict[
                        (ref,
                         coverage,
                         length)],
                    seq_m20_20_align_dict[
                        (ref,
                         coverage,
                         length)],
                    seq_0_5_align_dict[
                        (ref,
                         coverage,
                         length)],
                    seq_0_20_align_dict[
                        (ref,
                         coverage,
                         length)],
                    c_align_dict[
                        (ref,
                         coverage,
                         length)])
                counter += 1
    with open(filename + '_stats.csv', 'w') as f:
        f.write(csv)
    print(output)


def get_alignment_time(filename):
    if os.path.exists(DIR + 'ref{0}_c{1}_l{2}/'.format(ref, coverage, length) + filename):
        with open(DIR + 'ref{0}_c{1}_l{2}/'.format(ref, coverage, length) + filename,
                  'r') as f:
            align_time = float(f.read())
    else:
        align_time = -1
    return align_time


output = '\\begin{tabular}{|c|c|c|c|c|c|c|c|}\\hline\n\multicolumn{4}{|c|}{Data set} & \multicolumn{4}{c|}{Aligner}\\\\\nRef & Coverage & Length & Reads & SeqAn & SeqAn5 & SeqAn 20 & Calign \\\\\n\\hline\n'
seq_align_dict = OrderedDict()
seq_m20_20_align_dict = OrderedDict()
seq_m5_5_align_dict = OrderedDict()
seq_0_20_align_dict = OrderedDict()
seq_0_5_align_dict = OrderedDict()
seq_0_max_align_dict = OrderedDict()
c_align_dict = OrderedDict()
for ref in [1, 2, 3]:
    for coverage in coverages:
        for length in average_length_list:
            nr_of_reads = -1
            seq_align_time = get_alignment_time('seqan.time')
            seq_m20_20_align_time = get_alignment_time('seqan_-20_20.time')
            seq_m5_5_align_time = get_alignment_time('seqan_-5_5.time')
            seq_0_20_align_time = get_alignment_time('seqan_0_20.time')
            seq_0_5_align_time = get_alignment_time('seqan_0_5.time')
            seq_0_max_align_time = get_alignment_time('seqan_0_max.time')
            c_align_time = get_alignment_time('calign.time')
            with open(DIR + 'ref{0}_c{1}_l{2}/fasta_stats.txt'.format(ref, coverage, length),
                      'r') as f:
                f.readline()
                nr_of_reads = f.readline().split(': ')[1]
            output = output + '{0} & {1} & {2} & {3} & {4:.0f} & {5:.0f} & {6:.0f} & {7:.0f}\\\\\n'.format(ref,
                                                                                                           coverage,
                                                                                                           length,
                                                                                                           nr_of_reads,
                                                                                                           seq_align_time,
                                                                                                           seq_m5_5_align_time,
                                                                                                           seq_m20_20_align_time,
                                                                                                           c_align_time)
            seq_align_dict[(ref, coverage, length)] = seq_align_time
            seq_m20_20_align_dict[(ref, coverage, length)] = seq_m20_20_align_time
            seq_m5_5_align_dict[(ref, coverage, length)] = seq_m5_5_align_time
            seq_0_20_align_dict[(ref, coverage, length)] = seq_0_20_align_time
            seq_0_5_align_dict[(ref, coverage, length)] = seq_0_5_align_time
            seq_0_max_align_dict[(ref, coverage, length)] = seq_0_max_align_time
            c_align_dict[(ref, coverage, length)] = c_align_time
    output = output + '\\hline\n'
output = output + '\\end{tabular}'
with open(DIR + 'alignmenttime' + '_stats.tex', 'w') as f:
    f.write(output)
# print(output)

write_csv(DIR + 'alignmenttime')
# plot(DIR + 'alignmenttime')
# plotbanded(DIR + 'bandedalignmenttime')
