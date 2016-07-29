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
    ax1.plot([i for i, x in enumerate(seq_align_dict.values()) if x > 0],
             [x for x in seq_align_dict.values() if x > 0], 'ro-', label='SeqAn')
    ax1.plot([i for i, x in enumerate(seq20_align_dict.values()) if x > 0],
             list([x for x in seq20_align_dict.values() if x > 0]), 'go-', label='SeqAn20')
    ax1.plot([i for i, x in enumerate(c_align_dict.values()) if x > 0],
             list([x for x in c_align_dict.values() if x > 0]), 'bo-', label='Calign')
    # ax1.grid(True)
    ax1.set_xticks(range(len(seq_align_dict)))
    ax1.set_xticklabels([x[2] // 100 for x in list(seq_align_dict.keys())])

    legend = ax1.legend(['SeqAn', 'SeqAn20', 'Calign'], loc='lower right')

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
    ax1.plot([i for i, x in enumerate(seq5_align_dict.values()) if x > 0],
             [x for x in seq5_align_dict.values() if x > 0], 'ro-', label='SeqAn5')
    ax1.plot([i for i, x in enumerate(seq20_align_dict.values()) if x > 0],
             list([x for x in seq20_align_dict.values() if x > 0]), 'go-', label='SeqAn20')

    # ax1.grid(True)
    ax1.set_xticks(range(len(seq_align_dict)))
    ax1.set_xticklabels([x[2] // 100 for x in list(seq_align_dict.keys())])

    legend = ax1.legend(['SeqAn5', 'SeqAn20'], loc='lower right')

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


output = '\\begin{tabular}{|c|c|c|c|c|c|c|}\\hline\nRef & Coverage & Length & Reads & SeqAn & SeqAn 20 & Calign \\\\\n\\hline\n'
seq_align_dict = OrderedDict()
seq20_align_dict = OrderedDict()
seq5_align_dict = OrderedDict()
c_align_dict = OrderedDict()
for ref in [1, 2, 3]:
    for coverage in coverages:
        for length in average_length_list:
            nr_of_reads = -1
            if os.path.exists(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}.act'.format(ref, coverage, length)):
                with open(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}.act'.format(ref, coverage, length),
                          'r') as f:
                    seq_align_time = float(f.read())
            else:
                seq_align_time = -1
            if os.path.exists(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}_banded20.act'.format(ref, coverage,
                                                                                                      length)):
                with open(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}_banded20.act'.format(ref, coverage,
                                                                                                  length),
                          'r') as f:
                    seq20_align_time = float(f.read())
            else:
                seq20_align_time = -1
            if os.path.exists(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}_banded5.act'.format(ref, coverage,
                                                                                                     length)):
                with open(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}_banded5.act'.format(ref, coverage, length),
                          'r') as f:
                    seq5_align_time = float(f.read())
            else:
                seq5_align_time = -1
            if os.path.exists(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}.cact'.format(ref, coverage, length)):
                with open(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}.cact'.format(ref, coverage, length),
                          'r') as f:
                    c_align_time = float(f.read())
            else:
                c_align_time = -1
            with open(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}.stat'.format(ref, coverage, length),
                      'r') as f:
                f.readline()
                nr_of_reads = f.readline().split(': ')[1]
            output = output + '{0} & {1} & {2} & {3} & {4:.0f} & {5:.0f} & {6:.0f}\\\\\n'.format(ref, coverage, length,
                                                                                                 nr_of_reads,
                                                                                                 seq_align_time,
                                                                                                 seq20_align_time,
                                                                                                 c_align_time)
            seq_align_dict[(ref, coverage, length)] = seq_align_time
            seq20_align_dict[(ref, coverage, length)] = seq20_align_time
            seq5_align_dict[(ref, coverage, length)] = seq5_align_time
            c_align_dict[(ref, coverage, length)] = c_align_time
    output = output + '\\hline\n'
output = output + '\\end{tabular}'
with open(DIR + 'alignmenttime' + '_stats.tex', 'w') as f:
    f.write(output)
print(output)

plot(DIR + 'alignmenttime')
plotbanded(DIR + 'bandedalignmenttime')
