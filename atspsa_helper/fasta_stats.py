import glob
import math
from typing import List, Tuple

from atspsacore.parser import parse_fasta_with_id, parse_fasta
from config import *


def filter_included_reads(reads: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
    """
    filters all the reads which are included in another sequence
    :param reads: list of reads with ids of the form List[id, sequence]
    :return:
    """
    inclusions = []
    last_begin_index, last_end_index, _ = [int(x) for x in reads[0][0].split('/')[-1].split('_')]
    for i in range(1, len(reads)):
        next_begin_index, next_end_index, _ = [int(x) for x in reads[i][0].split('/')[-1].split('_')]
        if last_begin_index <= next_begin_index and last_end_index >= next_end_index:
            inclusions.append(i)
        else:
            last_begin_index, last_end_index = next_begin_index, next_end_index
    return [x for i, x in enumerate(reads) if i not in inclusions]


def compute_objective_function_value(filtered_reads: List[Tuple[str, str]]) -> int:
    objective_value = 0
    tar_list = []
    last_begin_index, last_end_index, _ = [int(x) for x in filtered_reads[0][0].split('/')[-1].split('_')]
    counter = 0
    for i in range(1, len(filtered_reads)):
        read_id = filtered_reads[i][0]
        next_begin_index, next_end_index, _ = [int(x) for x in read_id.split('/')[-1].split('_')]
        if last_end_index - next_begin_index >= MINIMAL_OVERLAP_SCORE:
            objective_value += last_end_index - next_begin_index
            tar_list.append(last_end_index - next_begin_index)
        else:
            objective_value += -BIG_M_WEIGHT
            counter += 1
            tar_list.append(-BIG_M_WEIGHT)
        last_begin_index, last_end_index = next_begin_index, next_end_index
    return objective_value


def compute_objective_function_value_of_fasta(fasta_file: str) -> int:
    reads = parse_fasta_with_id(fasta_file)
    filtered_reads = filter_included_reads(reads)
    return compute_objective_function_value(filtered_reads)


def write_fasta_stat(fasta_file, ref_file):
    with open(ref_file, 'r') as f:
        ref_seq = f.read()
        ref_length = len(ref_seq)

    reads = parse_fasta(fasta_file)
    read_lengths = [len(read) for read in reads]

    with open(fasta_file.split('reads.fasta')[0] + "fasta.stat", "w") as handle:
        real_coverage = sum([x for x in read_lengths]) / ref_length
        nr_reads = len(read_lengths)
        average = sum([x for x in read_lengths]) / len(read_lengths)
        median = sorted([x for x in read_lengths])[len(read_lengths) // 2]
        shortest = min([x for x in read_lengths])
        longest = max([x for x in read_lengths])
        objective_value = compute_objective_function_value_of_fasta(fasta_file)
        expected_gaps = len(reads) * math.exp(-real_coverage)
        unsequenced_part = math.exp(-real_coverage)
        handle.write('Real Coverage: {}\n'
                     'Number of reads: {}\n'
                     'Average read length: {}\n'
                     'Median read length: {}\n'
                     'Shortest read length" : {}\n'
                     'Longest read length: {}\n'
                     'Objective function value: {}\n'
                     'Expected gaps: {}\n'
                     'Unsequenced part: {}'.format(real_coverage, nr_reads, average, median, shortest, longest,
                                                   objective_value, expected_gaps,
                                                   unsequenced_part))
    return (real_coverage, nr_reads, average, median, shortest, longest, objective_value, expected_gaps)


def make_latex_stat_table(stats, ref_number: int):
    output = '\\begin{tabular}{|c|c|c|c|c|c|c|c|}\\hline\n& \\multicolumn{4}{c|}{Read Length} & & & \\\\\nCoverage & Average & Median & Shortest' \
             '& Longest & \\#Reads & Objective Value & Exp. Gaps \\\\\n\\hline\n'
    for coverage in coverages:
        for average_length in average_length_list:
            val = stats[(coverage, average_length)]
            output = output + '{0:.1f} & {1:.1f} & {2} & {3} & {4} & {5} & {6} & {7:.2f}\\\\\n'.format(val[0], val[2],
                                                                                                       val[3],
                                                                                                       val[4], val[5],
                                                                                                       val[1],
                                                                                                       val[6], val[7])
        output = output + '\\hline\n'
    output = output + '\\end{tabular}'
    with open(DIR + 'ref{}_stats.tex'.format(ref_number), 'w') as f:
        f.write(output)


if __name__ == '__main__':
    stats = {}
    for file in glob.glob(DIR + 'ref1_c*/reads.fasta'):
        coverage = int(file.split('_c')[1].split('_l')[0])
        average_length = int(file.split('_l')[1].split('/')[0])
        stats[(coverage, average_length)] = write_fasta_stat(file, ref1)
    make_latex_stat_table(stats, 1)
    stats = {}
    for file in glob.glob(DIR + 'ref2_c*/reads.fasta'):
        coverage = int(file.split('_c')[1].split('_l')[0])
        average_length = int(file.split('_l')[1].split('/')[0])
        stats[(coverage, average_length)] = write_fasta_stat(file, ref2)
    make_latex_stat_table(stats, 2)
