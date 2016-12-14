import glob
import os
import time
from collections import namedtuple
from typing import List, Set, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment

from atspsa_helper.create_lkh_full_matrix import create_full_atsp_via_scores_with_big_M, create_full_atsp_via_scores
from atspsa_helper.lkhrunner import run_lkh
from atspsacore import read_score_file, parser
from config import *

Node = namedtuple('Node', ['id', 'sequence', 'start_index', 'end_index', 'previous_node', 'next_node', 'length'])

circular = False


def align(seq1: str, seq2: str) -> str:
    """
    calls SeqAn and aligns two strings
    :param seq1:
    :param seq2:
    :return:
    """
    output = os.popen(
        '/home/andreas/GDrive/workspace/seqan-project/bin/seqconsensus --align {} {}'.format(seq1, seq2))
    return output.read()


def create_list_of_nodes(reads: List[Tuple[str, str]], tour: List[int], inclusions: List[int]) -> List[
    Node]:  # TODO remove?
    """
    creates a list of Nodes, a data structure containing useful data (see namedtuple above)
    List is sorted according to the tour
    :param reads:
    :param tour:
    :param inclusions:
    :return: nodes which are sorted according to the lkh output tour
    """
    nodes = []
    filtered_reads = [read for i, read in enumerate(reads) if i not in inclusions]
    # for i, read in enumerate(reads):
    for i, tour_elem in enumerate(tour):
        read = filtered_reads[tour_elem]
        start_index, end_index, length = read[0].split('/')[-1].split('_')
        if i < len(tour) - 1:
            next_node = tour[i + 1]
        else:
            next_node = -1
        if i > 0:
            previous_node = tour[i - 1]
        else:
            previous_node = -1
        nodes.append(
            Node(tour_elem, read[1], int(start_index), int(end_index), previous_node, next_node,
                 int(length)))
    nodes = sorted(nodes, key=lambda x: tour.index(x.id))
    return nodes


def compare_tours(tour, scores):
    """
    help method for controlling the tour result
    used for finding errors and unexpected behavior
    :param tour:
    :param scores:
    :return:
    """
    scores[scores < MINIMAL_OVERLAP_SCORE] = -BIG_M_WEIGHT
    should_be_tour_value = 0
    should_be_value_list = []
    lkh_tour_value = 0
    lkh_tour_value_list = []
    for i in range(len(tour) - 1):
        should_be_tour_value += scores[i][i + 1]
        should_be_value_list.append(scores[i][i + 1])
        lkh_tour_value += scores[tour[i]][tour[i + 1]]
        lkh_tour_value_list.append(scores[tour[i]][tour[i + 1]])
    if circular:
        should_be_tour_value += scores[-1][0]
        should_be_value_list.append(scores[-1][0])
        lkh_tour_value += scores[tour[-1]][tour[0]]
        lkh_tour_value_list.append(scores[tour[-1]][tour[0]])
    return 0


def remove_included_reads_in_scores(read_scores: np.ndarray, included_read_ids: List[int]) -> np.ndarray:
    """
    manipulates the scores matrix
    deletes all rows and columns of useless (i.e. contained in another sequence) reads
    :param read_scores: 2 dimensional np array containing the overlap scores of the reads
    :param included_read_ids: list of reads which are included by other reads (they will be removed
    :return: 2 dimensional np array without the rows and columns of the included_read_ids,
        note: pay attention to the indices; the id does not correspond to the index anymore
    """
    arr_with_reduced_cols = np.delete(read_scores, included_read_ids, 1)
    arr_with_reduced_rows = np.delete(arr_with_reduced_cols, included_read_ids, 0)
    return arr_with_reduced_rows


def check_for_gaps(nodes: List[Node], scores: np.ndarray, tour: List[int]):
    """
    TODO probably deprecated
    :param nodes:
    :param scores:
    :param tour:
    :return:
    """
    contig_counter = 0
    last_seq = nodes[0]
    contigs = [[last_seq.sequence]]
    inclusions = []
    for k in range(1, len(nodes)):
        output = align(last_seq.sequence, nodes[k].sequence)
        values = [int(x.split(': ')[1]) for x in output.split('\n')]
        if values[1] != values[0] + values[2] + values[3]:
            print('hi')
        if values[5] != values[4] + values[6] + values[7]:
            print('ho')

        if scores[last_seq.id][last_seq.next_node] == 0:
            print('Contig break')
            contig_counter += 1
            last_seq = nodes[k]
            contigs.append([last_seq.sequence])
            continue

        if values[2] == 0 and values[3] == 0 and values[6] == 0 and values[7] == 0:  # sequences are equal
            # i.e. no gaps at all
            if last_seq.sequence != nodes[k].sequence:
                print('What the heck?')
            pass
        elif values[2] == 0 and values[3] > 0 and values[6] >= 0 and values[7] == 0:  # normal overlap:
            # no leading gaps for first sequence
            # one or more trailing gaps for first sequence
            # zero or more leading gaps for second sequence
            # no trailing gaps for second sequence exists
            start_of_added_sequence = len(last_seq.sequence) - values[6]
            # xxxxxxyyyyyyyy------
            # ------yyyyyyyyzzzzzz
            #               ^- start of 'next' sequence
            # or
            # yyyyyyyy------
            # yyyyyyyyzzzzzz
            #         ^- start of 'next' sequence
            # leading gaps of sequence 2 = values[6]
            # the zzzzzzzzzz part is added to the list
            contigs[contig_counter].append(nodes[k].sequence[start_of_added_sequence:])
            last_seq = nodes[k]
        elif values[2] == 0 and values[3] == 0 and values[6] >= 0 and values[7] >= 0:  # inclusion of second sequence
            # xxxyyyyzzz
            # ---yyyy---
            # or
            # yyyyyzzzz
            # yyyyy----
            # or
            # xxxxxyyyy
            # -----yyyy
            # or
            # yyyyyyyyy     sequences are equal but actaully catched through the first condition
            # yyyyyyyyy
            # no gaps for first sequence
            # zero or more leading gaps for second sequence
            # one or more trailing gaps for second sequence
            # inclusions.append((last_seq.id, tour[k]))
            inclusions.append(tour[k])
        else:
            print('hmm')
    return contigs, inclusions


def make_new_lkh_run(filename: str, scores: np.ndarray, circular=False) -> None:  #TODO delete?
    """
    creates files for and runs LKH
    :param filename:
    :param scores:
    :param circular:
    :return:
    """
    create_full_atsp_via_scores_with_big_M(filename, scores, circular)
    run_lkh(filename + '.par')
    return


def look_for_included_reads(reads: List[str], scores: np.ndarray) -> List[int]:
    """
    filters all the reads which are included in another sequence
    :param reads: list of reads of the form List[sequence]
    :param scores: np.ndarray containing all the overlap scores
    :return:
    """
    inclusions = []
    for i, read in enumerate(reads):
        if len(read) == np.amax(scores, axis=0)[i]:
            inclusions.append(i)
    return inclusions


def detect_cycles(col_ind: np.ndarray) -> List[Set[int]]:
    """
    detects cycles of the assignment problem and stores them in a list of sets
    :param col_ind:
    :return: list of sets/cycles
    """
    elements_set = set(range(len(col_ind)))
    already_added_elements = set()
    current_element = len(col_ind) - 1  # elements_set.pop()  # elements_set[0]
    elements_set.remove(current_element)
    current_cycle = {current_element}
    # elements_set.remove(current_element)
    already_added_elements.add(current_element)
    cycles = []
    while elements_set:  # while set not empty
        new_element = col_ind[current_element]
        if new_element not in already_added_elements:  # is new_element==current_cycle[0] faster??
            current_cycle.add(new_element)
            already_added_elements.add(new_element)
            current_cycle.add(new_element)
            elements_set.remove(new_element)
            current_element = new_element
        else:
            new_element = elements_set.pop()  # [0]
            cycles.append(current_cycle)
            current_cycle = {new_element}
            already_added_elements.add(new_element)
            # elements_set.remove(new_element)
            current_element = new_element
    cycles.append(current_cycle)
    return cycles


def assemble_cycles(reads: List[str], cycles: List[Set[int]], col_ind: np.ndarray, scores: np.ndarray) -> List[str]:
    """
    method for assembling the cycles (i.e. List[int]) to contigs
    :param reads: list of sequences, no ids
    :param cycles: each set corresponds to a cycle
    :param col_ind: indicating the next(assigned) nodes in the assignment problem
    :param scores:
    :return: List[str]
    """
    contig_list = []
    for cycle in cycles:
        # find starting point in cycle
        min_score = BIG_M_WEIGHT
        if len(col_ind) - 1 in cycle:
            current_element = len(col_ind) - 1
        else:
            for elem in cycle:
                if min_score > scores[elem][col_ind[elem]]:
                    min_score = scores[elem][col_ind[elem]]
                    current_element = col_ind[elem]
        cycle.remove(current_element)
        if current_element == len(col_ind) - 1:
            current_contig = ['']
        else:
            current_contig = [reads[current_element]]
        while cycle:
            next_element = col_ind[current_element]
            if next_element == len(col_ind) - 1:
                cycle.remove(next_element)
                current_element = next_element
                continue
            if current_element == len(col_ind) - 1:
                begin_index = 0
            else:
                begin_index = scores[current_element, next_element]
            current_contig.append(reads[next_element][begin_index:])
            cycle.remove(next_element)
            current_element = next_element
        contig_list.append(''.join(current_contig))
    return contig_list


def munkres(reads: List[str], scores: np.ndarray) -> Tuple[List[str], float, int]:
    """
    applying the (negative) scores matrix to the munkres/linear_sum_assignment algorithm
    we detect the cycles and assemble them afterwards
    :param reads: only a list of the sequences! No ids!
    :param scores: actual weights have positive values (the values are gonna be inverted) and negative BIG_M_WEIGHT
    :return:
    """
    # scores[scores < MINIMAL_OVERLAP_SCORE] = 0
    # scores_copy = scores.copy()
    read_lengths = [len(x) for x in reads]
    summed_read_length = sum(read_lengths)
    manipulated_scores = manipulate_scores2ap(scores, read_lengths)
    start_time_ap = time.time()
    row_ind, col_ind = linear_sum_assignment(manipulated_scores)
    # row_ind, col_ind = linear_sum_assignment(manipulated_scores)
    ap_value = manipulated_scores[row_ind, col_ind].sum()
    end_time_ap = time.time()

    cycles = detect_cycles(col_ind)
    # ap_value += summed_read_length * (len(cycles) - 1)


    contig_list = assemble_cycles(reads, cycles, col_ind, scores)

    print(contig_list)
    print(end_time_ap - start_time_ap)
    return contig_list, end_time_ap - start_time_ap, ap_value


def manipulate_scores2atsp(scores: np.ndarray, read_lengths: List[int]) -> np.ndarray:
    n = len(read_lengths)
    scores = np.array(read_lengths).reshape(n, 1) - scores
    summed_read_length = sum(read_lengths)
    new_column = np.zeros((len(read_lengths), 1), dtype=np.int32)
    new_row = np.zeros((1, len(read_lengths) + 1), dtype=np.int32)
    scores = np.vstack((np.hstack((scores, new_column)), new_row))  # add column and row
    for i in range(n):
        # scores[i] = read_lengths[i] - scores[i]
        # for j in range(n):
        #     if i != j:
        #         scores[i][j] = read_lengths[i] - scores[i][j]
        scores[i][n] = read_lengths[i]
        scores[n][i] = summed_read_length
    return scores


def manipulate_scores2ap(scores: np.ndarray, read_lengths: List[int]) -> np.ndarray:
    n = len(read_lengths)
    scores = np.array(read_lengths).reshape(n, 1) - scores
    summed_read_length = sum(read_lengths)
    new_column = np.zeros((len(read_lengths), 1), dtype=np.int32)
    new_row = np.zeros((1, len(read_lengths) + 1), dtype=np.int32)
    scores = np.vstack((np.hstack((scores, new_column)), new_row))  # add column and row
    for i in range(n):
        # for j in range(n):
        #     if i != j:
        #         scores[i][j] = read_lengths[i] - scores[i][j]
        #     else:
        #         scores[i][j] = summed_read_length
        scores[i][i] = read_lengths[i]  # summed_read_length
        scores[i][n] = read_lengths[i]
        scores[n][i] = 0
    scores[n][n] = summed_read_length
    return scores


def assemble_tour(scores: np.ndarray, reads: List[str], tour: List[int]) -> List[str]:
    """
    assembles the reads to contigs given the data of scores, reads (only strings) and the tour
    the scores determine at which point the next read should start to be assembled
    :param scores:
    :param reads:
    :param tour:
    :return:
    """
    contig_list = []
    current_elem = tour[0]
    current_contig = [reads[current_elem]]
    for next_elem in tour[1:]:
        current_score = scores[current_elem][next_elem]
        if current_score > 0:
            current_contig.append(reads[next_elem][current_score:])
        else:
            contig_list.append(''.join(current_contig))
            current_contig = [reads[next_elem]]
        current_elem = next_elem
    contig_list.append(''.join(current_contig))
    return contig_list


def lkh_assemble(reads: List[Tuple[str, str]], scores: np.ndarray) -> Tuple[
    List[str], float, int]:
    """
    starting lkh run (i.e. calling make_new_lkh_run)
    parses the output tour and creates a List[Node]
    which is used for the assemble_tour method to assemble the contigs
    :param reads:
    :param scores:
    :return:
    """

    filtered_reads_lengths = [len(x[1]) for i, x in enumerate(reads)]
    read_strings = [x[1] for i, x in enumerate(reads)]
    manipulated_scores = manipulate_scores2atsp(scores, filtered_reads_lengths)
    create_full_atsp_via_scores('/home/andreas/GDrive/workspace/atspsa/footest', manipulated_scores)
    start_time_tsp = time.time()
    run_lkh('/home/andreas/GDrive/workspace/atspsa/footest.par')
    end_time_tsp = time.time()
    tour = parser.parse_tour('/home/andreas/GDrive/workspace/atspsa/footest.tour')
    with open('/home/andreas/GDrive/workspace/atspsa/footest.tour', 'r') as f:
        tsp_value = int(
            [x.split('COMMENT : Length = ')[1] for x in f.readlines() if x.startswith('COMMENT : Length = ')][0]) - sum(
            filtered_reads_lengths)
    tsp_contigs = assemble_tour(scores, read_strings, tour)

    print(tsp_contigs)
    # compare_tours(tour, scores)

    print(end_time_tsp - start_time_tsp)
    return tsp_contigs, end_time_tsp - start_time_tsp, tsp_value


def assemble(fasta_file: str, score_file: str, tsp_bool=True, ap_bool=True):
    """
    method for assembly via atsp/LKH or Hungarian Alogrithm
    :param fasta_file:
    :param score_file:
    :param tsp:
    :param ap:
    :return:
    """
    # read scores, reads
    nr_of_reads, scores = read_score_file.read_score_file(score_file)
    reads_with_id = parser.parse_fasta_with_id(fasta_file)
    print('Fasta file read.')
    # detect inclusions
    # inclusions = look_for_included_reads([read[1] for i, read in enumerate(reads)], scores)
    # reduce the scores matrix
    # reduced_scores = remove_included_reads_in_scores(scores, inclusions)

    # TSP
    if tsp_bool:
        lkh_contigs, lkh_time, lkh_value = lkh_assemble(reads_with_id, scores)
    else:
        lkh_contigs, lkh_time, lkh_value = [], -1.0, -1

    if ap_bool:
        # Assignment problem
        reads_str_only = [read[1] for i, read in enumerate(reads_with_id)]
        ap_contigs, ap_time, ap_value = munkres(reads_str_only, scores)
    else:
        ap_contigs, ap_time, ap_value = [], -1.0, -1
    with open(score_file.split('.score')[0] + '.assembly', 'w') as f:
        f.write('LKH: Circular: {}\nLKH_Contigs:\n'.format(circular))
        f.write('\n'.join(lkh_contigs))
        f.write('\nLKH_Objective_Value: {}\n'.format(lkh_value))
        f.write('LKH_Time: {}\n'.format(lkh_time))
        f.write('Assignment Problem:\nAP_Contigs:\n')
        f.write('\n'.join(ap_contigs))
        f.write('\nAP_Objective_Value: {}'.format(ap_value))
        f.write('\nAP_Time: {}'.format(ap_time))


if __name__ == "__main__":

    for file in sorted(glob.glob('/home/andreas/GDrive/workspace/sparsedata/ref2_c40_l100/calign*.score')):
        print(file)
        current_dir = os.path.dirname(os.path.realpath(file))
        fasta_file = current_dir + '/reads.fasta'
        t0 = time.time()
        assemble(fasta_file, file, False, True)
        t1 = time.time()
        print('Complete Time: ', t1 - t0)
