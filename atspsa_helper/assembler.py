import os
import time
from collections import namedtuple
from typing import List, Set, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment

from atspsa_helper.create_lkh_full_matrix import create_full_atsp_via_scores
from atspsa_helper.lkhrunner import run_lkh
from atspsacore import read_score_file, parser
from config import *

Node = namedtuple('Node', ['id', 'sequence', 'start_index', 'end_index', 'previous_node', 'next_node', 'length'])

circular = True


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


def create_list_of_nodes(reads: List[Tuple[str, str]], tour: List[int], inclusions: List[int]) -> List[Node]:
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
        start_index, end_index = read[0].split('/')[-1].split('_')
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
                 int(end_index) - int(start_index)))
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
    # scores[scores < 30] = -99999
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


def make_new_lkh_run(filename: str, scores: np.ndarray, circular=False) -> None:
    """
    creates files for and runs LKH
    :param filename:
    :param scores:
    :param circular:
    :return:
    """
    create_full_atsp_via_scores(filename, scores, circular)
    run_lkh(filename + '.par')
    return


def look_for_included_reads(reads: List[Tuple[str]], scores: np.ndarray) -> List[int]:
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
    current_element = elements_set.pop()  # elements_set[0]
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
        current_element = cycle.pop()
        current_contig = [reads[current_element]]
        while cycle:
            next_element = col_ind[current_element]
            current_contig.append(reads[next_element][scores[current_element, next_element]:])
            cycle.remove(next_element)
            current_element = next_element
        contig_list.append(''.join(current_contig))
    return contig_list


def munkres(reads: List[str], scores: np.ndarray) -> Tuple[List[str], float]:
    """
    applying the (negative) scores matrix to the munkres/linear_sum_assignment algorithm
    we detect the cycles and assemble them afterwards
    :param reads: only a list of the sequences! No ids!
    :param scores: actual weights have positive values (the values are gonna be inverted) and negative BIG_M_WEIGHT
    :return:
    """
    start_time_ap = time.time()
    # scores[scores < MINIMAL_OVERLAP_SCORE] = 0
    row_ind, col_ind = linear_sum_assignment(np.negative(scores))
    # print('Cost: ', scores[row_ind, col_ind].sum())
    cycles = detect_cycles(col_ind)
    contig_list = assemble_cycles(reads, cycles, col_ind, scores)
    end_time_ap = time.time()
    print(contig_list)
    print(end_time_ap - start_time_ap)
    return contig_list, end_time_ap - start_time_ap


def assemble_tour(nodes: List[Node], scores: np.ndarray) -> List[str]:
    """
    assembles the reads to contigs given the data of the nodes (and Node structure/namedtuple) and the scores
    the scores determine at which point the next read should start to be assembled
    TODO it is not clear so far how it should be handled with a banded alignment score
    :param nodes: sorted according to the tour (i.e. next_node)
    :param scores:
    :return: list of contigs
    """
    contig_list = []
    current_contig = [nodes[0].sequence]
    for i in range(0, len(nodes) - 1):
        if scores[nodes[i].id][nodes[i].next_node] > MINIMAL_OVERLAP_SCORE:
            current_contig.append(nodes[i + 1].sequence[scores[nodes[i].id][nodes[i].next_node]:])
        else:
            contig_list.append(current_contig)
            current_contig = [nodes[i + 1].sequence]
    contig_list.append(current_contig)
    return [''.join(contig) for contig in contig_list]


def lkh_assemble(reads: List[Tuple[str, str]], inclusions: List[int], scores: np.ndarray) -> Tuple[List[str], float]:
    """
    starting lkh run (i.e. calling make_new_lkh_run)
    parses the output tour and creates a List[Node]
    which is used for the assemble_tour method to assemble the contigs
    :param reads:
    :param inclusions:
    :param scores:
    :return:
    """
    start_time_tsp = time.time()
    make_new_lkh_run('/home/andreas/GDrive/workspace/atspsa/footest', scores, circular)
    tour = parser.parse_tour('/home/andreas/GDrive/workspace/atspsa/footest.tour', circular)
    nodes = create_list_of_nodes(reads, tour, inclusions)
    tsp_contigs = assemble_tour(nodes, scores)
    # contig_list, inclusions = check_for_gaps(nodes, scores, tour)
    # tsp_contigs = [''.join(contig) for contig in contig_list]
    print(tsp_contigs)
    # compare_tours(tour, scores)
    end_time_tsp = time.time()
    print(end_time_tsp - start_time_tsp)
    return tsp_contigs, end_time_tsp - start_time_tsp


def assemble(fasta_file: str, score_file: str):
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
    reads = parser.parse_fasta_with_id(fasta_file)
    reads_seqs = parser.parse_fasta(fasta_file)
    # detect inclusions
    inclusions = look_for_included_reads(reads_seqs, scores)
    # reduce the scores matrix
    reduced_scores = remove_included_reads_in_scores(scores, inclusions)

    # TSP
    lkh_contigs, lkh_time = lkh_assemble(reads, inclusions, reduced_scores)

    # Assignment problem
    reduced_reads = [read[1] for i, read in enumerate(reads) if i not in inclusions]
    ap_contigs, ap_time = munkres(reduced_reads, reduced_scores)
    with open(score_file.split('.score')[0] + '.assembly', 'w') as f:
        f.write('LKH: Circular: {}\nLKH_Contigs:\n'.format(circular))
        f.write('\n'.join(lkh_contigs))
        f.write('\nLKH_Time: {}\n'.format(lkh_time))
        f.write('Assignment Problem:\nAP_Contigs:\n')
        f.write('\n'.join(ap_contigs))
        f.write('\nAP_Time: {}'.format(ap_time))


if __name__ == "__main__":
    file = '/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c5_l700/calign.score'
    fasta_file = '/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c5_l700/shuffled_c5_l700.fasta'
    # for file in glob.glob('/home/andreas/GDrive/workspace/sparsedata/ref*/calign.score'):
    #     fasta_file = glob.glob(file.split('calign.score')[0] + '*.fasta')[0]
    assemble(fasta_file, file)
