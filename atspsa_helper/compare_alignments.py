import csv
import glob

import numpy as np

from config import *

QUALITY_THRESHOLDS = [0, 0.7, 0.9, 0.99]
ALIGNERS = ['SeqAn', 'SeqAn(0,max)', 'SeqAn(-5,5)', 'SeqAn(-20,20)', 'SeqAn(0,5)', 'SeqAn(0,20)', 'Calign']


def read_score_file(file, np_arr, k):
    with open(file, 'r') as f:
        f.readline()
        for line in f.readlines():
            i, j, val = [int(elem) for elem in line.split('\t')]
            np_arr[k][i][j] = val


def read_score_files(*args):
    with open(args[0], 'r') as f:
        nr_reads = int(f.readline())
    np_arr = np.zeros((len(args), nr_reads, nr_reads), dtype=np.uint16)
    for k, arg in enumerate(args):
        read_score_file(arg, np_arr, k)
    return np_arr


def compare_degree_of_nodes(np_arr):
    absolute_number_of_high_quality_edges = np.zeros((len(QUALITY_THRESHOLDS), len(np_arr)), dtype=np.int32)
    relative_number_of_high_quality_edges = np.zeros((len(QUALITY_THRESHOLDS), len(np_arr)), dtype=np.float32)
    for i, elem in enumerate(np_arr[0]):
        for j, threshold in enumerate(QUALITY_THRESHOLDS):
            foo = len(np_arr[0][i][np_arr[0][i] > (np.amax(np_arr[0][i]) * threshold)])
            for k in range(0, len(np_arr)):
                absolute_number_of_high_quality_edges[j][k] += len(
                    np_arr[k][i][np_arr[k][i] > (np.amax(np_arr[0][i]) * threshold)])
                relative_number_of_high_quality_edges[j][k] += len(
                    np_arr[k][i][np_arr[k][i] > (np.amax(np_arr[0][i]) * threshold)]) / foo
    for j, threshold in enumerate(QUALITY_THRESHOLDS):
        relative_number_of_high_quality_edges[j][0] = 1
        for k in range(1, len(np_arr)):
            relative_number_of_high_quality_edges[j][k] /= len(np_arr[0])
    # print(absolute_number_of_high_quality_edges)
    # print(relative_number_of_high_quality_edges)
    return absolute_number_of_high_quality_edges, relative_number_of_high_quality_edges


def write_csv(filename, arr) -> None:
    headers = ['Number', 'Aligner', *QUALITY_THRESHOLDS]
    rows = [[i, aligner, *[arr[j][i] for j, threshold in enumerate(QUALITY_THRESHOLDS)]] for i, aligner in
            enumerate(ALIGNERS)]
    with open(filename, 'w') as f:
        f_csv = csv.writer(f, delimiter='\t')
        f_csv.writerow(headers)
        # f_csv.writerows([[x[0], *x[1]] for x in zip(QUALITY_THRESHOLDS, arr[:])])
        f_csv.writerows(rows)


def write_average_edge_stats(read_length_string='') -> None:
    average_edge_stats = np.zeros((len(np_arr), len(QUALITY_THRESHOLDS)), dtype=float)
    counter = 0
    for file in glob.glob(DIR + '*l{}*/relative_edge_comparison.csv'.format(read_length_string)):
        with open(file, 'r') as f:
            f_csv = csv.reader(f, delimiter='\t')
            headers = next(f_csv)
            for i, row in enumerate(f_csv):
                for j, threshold in enumerate(QUALITY_THRESHOLDS):
                    average_edge_stats[i][j] += float(row[j + 2])
        counter += 1
    average_edge_stats /= counter
    with open(DIR + 'relative_edge_comparison{}.csv'.format(read_length_string), 'w') as f:
        f_csv = csv.writer(f, delimiter='\t')
        f_csv.writerow(headers)
        rows = [[i, aligner, *[average_edge_stats[i][j] for j, threshold in enumerate(QUALITY_THRESHOLDS)]] for
                i, aligner in enumerate(ALIGNERS) if aligner != 'SeqAn']
        f_csv.writerows(rows)


for ref in [1]:
    for coverage in [5, 20, 40]:
        for length in [100, 400, 700]:
            np_arr = read_score_files(
                DIR + 'ref{0}_c{1}_l{2}/seqan.score'.format(ref, coverage, length),
                DIR + 'ref{0}_c{1}_l{2}/seqan_0_max.score'.format(ref, coverage, length),
                DIR + 'ref{0}_c{1}_l{2}/seqan_-5_5.score'.format(ref, coverage, length),
                DIR + 'ref{0}_c{1}_l{2}/seqan_-20_20.score'.format(ref, coverage, length),
                DIR + 'ref{0}_c{1}_l{2}/seqan_0_5.score'.format(ref, coverage, length),
                DIR + 'ref{0}_c{1}_l{2}/seqan_0_20.score'.format(ref, coverage, length),
                DIR + 'ref{0}_c{1}_l{2}/calign.score'.format(ref, coverage, length))
            absolute_number_of_high_quality_edges, relative_number_of_high_quality_edges = compare_degree_of_nodes(
                np_arr)
            write_csv(DIR + 'ref{0}_c{1}_l{2}/absolute_edge_comparison.csv'.format(ref, coverage, length),
                      absolute_number_of_high_quality_edges)
            write_csv(DIR + 'ref{0}_c{1}_l{2}/relative_edge_comparison.csv'.format(ref, coverage, length),
                      relative_number_of_high_quality_edges)
write_average_edge_stats()
write_average_edge_stats('100')
write_average_edge_stats('400')
write_average_edge_stats('700')
