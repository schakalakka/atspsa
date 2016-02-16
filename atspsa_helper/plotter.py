import os
import pickle

import matplotlib.pyplot as plt
import numpy as np
import producestats

from config import *
# matplotlib.rcParams['backend'] = "PDF"



from atspsacore import stats


def plot(nr_of_reads, average_length, nr_best_edge_used):
    filename = "shuffled_{}_{}".format(nr_of_reads, average_length)
    plt.figure()
    plt.plot(nr_best_edge_used, 'k', nr_best_edge_used, "bo")
    plt.title(filename)
    plt.grid(True)
    plt.xticks(range(len(nr_best_edge_used)))
    # plt.yticks(xrange(0,max(nr_best_edge_used)+10, 5))
    plt.ylabel("How often is the n_th best edge used")
    plt.xlabel("n-th best edge")
    plt.savefig(filename + "_plot.pdf")


def plot_all_without_sum():
    DIR = "/run/media/andreas/INTENSO/fastas/artdata/"
    # nr_of_reads_list = [100]
    # average_length_list = [300]
    nr_of_reads_list = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    average_length_list = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    for nr_of_reads in nr_of_reads_list:
        for average_length in average_length_list:
            curr_DIR = DIR + "shuffled_{}_{}/".format(nr_of_reads, average_length)
            os.chdir(curr_DIR)
            filename = "shuffled_{}_{}".format(nr_of_reads, average_length)
            (nr_best_edge_used, scores) = stats.make_stat(filename)
            plt.plot(nr_best_edge_used, 'k', nr_best_edge_used, "bo")
            plt.title(filename)
            plt.grid(True)
            plt.xticks(range(len(nr_best_edge_used)))
            plt.yticks(range(0, max(nr_best_edge_used) + 10, 5))
            plt.ylabel("How often is the n_th best edge used")
            plt.xlabel("n-th best edge")
            plt.savefig(filename + "_plot.pdf")
            # plt.show()


def plot_all_with_sum():
    mydict = producestats.makestats()
    DIR = "/run/media/andreas/INTENSO/fastas/artdata/"
    curr_DIR = DIR + "plots/"
    os.chdir(curr_DIR)
    for (nr_of_reads, average_length), nr_best_edge_used in mydict.items():
        plot(nr_of_reads, average_length, nr_best_edge_used)


def make_align_matrix_plot():
    dirs = os.listdir(DIR)
    for current_elem in dirs:
        curr_DIR = DIR + current_elem
        os.chdir(curr_DIR)
        filename = current_elem + "_alignstat.pickle"
        with open(filename, 'rb') as f:
            dictmat = pickle.load(f)

        # get dimensions of matrix
        bar = max(max([key[0] for key in dictmat.keys()]), max([key[1] for key in dictmat.keys()]))
        mat = np.zeros((bar + 1, bar + 1))

        # fill matrix
        for key, value in dictmat.items():
            mat[key[0]][key[1]] = value

        fig = plt.figure()
        ax = fig.add_subplot(111)
        # plt.matshow(mat)
        cax = ax.matshow(mat)
        fig.colorbar(cax)

        plt.title(current_elem + ' max SW score')
        plt.grid(True)
        plt.ylabel("Sequence 2")
        plt.xlabel("Sequence 1")
        plt.savefig(current_elem + "_SW_plot.pdf")
        #plt.show()


def make_graph_matrix_plot():
    dirs = os.listdir(DIR)
    for current_elem in dirs:
        curr_DIR = DIR + current_elem
        os.chdir(curr_DIR)
        filename = current_elem + ".alignpickle"
        with open(filename, 'rb') as f:
            scores = pickle.load(f)

        # get dimensions of matrix
        bar = max(max([key[0] for key in scores.keys()]), max([key[1] for key in scores.keys()]))
        mat = np.zeros((bar + 1, bar + 1))

        # fill matrix
        for key, value in scores.items():
            mat[key[0]][key[1]] = value

        fig = plt.figure()
        ax = fig.add_subplot(111)
        # plt.matshow(mat)
        cax = ax.matshow(mat)
        fig.colorbar(cax)

        plt.title(current_elem + ' TSP graph ')
        plt.grid(True)
        # plt.xticks(range(len(nr_best_edge_used)))
        # plt.yticks(range(0, max(nr_best_edge_used) + 10, 5))
        plt.ylabel("Edges")
        plt.xlabel("Edges")
        plt.savefig(current_elem + "_graph_plot.pdf")
        #plt.show()


# plot_all_with_sum()
make_align_matrix_plot()
make_graph_matrix_plot()
