import makestat
import matplotlib.pyplot as plt

# matplotlib.rcParams['backend'] = "PDF"

import os
import sys

sys.path.insert(0, '/home/andreas/PycharmProjects/atspsa')

import stats


def plot(nr_of_reads, average_length, nr_best_edge_used):
    filename = "shuffled_{}_{}".format(nr_of_reads, average_length)
    plt.figure()
    plt.plot(nr_best_edge_used, 'k', nr_best_edge_used, "bo")
    plt.title(filename)
    plt.grid(True)
    plt.xticks(xrange(len(nr_best_edge_used)))
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
            plt.xticks(xrange(len(nr_best_edge_used)))
            plt.yticks(xrange(0, max(nr_best_edge_used) + 10, 5))
            plt.ylabel("How often is the n_th best edge used")
            plt.xlabel("n-th best edge")
            plt.savefig(filename + "_plot.pdf")
            # plt.show()


def plot_all_with_sum():
    mydict = makestat.makestats()
    DIR = "/run/media/andreas/INTENSO/fastas/artdata/"
    curr_DIR = DIR + "plots/"
    os.chdir(curr_DIR)
    for (nr_of_reads, average_length), nr_best_edge_used in mydict.items():
        plot(nr_of_reads, average_length, nr_best_edge_used)


plot_all_with_sum()
