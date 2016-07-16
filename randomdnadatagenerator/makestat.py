import os
import sys

from config import *

sys.path.insert(0, ATSPSA_DIR)
import stats


def makestats():
    # nr_of_reads_list = [100]
    # average_length_list = [100]
    nr_of_reads_list = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
    average_length_list = [100, 200, 300, 400, 500, 600, 700, 800, 900]

    mydict = {}

    with open("/run/media/andreas/INTENSO/fastas/artdata/stats.csv", "w") as f:
        f.write("nr_of_reads; average_length; nr_of_edges; best_edge_list \n")
        for nr_of_reads in nr_of_reads_list:
            for average_length in average_length_list:
                curr_DIR = DIR + "shuffled_{}_{}/".format(nr_of_reads, average_length)
                os.chdir(curr_DIR)
                filename = "shuffled_{}_{}".format(nr_of_reads, average_length)
                (nr_best_edge_used, scores) = stats.make_stat(filename)
                f.write("{}; {}; {}; {}\n".format(nr_of_reads, average_length, len(scores), nr_best_edge_used))
                mydict[(nr_of_reads, average_length)] = nr_best_edge_used

    mydict2 = {}
    mydict2[(0, 0)] = [sum(x) for x in zip(*[y for key, y in mydict.items()])]
    for average_length in average_length_list:
        mydict2[(0, average_length)] = [sum(x) for x in zip(*[mydict[(x, average_length)] for x in nr_of_reads_list])]

    for nr_of_reads in nr_of_reads_list:
        mydict2[(nr_of_reads, 0)] = [sum(x) for x in zip(*[mydict[(nr_of_reads, x)] for x in average_length_list])]

    mydict.update(mydict2)
    return mydict
