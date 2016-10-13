import os
import random

import matplotlib.pyplot as plt

from config import *


def create_reads(ref):
    with open(ref, 'r') as f:
        myseq = f.read()
        mylength = len(myseq)
        nuc_list = list(myseq)
    for coverage in coverages:
        for average_length in average_length_list:
            nr_of_reads = int(mylength * coverage / average_length)
            dist_till_next_read = int(mylength / nr_of_reads)
            foo = "mkdir {}".format(DIR + ref.split('.')[0] + "shuffled_c{}_l{}".format(coverage, average_length))
            os.system(foo)
            filename = DIR + ref.split('.')[0] + "shuffled_c{0}_l{1}/shuffled_c{0}_l{1}".format(coverage,
                                                                                                average_length)
            pre_list = []
            for i in range(nr_of_reads):
                startindex = random.randint(0, mylength)
                length = average_length + random.randint(-average_length * variance, average_length * variance)
                pre_list.append((startindex, length))
            pre_list = sorted(pre_list, key=lambda tup: tup[0])
            with open("{}.fasta".format(filename, nr_of_reads, average_length), "w") as handle:
                for pre_elem in pre_list:
                    startindex = pre_elem[0]
                    length = pre_elem[1]
                    endindex = int(min(startindex + length, mylength))
                    if endindex - startindex >= 30:
                        currentseq = nuc_list[startindex:endindex]
                        handle.write(
                            ">NC_005816_Yersinia_pestis_biovar_Microtus/{}/{}_{}\n".format(i + 1, startindex, endindex))
                        handle.write("".join(currentseq) + "\n")
                        # os.system("{}fasta2DB -v {} {}.fasta".format(DAZZ_DB, filename, filename))


def make_latex_stat_table(stats, ref):
    output = '\\begin{tabular}{|c|c|c|c|c|c|}\\hline\n& \\multicolumn{4}{c|}{Read Length} & \\\\\nCoverage & Average & Median & Shortest' \
             '& Longest & \\#Reads \\\\\n\\hline\n'
    for coverage in coverages:
        for average_length in average_length_list:
            val = stats[(coverage, average_length)]
            output = output + '{0:.1f} & {1:.1f} & {2} & {3} & {4} & {5} \\\\\n'.format(val[0], val[2], val[3],
                                                                                        val[4], val[5], val[1])
        output = output + '\\hline\n'
    output = output + '\\end{tabular}'
    with open(DIR + ref.split('.')[0] + '_stats.tex', 'w') as f:
        f.write(output)
    print(output)


def plot(prelist, filename, coverage, average_length):
    # plt.figure()
    plt.title("Coverage: {}, Average Length: {}, Reads: {}\nReads below average: {}, Reads above average: {}".format(
        coverage, average_length, len(prelist), len([x for x in prelist if x[1] < average_length]),
        len([x for x in prelist if x[1] > average_length])))
    # plt.hist([tup[1] for tup in prelist], bins=[0,average_length, maximum_length*average_length])
    plt.hist([tup[1] for tup in prelist], bins=50)
    plt.grid(True)
    # plt.xticks(xrange(len(list)))
    # plt.yticks(xrange(0,max(nr_best_edge_used)+10, 5))
    plt.ylabel("Frequency")
    plt.xlabel("Read Length")
    plt.savefig(filename + "_plot.pdf")
    plt.close()


def create_reads_shotgun_style(ref, ref_nr):
    '''
    create reads shotgun sequencing style
    :param ref:
    :return:
    '''
    stats = {}
    with open(ref, 'r') as f:
        ref_seq = f.read()
        ref_length = len(ref_seq)
        nuc_list = list(ref_seq)
    for coverage in coverages:
        for average_length in average_length_list:
            nr_of_shots = int(ref_length / average_length)
            os.system("mkdir {}".format(DIR + ref.split('.')[0] + "shuffled_c{}_l{}".format(coverage, average_length)))
            filename = DIR + ref.split('.')[0] + "shuffled_c{0}_l{1}/shuffled_c{0}_l{1}".format(coverage,
                                                                                                average_length)
            pre_list = []
            while sum([tup[1] for tup in pre_list]) < coverage * ref_length:
                hits = [random.randint(0, 20)]  # stores the hits of the shotgun sequence in on run
                for i in range(nr_of_shots):
                    hits.append(random.randint(0, ref_length))
                hits = sorted(hits)
                foolist = [(hits[i], hits[i + 1] - hits[i]) for i in range(len(hits) - 1)
                           if
                           minimum_length * average_length <= hits[i + 1] - hits[i] <= maximum_length * average_length]
                pre_list = pre_list + foolist
            pre_list = sorted(pre_list, key=lambda tup: tup[0])
            with open("{}.stat".format(filename), "w") as handle:
                handle.write('Real Coverage: {}\n'
                             'Number of reads: {}\n'
                             'Average read length: {}\n'
                             'Median read length: {}\n'
                             'Shortest read length" : {}\n'
                             'Longest read length: {}'.format(
                    sum([tup[1] for tup in pre_list]) / ref_length,
                    len(pre_list),
                    sum([tup[1] for tup in pre_list]) / len(pre_list),
                    sorted([tup[1] for tup in pre_list])[len(pre_list) // 2],
                    min([tup[1] for tup in pre_list]),
                    max([tup[1] for tup in pre_list])
                ))
                # save stats in dict for later creation of latex table
                stats[(coverage, average_length)] = (
                    sum([tup[1] for tup in pre_list]) / ref_length,
                    len(pre_list),
                    sum([tup[1] for tup in pre_list]) / len(pre_list),
                    sorted([tup[1] for tup in pre_list])[len(pre_list) // 2],
                    min([tup[1] for tup in pre_list]),
                    max([tup[1] for tup in pre_list])
                )
            with open("{}.fasta".format(filename), "w") as handle:
                for i, pre_elem in enumerate(pre_list):
                    startindex = pre_elem[0]
                    length = pre_elem[1]
                    endindex = int(min(startindex + length, ref_length))
                    if endindex - startindex >= 30:
                        currentseq = nuc_list[startindex:endindex]
                        handle.write(
                            ">{}/{}/{}_{}\n".format(references[ref_nr].replace(' ', '-'), i, startindex, endindex))
                        handle.write("".join(currentseq) + "\n")
                        # os.system("{}fasta2DB -v {} {}.fasta".format(DAZZ_DB, filename, filename))
            plot(pre_list, filename, coverage, average_length)
    make_latex_stat_table(stats, ref)


create_reads_shotgun_style(ref1, 0)
create_reads_shotgun_style(ref2, 1)
create_reads_shotgun_style(ref3, 2)
