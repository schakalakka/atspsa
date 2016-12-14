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
            foo = "mkdir {}".format(DIR + ref.split('.')[0] + "_c{}_l{}".format(coverage, average_length))
            os.system(foo)
            filename = DIR + ref.split('.')[0] + "_c{0}_l{1}/reads.fasta".format(coverage,
                                                                                 average_length)
            pre_list = []
            for i in range(nr_of_reads):
                startindex = random.randint(0, mylength)
                length = average_length + random.randint(-average_length * variance, average_length * variance)
                pre_list.append((startindex, length))
            pre_list = sorted(pre_list, key=lambda tup: tup[0])
            # write fasta
            with open(filename, "w") as handle:
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


def plot(prelist, sub_dir, coverage, average_length):
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
    plt.savefig(DIR + sub_dir + "plot.pdf")
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
            sub_dir = ref.split('.')[0] + "ref{0}_c{1}_l{2}/".format(ref_nr + 1, coverage, average_length)
            os.system("mkdir {}".format(DIR + sub_dir))
            meta_read_list = []  # list indicating the start and the length of a read -> (start_index,length)

            while sum(x[1] for x in meta_read_list) < coverage * ref_length:
                start_index = random.randint(0, ref_length)
                length = random.randint(average_length - average_length * variance,
                                        average_length + average_length * variance)
                if start_index + length <= ref_length and (start_index, length) not in meta_read_list:
                    meta_read_list.append((start_index, length))
            # sort by start index
            meta_read_list = sorted(meta_read_list, key=lambda x: -x[1])
            meta_read_list = sorted(meta_read_list, key=lambda x: x[0])
            # write fasta file
            with open(DIR + sub_dir + "reads.fasta", "w") as f:
                for i, elem in enumerate(meta_read_list):
                    startindex = elem[0]
                    length = elem[1]
                    endindex = startindex + length
                    currentseq = nuc_list[startindex:endindex]
                    f.write(
                        '>{}/{}/{}_{}_{}\n'.format(references[ref_nr].replace(' ', '-'), i, startindex, endindex,
                                                   length))
                    f.write(''.join(currentseq) + '\n')

            # write fasta stats
                    # with open(DIR + sub_dir + "fasta.stat", "w") as f:
            #     real_coverage = sum([x[1] for x in meta_read_list]) / ref_length
            #     nr_reads = len(meta_read_list)
            #     average = sum([x[1] for x in meta_read_list]) / len(meta_read_list)
            #     median = sorted([x[1] for x in meta_read_list])[len(meta_read_list) // 2]
            #     shortest = min([x[1] for x in meta_read_list])
            #     longest = max([x[1] for x in meta_read_list])
            #     target_value = compute_target_function_value_of_fasta(DIR+sub_dir+'reads.fasta')
                    #     f.write('Real Coverage: {}\n'
            #                  'Number of reads: {}\n'
            #                  'Average read length: {}\n'
            #                  'Median read length: {}\n'
            #                  'Shortest read length" : {}\n'
            #                  'Longest read length: {}\n'
            #                  'Target function value: {}'.format(real_coverage, nr_reads, average, median, shortest, longest, target_value))
            #     # save stats in dict for later creation of latex table
            #     stats[(coverage, average_length)] = (real_coverage, nr_reads, average, median, shortest, longest, target_value)

                    # plot(meta_read_list, sub_dir, coverage, average_length)
            # make_latex_stat_table(stats, ref)

# create_reads_shotgun_style(ref1, 0)
# create_reads_shotgun_style(ref2, 1)
# create_reads_shotgun_style(ref3, 2)
