import os
import random

from config import *

mylength = len(myseq)
nuc_list = list(myseq)


def create_pure_random_data():
    for coverage in coverages:
        for average_length in average_length_list:
            nr_of_reads = int(mylength * coverage / average_length)
            dist_till_next_read = int(mylength / nr_of_reads)

            foo = "mkdir {}".format(DIR + "shuffled_c{}_l{}".format(coverage, average_length))
            os.system(foo)
            filename = DIR + "shuffled_c{}_l{}/shuffled_c{}_l{}".format(coverage, average_length, coverage,
                                                                        average_length)
            with open("{}.fasta".format(filename, nr_of_reads, average_length), "w") as handle:
                for i in range(nr_of_reads):
                    while True:
                        startindex = random.randint(0, mylength)
                        length = average_length + random.randint(-average_length * 0.1, average_length * 0.1)
                        endindex = int(min(startindex + length, mylength))
                        if endindex - startindex >= 30:
                            break
                    currentseq = nuc_list[startindex:endindex]
                    handle.write(">NC_005816_Yersinia_pestis_biovar_Microtus/{}/0_{}\n".format(i + 1, len(currentseq)))
                    handle.write("".join(currentseq) + "\n")
            os.system("{}fasta2DB -v {} {}.fasta".format(DAZZ_DB, filename, filename))


def create_obvious_initial_data():
    for coverage in coverages:
        for average_length in average_length_list:
            nr_of_reads = int(mylength * coverage / average_length)
            dist_till_next_read = int(mylength / nr_of_reads)
            foo = "mkdir {}".format(DIR + "shuffled_c{}_l{}_n{}".format(coverage, average_length, nr_of_reads))
            os.system(foo)
            filename = DIR + "shuffled_c{}_l{}_n{}/shuffled_c{}_l{}_n{}".format(coverage, average_length, nr_of_reads,
                                                                                coverage, average_length, nr_of_reads)
            with open("{}.fasta".format(filename, nr_of_reads, average_length), "w") as handle:
                for i in range(nr_of_reads):
                    startindex = i * dist_till_next_read
                    length = average_length  # + random.randint(-average_length*0.1,average_length*0.1)
                    endindex = int(min(startindex + length, mylength))
                    if endindex - startindex >= 30:
                        currentseq = nuc_list[startindex:endindex]
                        handle.write(
                                ">NC_005816_Yersinia_pestis_biovar_Microtus/{}/{}_{}\n".format(i + 1, startindex,
                                                                                               endindex))
                        handle.write("".join(currentseq) + "\n")
            os.system("{}fasta2DB -v {} {}.fasta".format(DAZZ_DB, filename, filename))


def create_initial_data():
    for coverage in coverages:
        for average_length in average_length_list:
            nr_of_reads = int(mylength * coverage / average_length)

            dist_till_next_read = int(mylength / nr_of_reads)
            foo = "mkdir {}".format(DIR + "shuffled_c{}_l{}".format(coverage, average_length))
            os.system(foo)
            filename = DIR + "shuffled_c{}_l{}/shuffled_c{}_l{}".format(coverage, average_length, coverage,
                                                                        average_length)
            pre_list = []
            for i in range(nr_of_reads):
                startindex = random.randint(0, mylength)
                length = average_length + random.randint(-average_length * 0.1, average_length * 0.1)
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
            os.system("{}fasta2DB -v {} {}.fasta".format(DAZZ_DB, filename, filename))


create_pure_random_data()
# create_initial_data()
