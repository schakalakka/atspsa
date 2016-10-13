__author__ = 'andreas'

import os

from Bio import SeqIO
from Bio.Alphabet import generic_dna


def parse_fasta(fasta_file):
    print("Reading fasta file...")
    return [(str(seqRecord.seq)) for seqRecord in SeqIO.parse(fasta_file, "fasta", generic_dna)]
    print("Reading fasta file finished")


def parse_fasta_with_id(fasta_file):
    print("Reading fasta file...")
    return [(str(seqRecord.id), (str(seqRecord.seq))) for seqRecord in SeqIO.parse(fasta_file, "fasta", generic_dna)]
    print("Reading fasta file finished")


def parse_tour(tour_file, circular=False):
    if circular:
        add_extra_non_circular_node = 0
    else:  # adds an extra node to transform a path to a tour, here for substracting the index from the tour file
        add_extra_non_circular_node = 1
    print("Parsing tour.")
    tour = []
    with open(tour_file) as f:
        [tour.append(int(line.splitlines()[0]) - 1 - add_extra_non_circular_node) for line in f if
         ((line.splitlines()[0].isdigit()) and (int(line.splitlines()[0]) > add_extra_non_circular_node))]
    print("Parsing tour finished.")
    return tour


def run_lkh(LKH_EXE, lkh_par_file):
    os.system("{}LKH {}".format(LKH_EXE, lkh_par_file + ".par"))
