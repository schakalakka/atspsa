__author__ = 'andreas'

import os

import numpy as np
from Bio import SeqIO
from Bio.Alphabet import generic_dna


def parse_fasta(fasta_file):
    print("Reading fasta file...")
    return [(str(seqRecord.seq)) for seqRecord in SeqIO.parse(fasta_file, "fasta", generic_dna)]


def parse_fasta_with_id(fasta_file):
    print("Reading fasta file...")
    return np.array(
        [(str(seqRecord.id), (str(seqRecord.seq))) for seqRecord in SeqIO.parse(fasta_file, "fasta", generic_dna)])


def parse_tour(tour_file):
    print("Parsing tour.")
    tour = []
    with open(tour_file) as f:
        [tour.append(int(line.splitlines()[0]) - 1) for line in f if line.splitlines()[0].isdigit()]
    print("Parsing tour finished.")
    tour.remove(max(tour))
    return tour


def run_lkh(LKH_EXE, lkh_par_file):
    os.system("{}LKH {}".format(LKH_EXE, lkh_par_file + ".par"))
