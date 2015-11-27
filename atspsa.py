__author__ = 'andreas'

import createtsp
import db2score
import os
from Bio import SeqIO
from Bio.Alphabet import generic_dna

#configuartion

MHAP_JAR = "/home/andreas/mhap/mhap-1.6.jar"
LKH_EXE = "/home/andreas/lkh/LKH-2.0.7/LKH"
DAZZ_DB = "/home/andreas/PycharmProject/DAZZ_DB"
DALIGNER = "/home/andreas/PycharmProjects/DALIGNER/"


def parse_fasta(fasta_file):

    return [(str(seqRecord.seq)) for seqRecord in SeqIO.parse(file, "fasta", generic_dna)]


def parse_tour(tour_file):
    tour = []
    with open(file) as f:
        [tour.append(int(line.splitlines()[0]) - 2) for line in f if
         ((line.splitlines()[0].isdigit()) and (int(line.splitlines()[0]) > 1))]

    return tour


def run_lkh(lkh_par_file):
    os.system("{} {}".format(LKH_EXE, lkh_par_file))


def main():
    filename = "/home/andy/seqdata/ecoli.1.subreads/ecoli.1.subreads"

    DB = filename + ".db"
    FASTA = filename + ".fasta"
    LAS = filename + ".las"
    OVL = filename + ".ovl"
    LKH_PAR = filename + ".par"
    LKH_OUT = filename + ".tour"
    LKH_LIB = filename + ".atsp"

    reads = parse_fasta(FASTA)
    # mhapper.run_mhap(MHAP_JAR, FASTA_FILE, MHAP_OUT)
    # scores = mhapper.parse_mhap(MHAP_OUT)
    scores = db2score.read_ovl_file(OVL)
    createtsp.prepare_lkh(LKH_PAR, LKH_LIB, LKH_OUT, len(reads), db2score.get_scores_without_orientation(scores))
    # run_lkh()


if __name__ == "__main__":
    main()
