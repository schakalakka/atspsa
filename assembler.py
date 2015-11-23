__author__ = 'andreas'

import createtsp
import db2score
import os
from Bio import SeqIO
from Bio.Alphabet import generic_dna

#configuartion
FILE_NAME = "ecoli_1000"
DIR = "/run/media/andreas/INTENSO/fastas/ecoli_1000"
MHAP_JAR = "/home/andreas/mhap/mhap-1.6.jar"
LKH_EXE = "/home/andreas/lkh/LKH-2.0.7/LKH"
FASTA_FILE = "{}/{}.fasta".format(DIR, FILE_NAME)
MHAP_OUT = "{}/{}.mhap.out".format(DIR, FILE_NAME)
LKH_PAR = "{}/{}.par".format(DIR, FILE_NAME)
LKH_OUT = "{}/{}.tour".format(DIR, FILE_NAME)
LKH_LIB = "{}/{}.atsp".format(DIR, FILE_NAME)
MUNKRES_OUT = "{}/{}.munkres.out".format(DIR, FILE_NAME)
OVL_FILE = "{}/{}.db.ovl".format(DIR, FILE_NAME)

def parse_fasta(file=None):
    if file is None:
        file = FASTA_FILE
    return [(str(seqRecord.seq)) for seqRecord in SeqIO.parse(file, "fasta", generic_dna)]


def parse_tour(file=None):
    if file is None:
        file = LKH_OUT
    tour = []
    with open(file) as f:
        [tour.append(int(line.splitlines()[0]) - 2) for line in f if
         ((line.splitlines()[0].isdigit()) and (int(line.splitlines()[0]) > 1))]

    return tour

def run_lkh():
    os.system("{} {}".format(LKH_EXE, LKH_PAR))


def main():
    reads = parse_fasta(FASTA_FILE)
    # mhapper.run_mhap(MHAP_JAR, FASTA_FILE, MHAP_OUT)
    # scores = mhapper.parse_mhap(MHAP_OUT)
    scores = db2score.read_ovl_file(OVL_FILE)
    createtsp.prepare_lkh(LKH_PAR, LKH_LIB, LKH_OUT, len(reads), db2score.get_scores_without_orientation(scores))
    run_lkh()


if __name__ == "__main__":
    main()
