__author__ = 'andreas'

import createtsp
import db2score
import os
from Bio import SeqIO
from Bio.Alphabet import generic_dna

#configuartion

MHAP_JAR = None
LKH_EXE = None
DAZZ_DB = None
DALIGNER = None

def parse_fasta(fasta_file):
    print("Reading fasta file...")
    return [(str(seqRecord.seq)) for seqRecord in SeqIO.parse(fasta_file + ".fasta", "fasta", generic_dna)]
    print("Reading fasta file finished")


def parse_tour(tour_file):
    print("Parsing tour.")
    tour = []
    with open(tour_file + ".tour") as f:
        [tour.append(int(line.splitlines()[0]) - 2) for line in f if
         ((line.splitlines()[0].isdigit()) and (int(line.splitlines()[0]) > 1))]
    print("Parsing tour finished.")
    return tour


def run_lkh(lkh_par_file):
    os.system("{}LKH {}".format(LKH_EXE, lkh_par_file + ".par"))


def read_config():
    print("Read config.txt")
    config = {}
    with open("config.txt", "r") as f:
        for line in f:
            config[line.split("=")[0]] = line.split("=")[1].strip("\n")
    config["FILEPATH"] = config["DIR"] + config["FILENAME"]
    # if config["FILENAME"] is not "":
    #     filename = config["FILENAME"]
    #     config["DB"] = config["DIR"] + filename + ".db"
    #     config["FASTA"] = config["DIR"] + filename + ".fasta"
    #     config["LAS"] = config["DIR"] + filename + ".las"
    #     config["OVL"] = config["DIR"] + filename + ".ovl"
    #     config["LKH_PAR"] = config["DIR"] + filename + ".par"
    #     config["LKH_OUT"] = config["DIR"] + filename + ".tour"
    #     config["LKH_LIB"] = config["DIR"] + filename + ".atsp"
    return config


def global_config(config):
    global MHAP_JAR, LKH_EXE, DAZZ_DB, DALIGNER
    MHAP_JAR = config["MHAP_JAR"]
    LKH_EXE = config["LKH_EXE"]
    DAZZ_DB = config["DAZZ_DB"]
    DALIGNER = config["DALIGNER"]

def main():
    config = read_config()
    global_config(config)

    reads = parse_fasta(config["FILEPATH"])
    # mhapper.run_mhap(MHAP_JAR, FASTA_FILE, MHAP_OUT)
    # scores = mhapper.parse_mhap(MHAP_OUT)
    # makealignmentdb.align(config["DALIGNER"], config["DIR"], config["FILEPATH"])
    # makealignmentdb.create_ovl(config["DALIGNER"], config["FILEPATH"])
    scores = db2score.read_ovl_file(config["FILEPATH"])
    createtsp.prepare_lkh(config["FILEPATH"], len(reads),
                          db2score.get_scores_without_orientation(scores))
    run_lkh(config["FILEPATH"])


if __name__ == "__main__":
    main()
