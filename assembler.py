__author__ = 'andreas'


from Bio import SeqIO
from Bio.Alphabet import generic_dna
#import nwalign as nw
#import swalign
import os
import math

#configuartion
FILE_NAME = "SRR1536433_1000"
DIR = "/run/media/andreas/INTENSO/fastas/"
MHAP_JAR = "/home/andreas/mhap/mhap-1.6.jar"
LKH_EXE = "/home/andreas/lkh/LKH-2.0.7/LKH"
FASTA_FILE = "{}fastas/{}.fasta".format(DIR, FILE_NAME)
MHAP_OUT = "{}mhaps/{}.out".format(DIR, FILE_NAME)
LKH_PAR = "{}lkhs/{}.par".format(DIR, FILE_NAME)
LKH_OUT = "{}lkhs/{}.tour".format(DIR, FILE_NAME)
LKH_LIB = "{}lkhs/{}.atsp".format(DIR, FILE_NAME)

def parse_fasta(file=None):
    if file is None:
        file = FASTA_FILE
    return [(str(seqRecord.seq)) for seqRecord in SeqIO.parse(file, "fasta", generic_dna)]

def parse_mhap(file=None):
    if file is None:
        file = MHAP_OUT
    dd = {}
    #read mhap output
    with open(file) as f:
        for line in f:
            if line.startswith("#"):
                break
            else:
                line_arr = line.split(" ")
                i = int(line_arr[0])-1
                j = int(line_arr[1])-1
                dd[(i,j)] = int(float(line_arr[2])) #mhap counts from 1 to nr_of_reads
    return dd



def run_mhap():
    #run mhap
    os.system("java -server -Xmx32g -jar {} -s {} > {}".format(MHAP_JAR, FASTA_FILE, MHAP_OUT))

def prepare_lkh(reads, mhap_overlap):
    tsplib_par_string = "PROBLEM_FILE={}\nOUTPUT_TOUR_FILE={}".format(LKH_LIB, LKH_OUT)
    with open(LKH_PAR, "w") as f:
        f.write(tsplib_par_string)

    with open(LKH_LIB, "w") as f:
        # dimension +1 because of the special knot
        tsplib_string = "NAME: {}\nTYPE: ATSP \nCOMMENT: {}\nDIMENSION: {} \nEDGE_WEIGHT_TYPE: EXPLICIT \n" \
                        "EDGE_WEIGHT_FORMAT: FULL_MATRIX \nEDGE_WEIGHT_SECTION".format(FILE_NAME, FILE_NAME,
                                                                                       (len(reads) + 1))

        f.write(tsplib_string)

        # first line for the extra knot
        f.write('\n' + '\t'.join('0' for _ in range(len(reads) + 1)))

        for i in range(len(reads)):
            f.write('\n'+'\t'.join([str(-mhap_overlap.get((i,j), 0)) for j in range(-1, len(reads))]))

        f.write("\nEOF")






def run_lkh():
    os.system("{} {}".format(LKH_EXE, LKH_PAR))

def set_config(file=None, directory=None, jar=None, lkh=None):
    if file is not None:
        global FILE_NAME
        FILE_NAME = file

    if directory is not None:
        global DIR
        DIR = directory

    FASTA_FILE = "{}fastas/{}.fasta".format(DIR, FILE_NAME)
    MHAP_OUT = "{}mhaps/{}.out".format(DIR, FILE_NAME)
    LKH_PAR = "{}lkhs/{}.par".format(DIR, FILE_NAME)
    LKH_OUT = "{}lkhs/{}.tour".format(DIR, FILE_NAME)
    LKH_LIB = "{}lkhs/{}.atsp".format(DIR, FILE_NAME)

    if lkh is not None:
        global LKH_EXE
        LKH_EXE = lkh
    if jar is not None:
        global MHAP_JAR
        MHAP_JAR = jar
