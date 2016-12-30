import csv
import glob
import os
from typing import Tuple

from config import *


def read_alignment_time(time_file: str) -> float:
    if os.path.isfile(time_file):
        with open(time_file, 'r') as f:
            return float(f.read())
    else:
        return None


def read_assembly_times(time_file: str) -> Tuple[float, float]:
    if os.path.isfile(time_file):
        with open(time_file, 'r') as f:
            foo = f.read()
            lkh_time = float(foo.split('LKH_Time: ')[1].split('\n')[0])
            ap_time = float(foo.split('AP_Time: ')[1].split('\n')[0])
            if lkh_time > -1.0 and ap_time > -1.0:
                return lkh_time, ap_time

    return None, None


def write_time_csv(alignment=True):
    reads_to_times = {}
    for file in glob.glob(DIR + 'ref*00/fasta.stat'):
        subdir = file[:-10]
        length = int(subdir[-4:-1])
        with open(file, 'r') as f:
            nr_of_reads = int(f.read().split('Number of reads: ')[1].split('\n')[0])

            calign_time = read_alignment_time(subdir + 'calign.time')
            atsp_time, ap_time = read_assembly_times(subdir + 'calign.assembly')

            calign_time25 = read_alignment_time(subdir + 'calign_0_{}.time'.format(length // 4))
            atsp_time25, ap_time25 = read_assembly_times(subdir + 'calign_0_{}.assembly'.format(length // 4))

            calign_time50 = read_alignment_time(subdir + 'calign_0_{}.time'.format(length // 2))
            atsp_time50, ap_time50 = read_assembly_times(subdir + 'calign_0_{}.assembly'.format(length // 2))

            if alignment:
                if atsp_time:
                    atsp_time += calign_time
                if ap_time:
                    ap_time += calign_time
                if atsp_time25:
                    atsp_time25 += calign_time25
                if ap_time25:
                    ap_time25 += calign_time25
                if atsp_time50:
                    atsp_time50 += calign_time50
                if ap_time50:
                    ap_time50 += calign_time50

            if atsp_time and ap_time and calign_time and atsp_time25 and ap_time25 and calign_time25 and atsp_time50 and ap_time50 and calign_time50:
                reads_to_times[(nr_of_reads, length)] = {}
                reads_to_times[(nr_of_reads, length)]['catsp'] = atsp_time
                reads_to_times[(nr_of_reads, length)]['cap'] = ap_time
                reads_to_times[(nr_of_reads, length)]['catspquarter'] = atsp_time25
                reads_to_times[(nr_of_reads, length)]['capquarter'] = ap_time25
                reads_to_times[(nr_of_reads, length)]['catsphalf'] = atsp_time50
                reads_to_times[(nr_of_reads, length)]['caphalf'] = ap_time50

    csv_file = 'summed_time.csv' if alignment else 'assembly_time.csv'
    with open(DIR + csv_file, 'w') as f:
        f_csv = csv.writer(f, delimiter='\t')
        f_csv.writerow(['reads', 'catsp', 'cap', 'catspquarter', 'capquarter', 'catsphalf', 'caphalf'])
        for key, d in sorted(reads_to_times.items()):
            nr_of_reads, length = key[0], key[1]
            f_csv.writerow(
                [nr_of_reads, d['catsp'], d['cap'], d['catspquarter'], d['capquarter'], d['catsphalf'], d['caphalf']])


write_time_csv(True)
write_time_csv(False)
