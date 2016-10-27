import csv
import glob
from collections import namedtuple
from typing import List, Tuple

from config import *

HEADER = ['File', 'LKHContigs', 'LKHValue', 'LKHTime', 'APContigs', 'APValue', 'APTime', 'ActualObjectiveValue']
Assembly_Stats = namedtuple('Assembly_Stats', HEADER)

file = '/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c5_l700/calign.assembly'


def read_assembly_file(file: str) -> Tuple[List[str], float, List[str], float]:
    with open(file, 'r') as f:
        file_content_string = f.read()
        lkh_contigs = file_content_string.split('LKH_Contigs:\n')[1].split('\nLKH_Objective')[0].split('\n')
        lkh_value = int(file_content_string.split('LKH_Objective_Value: ')[1].split('\n')[0])
        lkh_time = float(file_content_string.split('LKH_Time: ')[1].split('\n')[0])
        ap_contigs = file_content_string.split('AP_Contigs:\n')[1].split('\nAP_Objective')[0].split('\n')
        ap_value = int(file_content_string.split('AP_Objective_Value: ')[1].split('\n')[0])
        ap_time = float(file_content_string.split('AP_Time: ')[1].split('\n')[0])
    with open(file.split('calign.assembly')[0] + 'fasta.stat', 'r') as f:
        file_content_string = f.read()
        actual_Objective_value = int(file_content_string.split('Objective function value: ')[1].split('\n')[0])
    return lkh_contigs, lkh_value, lkh_time, ap_contigs, ap_value, ap_time, actual_Objective_value


def write_assembly_stats(assembly_stats_list: List[Assembly_Stats]) -> None:
    with open('/home/andreas/GDrive/workspace/sparsedata/assembly_stats.csv', 'w') as f:
        f_csv = csv.writer(f, delimiter=',')
        f_csv.writerow(
            ['File', 'LKHContigs', 'LKHValue', 'LKHTime', 'APContigs', 'APValue', 'APTime', 'ActualObjectiveValue'])
        for elem in assembly_stats_list:
            f_csv.writerow(elem)


assembly_stats_list = []
for file in sorted(glob.glob('/home/andreas/GDrive/workspace/sparsedata/ref[1,2]_c*/calign.assembly')):
    lkh_contigs, lkh_value, lkh_time, ap_contigs, ap_value, ap_time, actual_Objective_value = read_assembly_file(file)
    file_sub_dir = file.split('/')[-2]  # example ref1_c5_l100
    ref_number = int(file_sub_dir.split('ref')[1].split('_')[0])
    coverage = int(file_sub_dir.split('_c')[1].split('_')[0])
    length = int(file_sub_dir.split('_l')[1])
    file = '{}-{}-{}'.format(references[ref_number - 1], coverage, length)
    assembly_stats_list.append(
        Assembly_Stats(file, len(lkh_contigs), lkh_value, lkh_time, len(ap_contigs), ap_value, ap_time,
                       actual_Objective_value))

write_assembly_stats(assembly_stats_list)
