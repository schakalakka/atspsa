import csv
import os
from collections import namedtuple
from typing import List, Dict

from config import *

HEADER = ['File', 'LKHContigs', 'LKHValue', 'LKHTime', 'APContigs', 'APValue', 'APTime', 'ActualObjectiveValue']
Assembly_Stats = namedtuple('Assembly_Stats', HEADER)

dir = '/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c5_l700/calign.assembly'


def read_assembly_file(file: str) -> List:
    if not os.path.isfile(file):
        return [-1, -1, -1, -1, -1, -1]
    with open(file, 'r') as f:
        file_content_string = f.read()
        if 'LKH_Contigs:\nLKH_Objective' in file_content_string:
            lkh_gaps = -1
        else:
            lkh_gaps = len(file_content_string.split('LKH_Contigs:\n')[1].split('\nLKH_Objective')[0].split('\n')) - 1
        lkh_value = int(file_content_string.split('LKH_Objective_Value: ')[1].split('\n')[0])
        lkh_time = float(file_content_string.split('LKH_Time: ')[1].split('\n')[0])
        if 'AP_Contigs:\nAP_Objective' in file_content_string:
            ap_gaps = -1
        else:
            ap_gaps = len(file_content_string.split('AP_Contigs:\n')[1].split('\nAP_Objective')[0].split('\n')) - 1
        ap_value = int(file_content_string.split('AP_Objective_Value: ')[1].split('\n')[0])
        ap_time = float(file_content_string.split('AP_Time: ')[1].split('\n')[0])

    return [lkh_value, lkh_gaps, lkh_time, ap_value, ap_gaps, ap_time]


def read_fasta_stats_file(file: str) -> Dict:
    with open(file, 'r') as f:
        file_content_string = f.read()
        actual_objective_value = int(file_content_string.split('Objective function value: ')[1].split('\n')[0])
        actual_gaps = int(file_content_string.split('Actual gaps: ')[1].split('\n')[0])
        no_of_reads = int(file_content_string.split('Number of reads: ')[1].split('\n')[0])
    return [no_of_reads, actual_objective_value, actual_gaps]


# def write_assembly_stats(assembly_stats_list: List[Assembly_Stats]) -> None:
#     with open('/home/andreas/GDrive/workspace/sparsedata/assembly_stats.csv', 'w') as f:
#         f_csv = csv.writer(f, delimiter=',')
#         f_csv.writerow(
#             ['File', 'LKHContigs', 'LKHValue', 'LKHTime', 'APContigs', 'APValue', 'APTime', 'ActualObjectiveValue'])
#         for elem in assembly_stats_list:
#             f_csv.writerow(elem)

def write_assembly_stats(statsdict: Dict) -> None:
    with open('/home/andreas/GDrive/workspace/sparsedata/assembly_stats.csv', 'w') as f:
        f_csv = csv.writer(f, delimiter=',')
        f_csv.writerow(
            ['Genome', 'Coverage', 'AvgLength', 'Reads', 'ActualValue', 'ActualGaps',
             'CalignLKHValue', 'CalignLKHGaps', 'CalignLKHTime',
             'CalignAPValue', 'CalignAPGaps', 'CalignAPTime',
             'CalignALKHValue', 'CalignALKHGaps', 'CalignALKHTime',
             'CalignAAPValue', 'CalignAAPGaps', 'CalignAAPTime',
             'CalignBLKHValue', 'CalignBLKHGaps', 'CalignBLKHTime',
             'CalignBAPValue', 'CalignBAPGaps', 'CalignBAPTime',
             ])
        for ref_name in [ref1_name, ref2_name, ref3_name]:
            for c in coverages:
                for length in average_length_list:
                    val = stats_dict[(ref_name, c, length)]
                    row = [ref_name, c, length]
                    row += val['Actual']
                    row += val['Calign']
                    row += val['Calign25']
                    row += val['Calign50']

                    f_csv.writerow(row)


def write_assembly_stats_tex(statsdict: Dict) -> None:
    with open('/home/andreas/GDrive/workspace/sparsedata/assembly_stats.tex', 'w') as f:
        for ref_name in [ref1_name, ref2_name, ref3_name]:
            if ref1_name == ref_name:
                dashline_active = ''
            else:
                dashline_active = '\\hdashline\n'
            f.write('{}\\bfseries {}\\\\\n'.format(dashline_active, ref_name))
            for c in coverages:
                f.write('$c = {}$\\\\\n'.format(c))
                for length in average_length_list:
                    val = stats_dict[(ref_name, c, length)]
                    row = [length]
                    row += [val['Actual'][0]]
                    row += ['']
                    row += val['Actual'][1:]
                    row += ['']
                    row += [*val['Calign'][0:2], '{0:.2f}'.format(val['Calign'][2]), *val['Calign'][3:5],
                            '{0:.2f}'.format(val['Calign'][5])]
                    row += ['']
                    row += [*val['Calign25'][0:2], '{0:.2f}'.format(val['Calign25'][2]), *val['Calign25'][3:5],
                            '{0:.2f}'.format(val['Calign25'][5])]
                    row += ['']
                    row += [*val['Calign50'][0:2], '{0:.2f}'.format(val['Calign50'][2]), *val['Calign50'][3:5],
                            '{0:.2f}'.format(val['Calign50'][5])]
                    f.write(' & '.join([str(x) for x in row]) + '\\\\\n')


def write_assembly_stats2(statsdict: Dict) -> None:
    with open('/home/andreas/GDrive/workspace/sparsedata/assembly_stats2.csv', 'w') as f:
        f_csv = csv.writer(f, delimiter=',')
        refs = [ref1_name, ref2_name]
        f_csv.writerow(range(len(refs) * 9))

        f_csv.writerow(
            [stats_dict[(ref_name, c, l)]['Actual'][0] for ref_name in refs for c in
             coverages for l in average_length_list])
        f_csv.writerow(
            [stats_dict[(ref_name, c, l)]['Actual'][1] for ref_name in refs for c in
             coverages for l
             in average_length_list])
        f_csv.writerow(
            [stats_dict[(ref_name, c, l)]['Actual'][2] for ref_name in refs for c in
             coverages for l
             in average_length_list])
        for foo in ['Calign', 'Calign25', 'Calign50']:
            for i in range(6):
                if i in [2, 5]:
                    f_csv.writerow(
                        ['{0:.2f}'.format(stats_dict[(ref_name, c, l)][foo][i]) for ref_name in refs for c in
                         coverages
                         for l in average_length_list])
                else:
                    f_csv.writerow(
                        [stats_dict[(ref_name, c, l)][foo][i] for ref_name in refs for c in
                         coverages
                         for l in average_length_list])


assembly_stats_list = []
stats_dict = {}
# for dir in sorted(glob.glob('/home/andreas/GDrive/workspace/sparsedata/ref[1,2,3]_c[5,20,40]*/')):
for ref_number in [1, 2, 3]:
    for coverage in coverages:
        for length in average_length_list:
            # file_sub_dir = dir.split('/')[-2]  # example ref1_c5_l100
            # ref_number = int(file_sub_dir.split('ref')[1].split('_')[0])
            ref_name = references[ref_number - 1]
            # coverage = int(file_sub_dir.split('_c')[1].split('_')[0])
            # length = int(file_sub_dir.split('_l')[1])
            dir = '/home/andreas/GDrive/workspace/sparsedata/ref{}_c{}_l{}/'.format(ref_number, coverage, length)
            stats_dict[(ref_name, coverage, length)] = {'Actual': read_fasta_stats_file(dir + 'fasta.stat'),
                                                        'Calign': read_assembly_file(dir + 'calign.assembly'),
                                                        'Calign25': read_assembly_file(
                                                            dir + 'calign_0_{}.assembly'.format(length // 4)),
                                                        'Calign50': read_assembly_file(
                                                            dir + 'calign_0_{}.assembly'.format(length // 2))}


            # dir = '{}-{}-{}'.format(references[ref_number - 1], coverage, length)
            # assembly_stats_list.append(
            #     Assembly_Stats(dir, len(lkh_contigs), lkh_value, lkh_time, len(ap_contigs), ap_value, ap_time,
            #                    actual_Objective_value))


def write_whole_stats() -> None:
    headers = ['CalignLKH', 'CalignAP', 'CalignALKH', 'CalignAAP', 'CalignBLKH',
               'CalignBAP']
    vals = {'CalignLKH': 0, 'CalignAP': 0, 'CalignALKH': 0, 'CalignAAP': 0, 'CalignBLKH': 0,
            'CalignBAP': 0}
    gaps = {'CalignLKH': 0, 'CalignAP': 0, 'CalignALKH': 0, 'CalignAAP': 0, 'CalignBLKH': 0,
            'CalignBAP': 0}
    both = {'CalignLKH': 0, 'CalignAP': 0, 'CalignALKH': 0, 'CalignAAP': 0, 'CalignBLKH': 0,
            'CalignBAP': 0}
    atspvsapval = {'CalignLKH': 0, 'CalignAP': 0, 'CalignALKH': 0, 'CalignAAP': 0, 'CalignBLKH': 0,
                   'CalignBAP': 0}
    atspvsap = {'CalignLKH': 0, 'CalignAP': 0, 'CalignALKH': 0, 'CalignAAP': 0, 'CalignBLKH': 0,
                'CalignBAP': 0}
    with open(DIR + 'assembly_stats.csv', 'r') as f:
        f_csv = csv.DictReader(f, delimiter=',')
        for row in f_csv:
            for elem in headers:
                if row['ActualValue'] == row[elem + 'Value']:
                    vals[elem] += 1
                if row['ActualGaps'] == row[elem + 'Gaps']:
                    gaps[elem] += 1
                if row['ActualValue'] == row[elem + 'Value'] and row['ActualGaps'] == row[elem + 'Gaps']:
                    both[elem] += 1
            if row['CalignLKHValue'] == row['CalignAPValue']:
                atspvsapval['CalignLKH'] += 1
                atspvsapval['CalignAP'] += 1
            if row['CalignALKHValue'] == row['CalignAAPValue']:
                atspvsapval['CalignALKH'] += 1
                atspvsapval['CalignAAP'] += 1
            if row['CalignBLKHValue'] == row['CalignBAPValue']:
                atspvsapval['CalignBLKH'] += 1
                atspvsapval['CalignBAP'] += 1
            if row['CalignLKHValue'] == row['CalignAPValue'] and row['CalignLKHGaps'] == row['CalignAPGaps']:
                atspvsap['CalignLKH'] += 1
                atspvsap['CalignAP'] += 1
            if row['CalignALKHValue'] == row['CalignAAPValue'] and row['CalignALKHGaps'] == row['CalignAAPGaps']:
                atspvsap['CalignALKH'] += 1
                atspvsap['CalignAAP'] += 1
            if row['CalignBLKHValue'] == row['CalignBAPValue'] and row['CalignBLKHGaps'] == row['CalignBAPGaps']:
                atspvsap['CalignBLKH'] += 1
                atspvsap['CalignBAP'] += 1
    with open(DIR + 'complete_stats.csv', 'w') as g:
        g_csv = csv.DictWriter(g, delimiter='&', fieldnames=headers)
        g_csv.writeheader()
        g_csv.writerow(vals)
        g_csv.writerow(gaps)
        g_csv.writerow(both)
        g_csv.writerow(atspvsapval)
        g_csv.writerow(atspvsap)


write_assembly_stats(stats_dict)
write_assembly_stats2(stats_dict)
write_assembly_stats_tex(stats_dict)
write_whole_stats()
