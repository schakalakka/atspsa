####################
# reads the alignmenttime_stats.csv
#   containing the time for every alignment method (Calign, calign25, Calign50)
#
# creates alignmentcomparison.csv
#   result is a table
#       each entry is the speed up (i.e. x times faster) of the alignment method compared to normal SeqAn
####################

import csv

from config import *

stats = []
with open(DIR + 'alignmenttime.csv', 'r') as f:
    f_csv = csv.reader(f, delimiter='\t')
    headers = next(f_csv)
    # Row = namedtuple('Row', headers)
    rowlist = []
    for line in f_csv:
        rowlist.append([*line])

    for i, aligner in enumerate(headers[5:]):
        current_stat = {}
        current_stat['Aligner'] = aligner
        current_stat['Number'] = i
        for length in average_length_list:
            current_rows = [x for x in rowlist if x[3] == str(length) and x[4] != '-1' and x[i + 5] != '-1']
            numerator = 0
            for row in current_rows:
                numerator += float(row[4]) / float(row[5 + i])
            current_stat[str(length)] = round(numerator / len(current_rows))
        stats.append(current_stat)

with open(DIR + 'alignmentcomparison.csv', 'w') as f:
    headers = ['Number', 'Aligner', *[str(x) for x in average_length_list]]
    f_csv = csv.DictWriter(f, headers, delimiter='\t')
    f_csv.writeheader()
    f_csv.writerows(stats)

print(stats)
