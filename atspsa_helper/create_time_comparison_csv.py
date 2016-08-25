import csv
from collections import namedtuple

from config import *

stats = []
with open(DIR + 'alignmenttime_stats.csv', 'r') as f:
    f_csv = csv.reader(f, delimiter='\t')
    headers = next(f_csv)
    Row = namedtuple('Row', headers)
    rowlist = []
    for line in f_csv:
        rowlist.append(Row(*line))

    for i, aligner in enumerate(headers[5:]):
        current_stat = {}
        current_stat['Aligner'] = aligner
        current_stat['Number'] = i
        for length in average_length_list:
            current_rows = [x for x in rowlist if x.Length == str(length) and x.SeqAn != '-1' and x[i + 5] != '-1']
            numerator = 0
            for row in current_rows:
                numerator += int(row.SeqAn) / int(row[5 + i])
            current_stat[str(length)] = int(numerator / len(current_rows))
        stats.append(current_stat)

with open(DIR + 'alignmentcomparison.csv', 'w') as f:
    headers = ['Number', 'Aligner', *[str(x) for x in average_length_list]]
    f_csv = csv.DictWriter(f, headers, delimiter='\t')
    f_csv.writeheader()
    f_csv.writerows(stats)

print(stats)
