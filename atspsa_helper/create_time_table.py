import os
import sys

sys.path.insert(0, "/home/andreas/GDrive/workspace/randomdnadatagenerator/")

DIR = "/home/andreas/GDrive/workspace/sparsedata/"

output = '\\begin{tabular}{|c|c|c|c|c|c|}\\hline\nRef & Coverage & Length & Reads & SeqAn & Calign \\\\\n\\hline\n'
for ref in [1, 2, 3]:
    for coverage in coverages:
        for length in average_length_list:
            if os.path.exists(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}.act'.format(ref, coverage, length)):
                with open(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}.act'.format(ref, coverage, length),
                          'r') as f:
                    seq_align_time = float(f.read())
            else:
                seq_align_time = -1
            if os.path.exists(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}.cact'.format(ref, coverage, length)):
                with open(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}.cact'.format(ref, coverage, length),
                          'r') as f:
                    c_align_time = float(f.read())
                with open(DIR + 'ref{0}shuffled_c{1}_l{2}/shuffled_c{1}_l{2}.calignscore'.format(ref, coverage, length),
                          'r') as f:
                    nr_of_reads = int(f.readline())
            else:
                c_align_time = -1
                nr_of_reads = -1
            output = output + '{0} & {1} & {2} & {3} & {4:.1f} & {5:.1f} \\\\\n'.format(ref, coverage, length,
                                                                                        nr_of_reads, seq_align_time,
                                                                                        c_align_time)
    output = output + '\\hline\n'
output = output + '\\end{tabular}'
with open(DIR + 'alignmenttime' + '_stats.tex', 'w') as f:
    f.write(output)
print(output)
