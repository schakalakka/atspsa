import os
from collections import namedtuple

from atspsacore import read_score_file, parser

Node = namedtuple('Node', ['id', 'sequence', 'start_index', 'end_index', 'previous_node', 'next_node', 'length'])

file = '/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c5_l700/calign.score'
tour_file = '/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c5_l700/calign_full.tour'
fasta_file = '/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c5_l700/shuffled_c5_l700.fasta'

nr_of_reads, scores = read_score_file.read_score_file(file)

tour = parser.parse_tour(tour_file)

reads = parser.parse_fasta_with_id(fasta_file)


def align(seq1, seq2):
    foo = os.popen(
        '/home/andreas/GDrive/workspace/seqan-project/bin/seqconsensus --align {} {}'.format(seq1, seq2))
    # print(foo.read())
    return foo


def consensus():
    foo = ' '.join(
        str(x) for x in tour)
    foo = os.popen(
        '/home/andreas/GDrive/workspace/seqan-project/bin/seqconsensus --consensus {} '.format(fasta_file) + ' '.join(
            str(x) for x in tour))
    # print(foo.read())
    return foo


nodes = []
for i, read in enumerate(reads):
    start_index, end_index = read[0].split('/')[-1].split('_')
    if tour.index(i) < len(tour) - 1:
        next_node = tour[tour.index(i) + 1]
    else:
        next_node = -1
    if tour.index(i) > 0:
        previous_node = tour[tour.index(i) - 1]
    else:
        previous_node = -1
    nodes.append(
        Node(i, read[1], int(start_index), int(end_index), previous_node, next_node, int(end_index) - int(start_index)))

nodes = sorted(nodes, key=lambda x: tour.index(x.id))


def check_for_gaps():
    last_seq = nodes[0].sequence
    result = [nodes[0].sequence]
    for i, node in enumerate(nodes):
        if i < len(nodes) - 1:
            out = align(last_seq, nodes[i + 1].sequence)
            output = out.read()
            foo = [x for x in output.split('\n')]
            values = [int(x.split(': ')[1]) for x in output.split('\n')]
            if values[1] != values[0] + values[2] + values[3]:
                print('hi')
            if values[5] != values[4] + values[6] + values[7]:
                print('ho')

            if values[2] == 0 and values[3] > 0 and values[6] > 0 and values[7] == 0:  # normal overlap:
                # leading gaps for first sequence exists
                # no trailing gaps for first sequence
                # no leading gaps for second sequence
                # trailing gaps for second sequence exists
                start_of_added_sequence = len(last_seq) - values[6]
                # xxxxxxyyyyyyyy------------
                # ------yyyyyyyyzzzzzzzzzzzz
                #               ^- start of next sequence
                # leading gaps of sequence 2 = values[6]
                # the zzzzzzzzzz part is added to the list
                result.append(nodes[i + 1].sequence[start_of_added_sequence:])
                last_seq = nodes[i + 1].sequence
            elif values[6] > 0 and values[7] > 0:  #
                pass
            elif values[2] > 0 and values[3] > 0:  #
                print('hä?')
            else:
                print('hmm')
    return result


bar = check_for_gaps()
print(''.join(bar))


# print(sorted(nodes, key=lambda x: tour.index(x.id)))
# for i, node in enumerate(nodes):
#     if i < len(nodes)-1:
#         align(nodes, i)
# print(nodes[0].sequence)

# consensus()
