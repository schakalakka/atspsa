import sys

sys.path.insert(0, '/home/andreas/workspace/align/')

import align


def read_align_file(filename):
    return 0


def write_align_file(filename, scores):
    with open(filename + '.align', 'w') as f:
        for key, elem in scores.items():
            f.write('{}\t{}\t{}\n'.format(key[0], key[1], elem))


def compute_scores(reads):
    # return {(i,j): align.Align(x, y, init=True).get_score() for i,x in enumerate(reads) for j,y in enumerate(reads)}
    scores = {}
    for i, x in enumerate(reads):
        for j, y in enumerate(reads):
            if i < j:
                print(i, j)
                scores[(i, j)] = align.Align(x, y, init=True).get_score()
    return scores
