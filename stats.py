__author__ = 'andreas'

from collections import defaultdict

import atspsa
import db2score


def make_stat(filename):
    # reads = atspsa.parse_fasta(filename)
    overlaps = db2score.get_scores_without_orientation(db2score.read_ovl_file(filename))
    tour = atspsa.parse_tour(filename)

    # make overlap_help dict
    ovl_sets = defaultdict(set)
    ovl_lists = defaultdict(list)
    [ovl_sets[i].add(overlaps[(i, j)]) for (i, j) in overlaps.keys()]
    for i in ovl_sets.keys():
        ovl_lists[i] = sorted(ovl_sets[i], key=int, reverse=True)

    # count how often was the best (second best...) edge used
    nr_best_edge_used = [0] * 20
    biggerthan = 0
    for i, elem in enumerate(tour):
        if i == len(tour) - 1:
            break
        next = tour[i + 1]
        try:
            if overlaps[(elem, next)] in ovl_lists[elem]:
                n = ovl_lists[elem].index(overlaps[(elem, next)])
                nr_best_edge_used[n + 1] += 1
        except KeyError:  # used expensive edge
            nr_best_edge_used[0] += 1
        except IndexError:
            print("n=", n)
            biggerthan += 1
    return (nr_best_edge_used, overlaps)
