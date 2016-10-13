def write_full_atsp(filename, nr_of_reads, scores, circular=False):
    if circular:
        add_extra_non_circular_node = 0
    else:  # adds an extra node to transform a path to a tour
        add_extra_non_circular_node = 1
    tsplib_par_string = "PROBLEM_FILE={}\nOUTPUT_TOUR_FILE={}\nTOUR_FILE={}".format(filename + ".atsp",
                                                                                    filename + ".tour",
                                                                                    filename + ".tour")
    with open(filename + ".par", "w") as f:
        f.write(tsplib_par_string)

    with open(filename + ".atsp", "w") as f:
        # dimension +1 because of the special knot
        tsplib_string = "NAME: {}\nTYPE: ATSP \nCOMMENT: {}\nDIMENSION: {} \nEDGE_WEIGHT_TYPE: EXPLICIT \n" \
                        "EDGE_WEIGHT_FORMAT: FULL_MATRIX \nEDGE_WEIGHT_SECTION".format(filename, filename,
                                                                                       (
                                                                                       nr_of_reads + add_extra_non_circular_node))

        f.write(tsplib_string)

        # first line for the extra knot if it is not a circular genome
        if not circular:
            f.write('\n' + '\t'.join('0' for _ in range(nr_of_reads + add_extra_non_circular_node)))

        for i in range(nr_of_reads):
            if not circular:
                row = ['0', *[str(-scores[i][j]) for j in range(0, nr_of_reads)]]
            else:
                row = [str(-scores[i][j]) for j in range(0, nr_of_reads)]
            f.write('\n' + '\t'.join(row))

        f.write("\nEOF")
