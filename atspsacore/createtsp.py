def write_full_atsp(filename, nr_of_reads, scores):
    tsplib_par_string = "PROBLEM_FILE={}\nOUTPUT_TOUR_FILE={}\nTOUR_FILE={}".format(filename + ".atsp",
                                                                                    filename + ".tour",
                                                                                    filename + ".tour")
    with open(filename + ".par", "w") as f:
        f.write(tsplib_par_string)

    with open(filename + ".atsp", "w") as f:
        # dimension +1 because of the special knot
        tsplib_string = "NAME: {}\nTYPE: ATSP \nCOMMENT: {}\nDIMENSION: {} \nEDGE_WEIGHT_TYPE: EXPLICIT \n" \
                        "EDGE_WEIGHT_FORMAT: FULL_MATRIX \nEDGE_WEIGHT_SECTION".format(filename, filename,
                                                                                       (nr_of_reads + 1))

        f.write(tsplib_string)

        # first line for the extra knot
        f.write('\n' + '\t'.join('0' for _ in range(nr_of_reads + 1)))

        for i in range(nr_of_reads):
            row = ['0', *[str(-scores[i][j]) for j in range(0, nr_of_reads)]]
            f.write('\n' + '\t'.join(row))

        f.write("\nEOF")
