def prepare_lkh(filename, nr_of_reads, scores):
    print("Preparing for LKH")
    tsplib_par_string = "PROBLEM_FILE={}\nOUTPUT_TOUR_FILE={}".format(filename + ".atsp", filename + ".tour")
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
            f.write('\n' + '\t'.join([str(-scores.get((i, j), 0)) for j in range(-1, nr_of_reads)]))

        f.write("\nEOF")
    print("Preparing for LKH finished.")
