def prepare_lkh(LKH_PAR, LKH_LIB, LKH_OUT, nr_of_reads, scores):
    tsplib_par_string = "PROBLEM_FILE={}\nOUTPUT_TOUR_FILE={}".format(LKH_LIB, LKH_OUT)
    with open(LKH_PAR, "w") as f:
        f.write(tsplib_par_string)

    with open(LKH_LIB, "w") as f:
        # dimension +1 because of the special knot
        tsplib_string = "NAME: {}\nTYPE: ATSP \nCOMMENT: {}\nDIMENSION: {} \nEDGE_WEIGHT_TYPE: EXPLICIT \n" \
                        "EDGE_WEIGHT_FORMAT: FULL_MATRIX \nEDGE_WEIGHT_SECTION".format(LKH_LIB, LKH_LIB,
                                                                                       (nr_of_reads + 1))

        f.write(tsplib_string)

        # first line for the extra knot
        f.write('\n' + '\t'.join('0' for _ in range(nr_of_reads + 1)))

        for i in range(nr_of_reads):
            f.write('\n' + '\t'.join([str(-scores.get((i, j), 0)) for j in range(-1, nr_of_reads)]))

        f.write("\nEOF")
