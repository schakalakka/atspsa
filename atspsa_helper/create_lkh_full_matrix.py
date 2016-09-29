from  atspsacore import createtsp, read_score_file


def create_full_atsp(file):
    current_file = file.split('.score')[0] + '_full'
    print(current_file)
    nr_of_reads, scores = read_score_file.read_score_file(file)
    scores[scores < 30] = -99999

    createtsp.write_full_atsp(current_file, nr_of_reads, scores)


# for file in glob.glob(DIR+'*/*.score'):
#     create_full_atsp(file)
create_full_atsp('/home/andreas/GDrive/workspace/sparsedata/ref1shuffled_c5_l700/calign.score')
