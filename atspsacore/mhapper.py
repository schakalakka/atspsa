import os


def run_mhap(MHAP_JAR, FASTA_FILE, MHAP_OUT):
    # run mhap
    os.system("java -server -Xmx32g -jar {} -s {} > {}".format(MHAP_JAR, FASTA_FILE, MHAP_OUT))


def parse_mhap(MHAP_OUT):
    dd = {}
    # read mhap output
    with open(MHAP_OUT) as f:
        for line in f:
            if line.startswith("#"):
                break
            else:
                line_arr = line.split(" ")
                i = int(line_arr[0]) - 1
                j = int(line_arr[1]) - 1
                dd[(i, j)] = int(float(line_arr[2]))  # mhap counts from 1 to nr_of_reads
    return dd
