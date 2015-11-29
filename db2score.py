def read_ovl_file(filename):
    print("Reading OVL file.")
    scores = {}
    with open(filename + ".ovl", "r") as f:
        for line in f:
            if line.startswith("P"):
                tail = int(line.split(" ")[1])
                head = int(line.split(" ")[2])
                orientation = line.split(" ")[3].split("\n")[0]
            elif line.startswith("C"):
                len1 = float(line.split(" ")[2]) - float(line.split(" ")[1])
                len2 = float(line.split(" ")[4]) - float(line.split(" ")[3])
            elif line.startswith("D"):
                diff = float(line.split(" ")[1])
                if (tail, head) in scores:
                    if scores[(tail, head)] > int((len1 + len2) / 2 - diff):
                        break
                scores[(tail, head)] = (int((len1 + len2) / 2 - diff), orientation)
    print("Reading OVL file finished")
    return scores


def get_scores_without_orientation(scores):
    return {key: scores[key][0] for key in scores.keys()}

