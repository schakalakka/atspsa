def read_ovl_file(file):
    scores = {}
    with open(file, "r") as f:
        for line in f:
            if line.startswith("P"):
                tail = int(line.split(" ")[1])
                head = int(line.split(" ")[2])
                orientation = line.split(" ")[3]
            elif line.startswith("C"):
                len1 = float(line.split(" ")[2]) - float(line.split(" ")[1])
                len2 = float(line.split(" ")[4]) - float(line.split(" ")[3])
            elif line.startswith("D"):
                diff = float(line.split(" ")[1])
                if (tail, head) in scores:
                    if scores[(tail, head)] > diff / (len1 + len2):
                        break
                scores[(tail, head)] = (int(100 * diff / (len1 + len2)), orientation)
    return scores


def get_scores_without_orientation(scores):
    return {key: scores[key][0] for key in scores.keys()}


read_ovl_file("/run/media/andreas/INTENSO/fastas/ecoli_1000/ecoli_1000.db.ovl")
