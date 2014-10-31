__author__ = 'tanya'

import sys

def getIdealProteins(fileName):
    fin = open(fileName, "r")
    lines = fin.readlines()
    res = []
    for line in lines:
        if line.startswith(">"):
            targetName = line.split(" ")[0][1:]
            res.append("")
        else:
            if line.endswith("\n"):
                res[-1] += line[:-1]
            else:
                res[-1] += line
    fin.close()
    return res

def getAugustusProteins(fileName):
    fin = open(fileName, "r")
    lines = fin.readlines()
    res = []
    nowProtein = False
    for line in lines:
        if line.startswith("# protein sequence = ["):
            nowProtein = True
            res.append("")
            if line.endswith("\n"):
                res[-1] += line[len("# protein sequence = ["):-1]
            else:
                res[-1] += line[len("# protein sequence = ["):]
        elif line.startswith("# end gene"):
            nowProtein = False
        elif nowProtein:
            if line.endswith("]\n"):
                res[-1] += line[2:-2]
            elif line.endswith("\n"):
                res[-1] += line[2:-1]
            else:
                res[-1] += line[2:]
    fin.close()
    return res


idealProteinsFile = sys.argv[2]
augustusProteinsFile = sys.argv[1]

idealProteins = getIdealProteins(idealProteinsFile)
augustusProteins = getAugustusProteins(augustusProteinsFile)

acc = 0
for it in augustusProteins:
    for it2 in idealProteins:
        if it == it2:
            acc += 1
            break

print("Total accuracy: " + str(acc*1.0/len(augustusProteins)))
print("Total proteins: " + str(acc))
