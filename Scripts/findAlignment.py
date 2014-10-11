__author__ = 'tanya'

import sys
import re

def findSeqCoordinates(genome, begin, seq, dir):
    maxMove = 20
    begin_old = begin
    if dir == "left":
        while maxMove > 0 and begin >= 0:
            if genome[begin:begin+3] in seq:
                break
            else:
                begin -= 1
            maxMove -= 1
        if begin >= 0 and genome[begin:begin+3] in seq:
            return [begin, genome[begin: begin_old]]
        else:
            return [-1,""]
    else:
        while maxMove > 0 and begin < len(genome):
            if genome[begin-3:begin] in seq:
                break
            else:
                begin += 1
            maxMove -= 1
        if begin < len(genome) and genome[begin-3:begin] in seq:
            return [begin, genome[begin_old: begin]]
        else:
            return [-1,""]


def getSeq(genome, begin, end):
    ind = 0
    seq = genome[begin : end]
    if begin > end:
        seq = genome[end : begin]
        code = {"A":"T", "C":"G", "G":"C", "T":"A"}
        seq = "".join([code[x] for x in seq])
        seq = seq[::-1]

    return seq

def getCompliment(seq):
    code = {"A":"T", "C":"G", "G":"C", "T":"A"}
    seq = "".join([code[x] for x in seq])
    seq = seq[::-1]
    return seq

def readGenome(fileName):
    inputGenome = open(sys.argv[1],"r")
    inputGenome.readline()
    genome = ""
    for line in inputGenome:
        if line[-1] == "\n":
            genome += line[:-1]
        else:
            genome += line
    inputGenome.close()
    return genome

def getNucSeq(ln):
    if ln[-1] == "\n":
        return ln[:-1]
    else:
        return ln

def getAlignmentInfo(ln):
    lst = ln.split(" ")
    name = lst[0]
    length = int(lst[1])
    begin = int(lst[2][1:])
    end = int(lst[4][:-1])
    score = float(lst[5][:-1])
    return [name, length, begin, end, score]

def checkCodons(seq, codons):
    for i in range(0,len(seq)-3, 3):
        if seq[i:i+3] in codons:
            print(seq[i:i+3])
            print(i)
            return i
    return -1

def testGoodSeq(seq, genome, beginOld, endOld):
    if beginOld < endOld:
        startCodons = {"ATG"}
        stopCodons = {"TAG", "TAA", "TGA"}
        [begin, prefix] = findSeqCoordinates(genome, beginOld, startCodons, "left")
        if begin == -1:
            print("FAILED: start codon not found")
            return [seq, -1, -1]
        [end, suffix] = findSeqCoordinates(genome, endOld, stopCodons, "right")
        if end == -1:
            print("FAILED: stop codon not found")
            return [seq, -1, -1]
        seq = prefix + seq +suffix[3:]
    else:
        startCodons = {"CAT"}
        stopCodons = {"CTA", "TTA", "TCA"}
        [begin, prefix] = findSeqCoordinates(genome, beginOld, startCodons, "right")
        if begin == -1:
            print("FAILED: start codon not found")
            return [seq, -1, -1]
        [end, suffix] = findSeqCoordinates(genome, endOld, stopCodons, "left")
        if end == -1:
            print("FAILED: stop codon not found")
            return [seq, -1, -1]

        seq = getCompliment(prefix) +seq + getCompliment(suffix)[3:]


    [div, mod] = divmod(len(seq), 3)
    if mod != 0:
       print("FAILED: sequence len % 3 not 0")
       print("Real: ("+str(begin) +" - " + str(end) +")")
       return [seq, -1, -1]

    stopCodons = {"TAG", "TAA", "TGA"}
    if checkCodons(seq[:-3], stopCodons) != -1:
        print("FAILED: stop codon in sequence")
        print("Real: ("+str(begin) +" - " + str(end) +")")
        return [seq, -1, -1]

    print("OK: ("+str(begin) +" - " + str(end) +")")
    return [seq, begin, end]


genome = readGenome(sys.argv[1])

inputExRes = open(sys.argv[2], "r")

resFile = open(sys.argv[3], "w")
ln = inputExRes.readline()
while ln != "":
    while not ln.startswith("p") and ln != "":
        ln = inputExRes.readline()
    if ln == "":
        break
    print(ln[:-1])
    resFile.write(ln)
    name, length, begin, end, score = getAlignmentInfo(ln)
    nucs = {"A","C","G","T"}
    inputExRes.readline()
    ln = inputExRes.readline()
    seq = getNucSeq(ln)[2:]
    ln = inputExRes.readline()
    while ln[0] in nucs:
        seq += getNucSeq(ln)
        ln = inputExRes.readline()
    [seq, begin, end] = testGoodSeq(seq, genome, begin, end)
    resFile.write("Real: ("+str(begin) +" - " + str(end) +")" +"\n")
    resFile.write(seq+"\n")

inputExRes.close()
resFile.close()

