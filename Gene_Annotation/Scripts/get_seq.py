__author__ = 'tanya'

import sys

input = open(sys.argv[1],"r")
begin = int(sys.argv[2])
end = int(sys.argv[3])


input.readline()
ind = 0
seq = ""
rc =False
code = {"A":"A", "C":"C", "G":"G", "T":"T"}
if begin > end:
    rc = True
    code = {"A":"T", "C":"G", "G":"C", "T":"A"}
    h = begin
    begin = end
    end = h

for line in input.readlines():
    ln = line[:-1]
    for it in ln:
        ind += 1
        if ind >= begin and ind <= end:
            seq += code[it]
    if ind > end:
        break

input.close()

if rc:
    seq = seq[::-1]

print(seq)