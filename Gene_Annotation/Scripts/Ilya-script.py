
def DNAtoProtein (DNAseq, frame, letter):
    """Converts DNA sequence into protein sequence.
    frame: +1, +2, +3, -1, -2 or -3
    letter: 1 or 3"""
    assert frame >= - 3 and frame <= 3 and frame != 0, 'frame must be +1, +2, +3, -1, -2 or -3'
    assert letter == 1 or letter == 3, 'letter must be 1 or 3'
    resultDNA = DNAframe(DNAseq, frame) 
    resultProt = readCode(resultDNA)
    if letter == 1:
        return resultProt
    else:
        return protThree(resultProt)

def revSeq(seq):
    """Returns any sequence in the reverse order"""
    rev = ''
    for n in seq:
        rev = n + rev
    return rev

def DNAcompl(DNAseq):
    """Returns the compltement strand of DNA"""
    #print 'DNAcompl called with', DNAseq
    compl = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
             'a': 't', 't': 'a', 'g': 'c', 'c': 'g', ' ': ' '}
    complStrand = ''
    for nuc in DNAseq:
        complStrand += compl[nuc]
    return complStrand
    
def DNAframe(DNAseq, startPoint):
    """Returns DNA sequence in six different reading frames"""
    assert startPoint != 0, 'startPoint must be +1, +2, +3, -1, -2 or -3'
    #print 'DNAframe called with', DNAseq, startPoint
    if startPoint > 0:
        return DNAseq[startPoint-1:]
    else:
        newStrand = DNAcompl(revSeq(DNAseq))
        startPoint *= (-1)
        return DNAframe(newStrand, startPoint)
    
def readCode(DNAseq):
    """Reads DNA code and returns a protein sequence written in the one-letter way"""
    DNAcode = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
               'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
               'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
               'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
               'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
               'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
               'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
               'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
               'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
               'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
               'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
               'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
               'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
               'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
               'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
               'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
    protSeq = ''
    triplNum = len(DNAseq)/3
    for i in range(triplNum):
        tripl = DNAseq[i*3: i*3+3]
        #print tripl
        aa = DNAcode[tripl]
        #print aa,
        protSeq += aa
    return protSeq

def protThree(protSeq):
    """Converts the one-letter way of writing protein sequence into the three-latter way"""
    aaNames = {'F': 'Phe', 'L': 'Leu', 'S': 'Ser', 'Y': 'Tyr', '*': 'STP',
               'W': 'Trp', 'P': 'Pro', 'H': 'His', 'Q': 'Gln', 'R': 'Arg',
               'I': 'Ile', 'M': 'Met', 'T': 'Thr', 'N': 'Asn', 'K': 'Lys',
               'V': 'Val', 'A': 'Ala', 'D': 'Asp', 'E': 'Glu', 'G': 'Gly',
               'C': 'Cys'}
    newSeq = ''
    for i in protSeq:
        newSeq += aaNames[i]
    return newSeq

def showDNA(name, seq, length, numGroup):
    """Prints DNA sequnces in GenBank-like manner"""
    print '>', name
    print
    finalString = ''
    stringLength = length*numGroup                      # calculates the length of 1 string
    nucNum = 1                                        
    longestNum = len(str(len(seq)))                     # calculates the length of the biggest nucleotide number
    while len(seq) > 0:
        DNAseq = DNAstring(seq[:stringLength], length)  # slices the original sequence into pieces of particular length and sends a message to DNAstring
        seq = seq[stringLength:]                        # removes the pieces from the original sequence
        gap = longestNum - len(str(nucNum))+1           # calculates the length of gaps between the nucleotide number and the DNA sequence of a string
        finalString = str(nucNum) + ' '*gap + DNAseq    # assembling strings
        print finalString
        nucNum += stringLength                          # calculates nucleotide number

def DNAstring(seq, length):
    """Slices DNA sequences into pieces of particular length devided by gaps"""
    newString = ''
    while len(seq) > 0:
        newString = newString + seq[:length] + ' '
        seq = seq[length:]
    return newString

def show_dsDNA(name, seq, length, numGroup):
    print '>', name
    print
    senseString = ''
    antisenseString = ''
    stringLength = length*numGroup                      # calculates the length of 1 string
    nucNum = 1                                        
    longestNum = len(str(len(seq)))                     # calculates the length of the biggest nucleotide number
    while len(seq) > 0:
        DNAseq = DNAstring(seq[:stringLength], length)  # slices the original sequence into pieces of particular length and sends a message to DNAstring
        compl = DNAcompl(DNAseq)                        # generates the second string
        gap = longestNum - len(str(nucNum))+1           # calculates the length of gaps between the nucleotide number and the DNA sequence of a string
        senseString = str(nucNum) + ' '*gap + DNAseq    # assembling the first strings
        antisenseString = ' '*len(str(nucNum)) + ' '*gap + compl # assembling the second string
        print senseString
        print antisenseString
        print
        nucNum += stringLength                          # calculates nucleotide number
        seq = seq[stringLength:]                        # removes the pieces from the original sequence

def findSeq(find_seq, in_seq):
    """Finds any sequence in another sequence"""
    seqLen = len(find_seq)
    found = False
    positions = []
    itemNum = 0
    for i in in_seq:
        if find_seq == in_seq[itemNum : itemNum+seqLen]:
            found = True
            positions.append(itemNum)
        itemNum += 1
    if found == True:
        return positions
    else: return None
                
##def locatePrimers(first, second, seq):
##    Forward = findSeq(first, seq)
##    R = revSeq(DNAcompl(second))
##    Reverse = findSeq(R, seq)
##    if Forward == None:
##        print 'Primer', first, 'hasn\'t been located'
##    else:
##        for i in Forward: print 'Primer', first, 'has been located at position', i        
##    if Reverse == None:
##        print 'Primer', second, 'hasn\'t been located'
##    else:
##        for i in Reverse:
##            print 'Primer', second, 'has been located at position', i
    
    
gene = 'CGCAGTGACGGTGAAGGTGGAACTGAGCTGCGCACCCGGGGATCTCGATGCCGTCCTCATCCTGCAGGGTCCCCCCTACGTGTCCTGGCGTGACGGTGAAGGTGGAACTGAGCTGCGCAGCGCACGCAGTGACGGTGAAGGTGGAACTGAGCTGCGCACCCGGGGATCTCGATGCCGTCCTCATCCTGCAGGGTCCCCCCTACGTGTCCTGGCGTGACGGTGAAGGTGGAACTGAGCTGCGCAGCGCATTCACCCGACCCAGAAATTAGCTATCTGATAAGCAAAATAACATTTTCAGAATCTTATCGAAGCTGTTAATTACGATACGAATCCGGTCGTCGTTCACTGCTCTGCTGGAGTTGGTCGTTCTGGAACTATTGTGGGAATTTCTTTGATTATGGATAAGATGATTCAGGGAGTGAGTTTGAAATTTTTCATGTCTTTCGTTTCAATTTGCTATATTTTTCAGATTAATTGCAAAGACATGAAAAAGTTAGTCGAGGAAATCCGCAATCAGCGTCACTATGCCATTCAGACGGAGGCGGTATGCAAATTTCCAGTTATTATTATAAAATAGAAGTCAATTTCAATTTCAGCAATACATGTACATTCACCGTGTTCTCCTCGAATACTTCTTGGAATTGCACAAAGAGACGTACGAAGGATTGTTGCTGACAAAGAATTACGAGGAAAAGTATGCAAAGTGGCTGGACGATTATAACACCTACGCAGAAGCCAACGATCGTACTCAAGCAGCTAAGGCGGCCCGCGCAGCCTTGT'

DNAtoProtein(gene, 1, 3)



