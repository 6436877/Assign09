import re

# SeqCal module
def gcContent(seq):
    # G+C/(A+T+G+C)
    return (countBase(seq, 'G') + countBase(seq, 'C'))/len(seq)

def atContent(seq):
    # A+T/(A+T+G+C)
    return (countBase(seq, 'A') + countBase(seq, 'T'))/len(seq)

def countBase(seq, base):
    return seq.count(base.upper())

def countBasesDict(seq):
    basesM = {}
    for base in seq:
        basesM[base] = basesM.get(base, 0)+1
    return basesM

# SeqPattern module
def cpgSearch(seq):
    cpgs = []
    for m in re.finditer(r'CG', seq, re.I):
        cpgs.append((m.group(), m.start(), m.end()))
    return cpgs

def enzTargetsScan(seq, enz):
    resEnzyme = dict(EcoRI='GAATTC', BamHI='GGATCC', 
                 HindIII='AAGCTT',AccB2I='[AG]GCGC[CT]',
                 AasI='GAC[ATCG][ATCG][ATCG][ATCG][ATCG][ATCG]GTC',
                 AceI='GC[AT]GC')
    
    out = []
    if enz in resEnzyme:
        for m in re.finditer(resEnzyme.get(enz,),seq):
            out.append((m.group(0),m.start(),m.end()))
    return out

def reverseSeq(seq):
    return seq[::-1]

def complementSeq(seq):
    compl = {"A": "T", "T": "A",
             "G": "C", "C": "G"}
    complementary = "".join([ compl[base] for base in seq ])
    return complementary

def reverseComplementSeq(seq):
    revComp = complementSeq(reverseSeq(seq))
    return revComp

def dna2rna(seq):
    return seq.replace("T","U")

def dna2protein(seq):
    DNA_Codons = loadCodons()
    protein = ""
    for i in range(0,len(seq),3):
        dna = seq[i:i+3]
        protein += DNA_Codons.get(dna, "")
    return protein

def loadCodons():
    DNA_Codons = {
        # 'M' - START, '_' - STOP
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "TGT": "C", "TGC": "C",
        "GAT": "D", "GAC": "D",
        "GAA": "E", "GAG": "E",
        "TTT": "F", "TTC": "F",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        "CAT": "H", "CAC": "H",
        "ATA": "I", "ATT": "I", "ATC": "I",
        "AAA": "K", "AAG": "K",
        "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "ATG": "M",
        "AAT": "N", "AAC": "N",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
        "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
        "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "TGG": "W",
        "TAT": "Y", "TAC": "Y",
        "TAA": "_", "TAG": "_", "TGA": "_"
    }
    return DNA_Codons

# Input
seq = 'ATGGGccGTAGAATTCTTGCaaGCCCGT'
seq = seq.upper()
print("Transcription: ", dna2rna(seq))
print("Transcription-revcomp: ", dna2rna(reverseComplementSeq(seq)))
print("Translation: ", dna2protein(seq))
print("Translation-revcomp: ", dna2protein(reverseComplementSeq(seq)))
print("GC Content:", gcContent(seq))
print("Count Bases: ", countBasesDict(seq))
print("Count Bases-revcomp: ", countBasesDict(reverseComplementSeq(seq)))
print("Search EcoRI: ", enzTargetsScan(seq, 'EcoRI'))
print("Search EcoRI-revcomp: ", enzTargetsScan(reverseComplementSeq(seq), 'EcoRI'))
