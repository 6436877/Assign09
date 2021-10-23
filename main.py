# Input
import re

from seqbio.calculation.SeqCal import countBasesDict, gcContent
from seqbio.pattern.SeqPattern import enzTargetsScan, cpgSearch, reverseSeq, complementSeq, reverseComplementSeq, dna2protein, dna2rna

def test():
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


def main():
    from argparse import ArgumentParser
    print("HEY")

    parser = ArgumentParser(prog="myseq", description="work with seq")
    args = parser.parse_args()
    if args.command == "gcContent":
        if args.seq == None:
            print("Error")
            exist(parser.parse_args(['gcContent','-h']))
    elif args.command == "enzTargetScan":
        if args.seq == None or args.enz == None:
            print("Error:invalid input")
            exist(parser.parse_args(["enzTargetScan","-h"]))
        print(args.seq, gcContent(args.seq),enzTargetsScan(args.seq, args.enz))

    subparsers = parser.add_subparsers(
        title="commands", description="Pls choose command below:",
        dest='command'
    )
    subparsers.required = True
    cgc_command = subparsers.add_parser('gcContent', help='Calculate GC content')
    cgc_command.add_argument('-s','--seq', type=str, default=None,help="Provide seq")
    cgc_command.add_argument('-e','--enz', type=str, default=None,help="Enz name")
    
    renz_command = subparsers.add_parser('enzTargetsScan', help='Find res enz')
    renz_command.add_argument('-s','--seq', type=str, default=None,help="Provide seq")
    renz_command.add_argument('-e','--enz', type=str, default=None,help="Enz name")

    
    return parser
if __name__ == "__main__":
    test()

