from Bio import SeqIO

def parse() -> str:
    sequence = ""
    record = SeqIO.read('source/sequence.fasta','fasta')
    sequence = record.seq

    return sequence
