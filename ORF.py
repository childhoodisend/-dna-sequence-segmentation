def isstart(codon):
    """Check if triplet is start codon"""

    startcodons = ["ATG", "GTG", "TTG"]

    if codon in startcodons:
        return True
    else:
        return False


def isstop(codon):
    """Check if triplet is stop codon"""

    stopcodons = ["TAA", "TAG", "TGA"]

    if codon in stopcodons:
        return True
    else:
        return False


def revcompdna(DNAseq):
    """Get reverse complement DNA"""

    # complements = {'A': 'T', 'T': 'A', 'U': 'A', 'C': 'G', 'G': 'C', 'N': 'N', 'R': 'Y', 'Y': 'R', 'M': 'K', 'K': 'M', 'W': 'W', 'S': 'S', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    if not isinstance(DNAseq, str):
        assert False, "Wrong DNA sequence format (type)"

    try:
        return "".join([complements[base] for base in reversed(DNAseq.upper())])
    except KeyError:
        assert False, "Wrong DNA sequence format (code)"


def find_codons(DNAseq, frame=1):
    """Get all tripletes based on chosen frame"""

    frames = [1, 2, 3, -1, -2, -3]

    collector = []

    if isinstance(DNAseq, str) and frame in frames:
        if (frame < 0):
            DNAseq = revcompdna(DNAseq)
            frame = abs(frame)
        else:
            DNAseq = (DNAseq)
        triplet = ""
        for nuc in DNAseq[frame - 1:]:
            triplet = triplet + nuc
            if len(triplet) == 3:
                collector.append(triplet)
                triplet = ""
        return collector
    else:
        assert False, "Wrong input format"

class ORF(object):
    """Class object for Open Reading Frame"""

    def __init__(self, seq, strain, frame, start, stop):
        self.seq = seq
        self.length = len(seq)
        self.strain = strain
        self.frame = frame
        self.start = start
        self.startcod = seq[:3]
        self.stop = stop
        self.stopcod = seq[-3:]
        self.annotation = ""

    def report(self):
        print("Последовательность ORF:", self.seq)
        print("Длина последовательности:", self.length, "п.н.")
        print("Рамка считывания:", self.frame)
        print("Стартовая позиция:", self.start)
        print("Стоповая позиция:", self.stop)


def ORF_finder_pro(DNAseq, minlen=10):
    frames = [1, 2, 3, -1, -2, -3]
    # started = []
    ORFs = []
    final_result = []

    for frame in frames:
        codons = find_codons(DNAseq, frame=frame)

        for i, codon in enumerate(codons):

            if isstart(codon):
                for orf in ORFs:
                    orf[0].append(codon)
                if frame > 0:
                    ORFs.append([[codon], (i + 1) * 3 - 2 + frame - 1])
                else:
                    ORFs.append([[codon], len(DNAseq) - (i + 1) * 3 + 3 + frame + 1])
            elif isstop(codon):
                for orf in ORFs:
                    orf[0].append(codon)
                    result = "".join(orf[0])
                    if len(result) >= minlen and isstart(result[:3]):
                        if frame > 0:
                            strain = "+"
                            startpos = orf[1]
                            stoppos = startpos + len(result) - 1
                            final_result.append(ORF(result, strain, frame, startpos, stoppos))
                        else:
                            strain = "-"
                            startpos = orf[1]
                            stoppos = startpos - len(result) + 1
                            final_result.append(ORF(result, strain, frame, stoppos, startpos))

                ORFs = []
            else:
                for orf in ORFs:
                    orf[0].append(codon)
    return final_result