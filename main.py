import parse_gff
import parse_fasta
import ORF
import composition

seq = (str(parse_fasta.parse()))
print("Sequence size -> {}\n".format(len(seq)))
db = parse_gff.parse()

genes = list(db.features_of_type('gene'))
print("Gene size -> {}\n".format(len(genes)))

start_stop = []
for gen in genes:
    start_stop.append([gen.start, gen.stop])

subseq = []
for ss in start_stop:
    subseq.append(seq[ss[0] - 1:ss[1]])

ORFs = ORF.ORF_finder_pro(seq, 10)
print("ORF size -> {}\n".format(len(ORFs)))

counter = 0
for orf in ORFs:
    if orf.strain == '-':
        if ORF.revcompdna(orf.seq) in subseq:
            counter += 1
    if orf.strain == '+':
        if orf.seq in subseq:
            counter += 1

print("Find subseq in ORFs -> {} / {}\n".format(counter, len(genes)))

composition.plot(subseq, ORFs)