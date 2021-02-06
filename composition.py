import seaborn as sns
import matplotlib.pyplot as plt

def plot(DNA: list, ORFs: list):
    length_DNA = [len(i) for i in DNA]
    DNA_length_counter = [length_DNA.count(x) for x in set(length_DNA)]
    length_ORFs = [len(i.seq) for i in ORFs]
    ORFs_length_counter = [length_ORFs.count(x) for x in set(length_ORFs)]

    sns.set_theme(style="whitegrid")
    sns.distplot(DNA_length_counter)
    sns.distplot(ORFs_length_counter)
    plt.show()
    return 0




def composition(DNA: str) -> []:
    nucleotids = "ACGT"
    DNA_freq = {x: DNA.count(x) for x in nucleotids}
    return DNA_freq