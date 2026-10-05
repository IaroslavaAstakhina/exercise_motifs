from Bio import motifs
from Bio.Seq import Seq
from Bio import SeqIO
import random
from itertools import product


#______________ TASK 1  -  Warm-up functions__________________________________________________________________
# Write the following five functions: count matrix, score, consensus,  Hamming distance and total distance
BASES = "ACGT"
lecture_dna = [
    "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
    "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
    "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
    "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
    "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
    "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
]
def count_matrix(motifs):
    l = len(motifs[0])
    counts = {base: [0] * l for base in BASES}
    for motif in motifs:
        for i, base in enumerate(motif):
            counts[base][i] += 1
    return counts

def score(motifs):
    counts = count_matrix(motifs)
    l = len(motifs[0])
    return sum(max(counts[base][i] for base in BASES) for i in range(l))

def consensus(motifs):
    counts = count_matrix(motifs)
    l = len(motifs[0])
    result = ""
    for i in range(l):
        best = "A"
        for base in BASES:
            if counts[base][i] > counts[best][i]:
                best = base
        result += best
    return result

def hamming_distance(a,b):
    return sum(1 for x, y in zip(a, b) if x != y)

def total_distance(pattern, sequences):
    l = len(pattern)
    total = 0
    for seq in sequences:
        total += min(
            hamming_distance(pattern, seq[i:i + l])
            for i in range(len(seq)-l + 1)
        )
    return total


red = ["TAAGTT", "TGAATT", "GGAGTG", "CGAGTC", "TGTGTG", "TGGGTC"]  # slide 19
best = ["AGATAG", "AGATAG", "AGATAG", "AGACAG", "AGATAG", "AGGTAG"]


print(score(red))                                # 26
print(consensus(best), score(best))              # AGATAG 34
print(hamming_distance("TAAGTT", "TGAATT"))      # 2
print(total_distance("TGCGTT", lecture_dna))     # 13


#_____________________ TASK 2 - MotifProfile class ___________________________________________________________________________

class MotifProfile:
    def __init__(self, motifs, pseudocount=1):
        self.l = len(motifs[0])
        t = len(motifs)
        counts = count_matrix(motifs)
        self.ppm = {
            base: [(c + pseudocount) / (t + 4 * pseudocount) for c in counts[base]]
            for base in BASES
        }

    def lmer_probability(self, lmer):
        p = 1.0
        for i, base in enumerate(lmer):
            p *= self.ppm[base][i]
        return p

    def most_probable_lmer(self, sequence):
        best_lmer = None
        best_p = -1.0
        for i in range(len(sequence) - self.l + 1):
            lmer = sequence[i:i + self.l]
            p = self.lmer_probability(lmer)
            # строго больше, поэтому при ничьей остаётся самое левое окно
            if p > best_p:
                best_p = p
                best_lmer = lmer
        return best_lmer

    def consensus(self):
        result = ""
        for i in range(self.l):
            best = "A"
            for base in BASES:
                if self.ppm[base][i] > self.ppm[best][i]:
                    best = base
            result += best
        return result

profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
print ( " ____________ TASK 2  ____________")
print(profile.l)                                       # 7
print(profile.ppm["A"])                                # [0.5, 0.25, 0.125, 0.125, 0.25, 0.125, 0.5]
print(profile.consensus())                             # ATGCGTA
print(round(profile.lmer_probability("ATGCGTA"), 4))   # 0.0122

two = MotifProfile(["GTAC", "TTAA"])
print(two.most_probable_lmer("ACTGGATGACCC"))          # TGAC
print(round(two.lmer_probability("TGAC"), 4))          # 0.0093

bio = motifs.create([Seq(site) for site in ["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"]])
bio.pseudocounts = 1
print( "Porovnani s BioPythonem ")
print(bio.consensus)     # ATGCGTA
print(bio.pwm["A"])


#_______________________ TASK 3 -  MotifFinder class and random module ___________________________________________

class MotifFinder:
    def __init__(self, sequences, l, seed=None):
        self.sequences = sequences
        self.l = l
        self.rng = random.Random(seed)
        self.windows = [
            [seq[i:i + l] for i in range(len(seq) - l + 1)]
            for seq in sequences
        ]

    def total_distance(self, pattern):
        return sum(
            min(hamming_distance(pattern, w) for w in seq_windows)
            for seq_windows in self.windows
        )

    def median_string(self):
        best_pattern = None
        best_dist = float("inf")
        for letters in product(BASES, repeat=self.l):
            pattern = "".join(letters)
            d = self.total_distance(pattern)
            if d < best_dist:
                best_dist = d
                best_pattern = pattern
        return best_pattern, best_dist


    def randomized_search(self):
        motifs = [self.rng.choice(w) for w in self.windows]
        current_score = score(motifs)
        while True:
            profile = MotifProfile(motifs, pseudocount=1)
            new_motifs = [profile.most_probable_lmer(s) for s in self.sequences]
            new_score = score(new_motifs)
            if new_score > current_score:
                motifs, current_score = new_motifs, new_score
            else:
                return motifs, current_score

    def best_of(self, runs):
        best_motifs, best_score = None, -1
        for _ in range(runs):
            motifs, sc = self.randomized_search()
            if sc > best_score:
                best_motifs, best_score = motifs, sc
        return best_motifs, best_score


#Check on the lecture data. MotifFinder(lecture_dna, 6) must give the median string AGATAG with total
# distance 2 (the answer from the lecture), and best_of(100) should find motifs with the same consensus and score 34
finder = MotifFinder(lecture_dna, 6, seed=1)
print( "______   TASK 3 _______")
print(finder.median_string())            # ('AGATAG', 2)
m, sc = finder.best_of(100)
print(consensus(m), sc)                  # AGATAG 34

#Find the planted motif. Run median_string() and best_of(100) on planted_motif.fasta with l = 7.
# Print the l-mers that best_of(100) picked in each sequence and compare their consensus with the median string.

sequences = [str(r.seq) for r in SeqIO.parse("planted_motif.fasta", "fasta")]

planted = MotifFinder(sequences, 7, seed=1)
median, dist = planted.median_string()
print("median string:", median, dist)    # GCTAAAG 10

best_motifs, best_score = planted.best_of(100)
for m in best_motifs:
    print(m)
print("consensus:", consensus(best_motifs), "score:", best_score)   # score 60
print("agree:", consensus(best_motifs) == median)
print("_______________________________________________")
# How many restarts do you need? On planted_motif.fasta
def success_rate_single(finder, target, n=200):
    ok = 0
    for _ in range(n):
        m, _ = finder.randomized_search()
        if consensus(m) == target:
            ok += 1
    return ok / n


def success_rate_best_of(finder, target, R, repeats=50):
    ok = 0
    for _ in range(repeats):
        m, _ = finder.best_of(R)
        if consensus(m) == target:
            ok += 1
    return ok / repeats


def restart_table(finder, target):
    p = success_rate_single(finder, target)
    print(f"p (one run) = {p:.3f}")
    print(f"{'R':>3} {'measured':>9} {'predicted':>10}")
    for R in (5, 10, 20, 50):
        measured = success_rate_best_of(finder, target, R)
        predicted = 1 - (1 - p) ** R
        print(f"{R:>3} {measured:>9.2f} {predicted:>10.2f}")


print("planted:")
restart_table(planted, median)

print("lecture:")
restart_table(MotifFinder(lecture_dna, 6, seed=2), "AGATAG")



l = planted.l
found = motifs.create([Seq(m) for m in best_motifs])
found.pseudocounts = 1
pssm = found.pssm


def scan(fpr):
    threshold = pssm.distribution().threshold_fpr(fpr)
    total = found_sites = reverse = 0
    for idx, sequence in enumerate(sequences):
        for position, hit_score in pssm.search(Seq(sequence), threshold=threshold):
            total += 1
            if position < 0:
                reverse += 1
                start = len(sequence) + position

            elif sequence[position:position + l] == best_motifs[idx]:
                found_sites += 1
    expected = len(sequences) * (len(sequences[0]) - l + 1) * 2 * fpr
    print(f"fpr={fpr}: threshold={threshold:.2f}, hits={total}, "
          f"found sites={found_sites}, reverse={reverse}, expected by chance={expected:.1f}")


scan(0.01)
scan(0.001)
