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