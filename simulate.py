import random
random.seed(42)

L = 50000
ref = [random.choice("ACGT") for _ in range(L)]

# ساخت ۲۰ جهش SNP مصنوعی (ژنوتیپ هتروزیگوت)
positions = sorted(random.sample(range(1000, L - 1000), 20))
truth = {}
for p in positions:
    alt = random.choice([b for b in "ACGT" if b != ref[p]])
    truth[p] = alt

# ذخیره ژنوم مرجع
with open("ref.fa", "w") as f:
    f.write(">chrSim\n")
    s = "".join(ref)
    for i in range(0, L, 60):
        f.write(s[i:i+60] + "\n")

# ذخیره جهش‌های واقعی (برای مقایسه بعدی)
with open("truth.txt", "w") as f:
    for p, a in truth.items():
        f.write(f"chrSim\t{p+1}\t{ref[p]}\t{a}\n")

# دو هاپلوتایپ: یکی مرجع، یکی دارای جهش‌ها
hap1 = ref[:]
hap2 = ref[:]
for p, a in truth.items():
    hap2[p] = a
haps = ["".join(hap1), "".join(hap2)]

comp = {"A": "T", "C": "G", "G": "C", "T": "A"}
def rc(s):
    return "".join(comp[c] for c in reversed(s))

def mutate_err(s):
    s = list(s)
    qual = []
    for i in range(len(s)):
        if random.random() < 0.005:
            s[i] = random.choice("ACGT")
            qual.append("5")
        else:
            qual.append("I")
    return "".join(s), "".join(qual)

N = 25000          # تعداد جفت‌ریدها (~۵۰۰x پوشش روی ۵۰ کیلوباز)
RL = 100
with open("sample_R1.fastq", "w") as r1, open("sample_R2.fastq", "w") as r2:
    for n in range(N):
        h = random.choice(haps)
        frag = random.randint(250, 400)
        start = random.randint(0, L - frag)
        fragment = h[start:start + frag]
        s1, q1 = mutate_err(fragment[:RL])
        s2, q2 = mutate_err(rc(fragment)[:RL])
        r1.write(f"@read{n}/1\n{s1}\n+\n{q1}\n")
        r2.write(f"@read{n}/2\n{s2}\n+\n{q2}\n")

print("Done: ref.fa, truth.txt, sample_R1.fastq, sample_R2.fastq")
