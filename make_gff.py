genes = [
    ("gene1", 1000, 4998, "+"),
    ("gene2", 7000, 9000, "-"),
    ("gene3", 20000, 22999, "+"),
    ("gene4", 26000, 28999, "-"),
    ("gene5", 31000, 33999, "+"),
    ("gene6", 38000, 40999, "+"),
]
with open("genes.gff", "w") as f:
    f.write("##gff-version 3\n")
    for g, s, e, st in genes:
        f.write(f"chrSim\tsim\tgene\t{s}\t{e}\t.\t{st}\t.\tID={g};Name={g}\n")
        f.write(f"chrSim\tsim\tmRNA\t{s}\t{e}\t.\t{st}\t.\tID={g}.t1;Parent={g}\n")
        f.write(f"chrSim\tsim\texon\t{s}\t{e}\t.\t{st}\t.\tID={g}.e1;Parent={g}.t1\n")
        f.write(f"chrSim\tsim\tCDS\t{s}\t{e}\t.\t{st}\t0\tID={g}.c1;Parent={g}.t1\n")
