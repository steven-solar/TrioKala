import sys

MAX_SCORE = 999

gene_list = open(sys.argv[1])
alias_file = open(sys.argv[2])
links_file = open(sys.argv[3])

genes = set()
aliases = dict()

for line in gene_list:
    genes.add(line.strip().split()[0])

for line in alias_file:
    line_split = line.strip().split('\t')
    if line_split[1] in genes:
        aliases[line_split[0]] = line_split[1]

printed_genes = set()

for line in links_file:
    line_split = line.strip().split()
    if line_split[0] in aliases and line_split[1] in aliases:
        if int(line_split[-1]) > MAX_SCORE / 2:
            if not aliases[line_split[0]] in printed_genes and not aliases[line_split[1]] in printed_genes:
                if aliases[line_split[0]] != aliases[line_split[1]]:
                    print(aliases[line_split[0]] + " - " + aliases[line_split[1]] + ": " + line_split[-1])
                    printed_genes.add(aliases[line_split[0]])
                    printed_genes.add(aliases[line_split[1]])
