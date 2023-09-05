import sys

layout = open(sys.argv[1], 'r')
vcf = open(sys.argv[2], 'r')

layout_dict = dict()
for line in layout:
    if line[:4] == 'path':
        line_split = line.strip().split(' ')
        layout_dict[line_split[2]] = line_split[1]

missing_nodes = set()
for line in vcf:
    if line[0] == '#':
        print(line.strip())
    else:
        line_split = line.strip().split('\t')
        try:
            line_split[0] = layout_dict[line_split[0]]
            print('\t'.join(line_split))
        except:
            missing_nodes.add(line_split[0])