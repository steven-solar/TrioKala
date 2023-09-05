#!/usr/bin/env python
from __future__ import print_function
from __future__ import division

import sys
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_vcf(in_vcf_fp, out_vcf):
    in_vcf = open(in_vcf_fp, 'r')
    vcf_dict = dict()
    for line in in_vcf:
        if line[0] == '#':
            out_vcf.write(line)
        else:
            line_split = line.strip().split('\t')
            if 'denovo' in line_split[0]:
                node = line_split[0].split('_')[0]
            else:
                node = line_split[0]
            if node in vcf_dict:
                vcf_dict[node].append(line_split)
            else:
                vcf_dict[node] = [line_split]
    return vcf_dict

def rev_comp(seq):
    comps = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    rev = seq[::-1]
    out = ''
    for b in rev:
        out += comps[b]
    return out

def my_segment(s, flag):
    if flag == '-':
        return s.reverse_complement()
    elif flag == '+':
        return s
    else:
        sys.exit(1)
        return ''

def swap(c):
    if c == '-':
        return '+'
    elif c == '+':
        return '-'
    assert False

def split_segment_desc(s_d):
    return (s_d[:-1], s_d[-1])

def read_cigars(gfa_fn):
    cigars = dict()
    with open(gfa_fn) as gfa:
        for l in gfa:
            if l[0] == 'L':
                l=re.sub('tig0+', '', l)
                sp=l.split('\t')
                cigars['%s%s|%s%s' % (sp[1], sp[2], sp[3], sp[4])] = sp[5]
                cigars['%s%s|%s%s' % (sp[3], swap(sp[4]), sp[1], swap(sp[2]))] = sp[5]
    return cigars

def get_cigar(cigars, s_d_1, s_d_2):
    q = s_d_1 + "|" + s_d_2
    return cigars[q] if q in cigars else None

def trim_len(cigar, second=True):
    pos = 0
    for match in re.finditer('(\d+)M', cigar):
        pos += int(match.group(1))
    for match in re.finditer('(\d+)I' if second else '(\d+)D', cigar):
        pos += int(match.group(1))
    return pos

def make_sequence(cigars, contig_dict, segment_order, trim_flanks_to):
    print(segment_order)
    segment_descs = [s.strip() for s in segment_order.split(',')]
    print(segment_descs)
    segments=[my_segment(contig_dict[s_d[:-1]], s_d[-1]) for s_d in segment_descs]

    tot = 0
    if trim_flanks_to < 0:
        final = segments[0]
        tot += len(segments[0])
    else:
        final = segments[0][-trim_flanks_to:]
        tot += trim_flanks_to
        print('Here taking', trim_flanks_to)

    for i in range(0, len(segments)-1):
        cigar = get_cigar(cigars, segment_descs[i], segment_descs[i+1])
        print(segment_descs[i], segment_descs[i+1])
        if (cigar is None):
           print("Cigar not known between %s and %s\n"%(segment_descs[i], segment_descs[i+1]))
        assert cigar is not None
        print('================== Basic step')
        print(segment_descs[i], segment_descs[i+1], cigar)
        pos = trim_len(cigar)
        print('From cigar', pos)
        l_before = tot
        final += segments[i+1][pos:]
        tot += (len(segments[i+1]) - pos)
        print('Added coordinates: [%d - %d)' % (l_before, tot))

    if trim_flanks_to > 0 and len(segments[-1]) > trim_flanks_to:
        final = final[:-(len(segments[-1]) - trim_flanks_to)]
        print('Here trimming', len(segments[-1]) - trim_flanks_to)

    return final, False
    #cigar = read_cigar(sys.argv[3], segment_descs[-1], segment_descs[0])
    #if cigar is not None:
    #    print('================== Circular step')
    #    print(segment_descs[-1], segment_descs[0], cigar)
    #    pos = trim_len(cigar)
    #    print('From cigar', pos)
    #    final = final[pos:]
    #    return final, True
    #else:
    #    return final, False

def make_line(cigars, contig_dict, segment_order, trim_flanks_to, oldvcf_line_split, gaps):
    if 'denovo' in oldvcf_line_split[0]:
        oldvcf_node = oldvcf_line_split[0].split('_')[0]
    else:
        oldvcf_node = oldvcf_line_split[0]
    oldvcf_pos, oldvcf_ref, oldvcf_alt = int(oldvcf_line_split[1]), oldvcf_line_split[3], oldvcf_line_split[4]
    print(segment_order)
    segment_descs = [s.strip() for s in segment_order.split(',')]
    print(segment_descs)
    segments=[my_segment(contig_dict[s_d[:-1]], s_d[-1]) for s_d in segment_descs]
    print(segments)
    print(gaps)
    tot = 0
    if trim_flanks_to < 0:
        final = segments[0]
        tot += len(segments[0])
    else:
        final = segments[0][-trim_flanks_to:]
        tot += trim_flanks_to
        print('Here taking', trim_flanks_to)
    # print(segment_descs)
    for i in range(0, len(segments)-1):
        use_gap = False
        cigar = get_cigar(cigars, segment_descs[i], segment_descs[i+1])
        # print(segment_descs[i], segment_descs[i+1])
        if (cigar is None):
            if (segment_descs[i], segment_descs[i+1]) in gaps:
                pos = -1 * 1 # gaps[(segment_descs[i], segment_descs[i+1])]
                use_gap = True
            else:
                print("Cigar or gap not known between %s and %s\n"%(segment_descs[i], segment_descs[i+1]))
        assert use_gap or cigar is not None
        print('================== Basic step')
        if not use_gap:
            print(segment_descs[i], segment_descs[i+1], cigar)
            pos = trim_len(cigar)
        else:
            print(segment_descs[i], segment_descs[i+1], 'N' + str(-1 * pos) + 'N')
        if i == 0 and segment_descs[i][:-1] == oldvcf_node:
            if segment_descs[i][-1] == '+':
                newvcf_pos = oldvcf_pos
                newvcf_ref, newvcf_alt = oldvcf_ref, oldvcf_alt
            if segment_descs[i][-1] == '-':
                newvcf_pos = len(segments[i]) - oldvcf_pos
                newvcf_ref, newvcf_alt = rev_comp(oldvcf_ref), rev_comp(oldvcf_alt)
        if segment_descs[i+1][:-1] == oldvcf_node:
            if segment_descs[i+1][-1] == '+':
                newvcf_pos = tot - pos + oldvcf_pos
                newvcf_ref, newvcf_alt = oldvcf_ref, oldvcf_alt
            if segment_descs[i+1][-1] == '-':
                newvcf_pos = tot - pos + len(segments[i+1]) - oldvcf_pos
                newvcf_ref, newvcf_alt = rev_comp(oldvcf_ref), rev_comp(oldvcf_alt)
        print('From cigar', pos)
        l_before = tot
        if use_gap:
            final += 'N' * (-1 * pos) + segments[i+1]
        else:
            final += segments[i+1][pos:]
        tot += (len(segments[i+1]) - pos)
        print('Added coordinates: [%d - %d)' % (l_before, tot))
    print(segment_descs[0][:-1] == oldvcf_node)
    print(len(segment_descs) == 1)
    if segment_descs[0][:-1] == oldvcf_node and len(segment_descs) == 1:
        if segment_descs[0][-1] == '+':
            newvcf_pos = oldvcf_pos
            newvcf_ref, newvcf_alt = oldvcf_ref, oldvcf_alt
        if segment_descs[0][-1] == '-':
            newvcf_pos = len(segments[0]) - oldvcf_pos
            # Here is the problem, segments are in HPC space and vcf positions are in post-consensus uncompressed space
            # Lift vcf pos to HPC space using chains from step 3-windows/
            print('subtract: ')
            print(len(segments[0]))
            print(oldvcf_pos)
            print(newvcf_pos)
            newvcf_ref, newvcf_alt = rev_comp(oldvcf_ref), rev_comp(oldvcf_alt)
        
    if trim_flanks_to > 0 and len(segments[-1]) > trim_flanks_to:
        final = final[:-(len(segments[-1]) - trim_flanks_to)]
        print('Here trimming', len(segments[-1]) - trim_flanks_to)
    
    return newvcf_pos, newvcf_ref, newvcf_alt

def transform_dir_seg(s):
    if s[0] == '>':
        return s[1:] + '+'
    elif s[0] == '<':
        return s[1:] + '-'
    assert(False)

parser = argparse.ArgumentParser(description="Join sequences of GFA nodes")
parser.add_argument("segments", help="FASTA file with graph segment sequences")
parser.add_argument("paths", help="File detailing paths (name id1[+-],id2[+-],...)")
parser.add_argument("gfa", help="GFA file specifying graph structure")
parser.add_argument("result", help="Output FASTA file")
parser.add_argument("vcf", help="Input vcf file to be transformed")
parser.add_argument("out_vcf", help="Output vcf file")

parser.add_argument("--trim-flanks-to", type=int, default=-1, help="Trim flanking segments to specified value (disabled if negative)")
args = parser.parse_args()

if args.trim_flanks_to > 0:
    print('Will trim flanking segments to %dbp' % args.trim_flanks_to)

print('Collecting sequences from', args.segments)
contig_dict = dict()
for r in SeqIO.parse(open(args.segments, "r"), "fasta"):
    contig_dict[re.sub('tig0+', '', r.name)] = r.seq

print('Reading links from', args.gfa)
cigars = read_cigars(args.gfa)

print('Reading vcf', args.vcf)
out_vcf = open(args.out_vcf, 'w')
vcf_dict = read_vcf(args.vcf, out_vcf)

directed_seg_pattern=re.compile(r"[<>]\w+-\d+")
gap_pattern=re.compile(r"(?=([<>]\w+-\d+)\[N(\d+)N\]([<>]\w+-\d+))")

records = []
i = 0
for line in open(args.paths, 'r'):
    i += 1
    print('=======================================')
    line_split = line.split()
    if line_split[0] == 'name':
        continue
    l = line_split[1]
    print('Processing ', l)
    if len(l) == 0:
        continue

    if l[0] == '>' or l[0] == '<':
        nodes = ','.join([transform_dir_seg(s) for s in directed_seg_pattern.findall(l)])
        gap_regexs = gap_pattern.findall(l)
        gaps = dict()
        for g in gap_regexs:
            g_list = list(g)
            gaps[(transform_dir_seg(g[0]), transform_dir_seg(g[2]))] = int(g[1])
    for oldvcf_node in [s.strip() for s in nodes.split(',')]:
        # print(oldvcf_node[:-1])          
        if oldvcf_node[:-1] in vcf_dict:
            print('MATCH')
            for vcf_line_split in vcf_dict[oldvcf_node[:-1]]:
                newvcf_pos, newvcf_ref, newvcf_alt = make_line(cigars, contig_dict, nodes, args.trim_flanks_to, vcf_line_split, gaps)
                name = line_split[0]
                print('line:')
                print(name)
                print(str(newvcf_pos))
                print(vcf_line_split[2])
                print(vcf_line_split[3] + ' ' + newvcf_ref)
                print(vcf_line_split[4])
                print(newvcf_alt)
                print(vcf_line_split[5:])
                print('\t'.join([name, str(newvcf_pos)] + [vcf_line_split[2]] + [newvcf_ref, newvcf_alt] + vcf_line_split[5:]))
                out_vcf.write('\t'.join([name, str(newvcf_pos)] + [vcf_line_split[2]] + [newvcf_ref, newvcf_alt] + vcf_line_split[5:]) + '\n')

    #seq, circular = make_sequence(cigars, contig_dict, l, args.trim_flanks_to)
    #name = 'curated%d length=%d circular=%s' % (i, len(seq), 'true' if circular else 'false')
    name = line_split[0]
    print('Writing', name)
    #records.append(SeqRecord(seq, id=name, description=''))

#SeqIO.write(records, open(args.result, 'w'), 'fasta')
