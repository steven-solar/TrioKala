from Bio import SeqIO
import subprocess
import shlex

"""
Iterate through all vcfs
Grab compressed de novo node sequence
Grab nondenovo uncompressed consensus path sequence
Compress path sequence
Align de novo node to path - get coordinates in bed format
Run my mapper to uncompress the bed file
filter for variants within this range
"""

def get_record(id, record_dict):
    for r_id in record_dict:
        if r_id == id:
            return record_dict[r_id]
    return None

pairs_file = open('/data/Phillippy/projects/genochondromatosis/steven/verkko/assembly_trio/analysis/compare_haps/pairs.txt', 'r')
compressed_fasta_file = '/data/Phillippy/projects/genochondromatosis/steven/verkko/assembly_trio/2-processGraph/unitig-unrolled-hifi-resolved.fasta'
paths = open('/data/Phillippy/projects/genochondromatosis/steven/verkko/assembly_trio/analysis/consensus_haps/6-layoutContigs/consensus_paths.txt', 'r')

compressed_records = SeqIO.to_dict(SeqIO.parse(compressed_fasta_file, 'fasta'))

for line in pairs_file:
    line_split = line.strip().split('\t')
    denovo_id = line_split[0]
    denovo = line_split[1]
    nondenovo = line_split[2]
    print(denovo + '_' + nondenovo)
    denovo_record = compressed_records[denovo_id]
    SeqIO.write(denovo_record, 'denovo_seqs/' + denovo_id + '.fasta', 'fasta')
    SeqIO.write(path_record, 'path_seqs/' + denovo_id + '.fasta', 'fasta')

