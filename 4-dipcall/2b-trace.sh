#!/bin/bash

# Track nodes through verkko steps to rukki contigs
python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/combine_vcf.py \
    all$2compressed.vcf \
    $1

python path_to_node_vcf.py \
    ../../2-processGraph/unitig-unrolled-hifi-resolved.fasta \
    ../1-generate-paths/comparison-paths.txt \
    ../../2-processGraph/unitig-unrolled-hifi-resolved.gfa \
    utig1$2compressed.fasta \
    all$2compressed.compressed.vcf \
    utig1$2compressed.vcf

python join_segments_vcf.py \
    ../../2-processGraph/unitig-unrolled-hifi-resolved.fasta \
    ../../5-untip/unitig-mapping-3.txt \
    ../../2-processGraph/unitig-unrolled-hifi-resolved.gfa \
    utig3$2fasta \
    utig1$2compressed.vcf \
    utig3$2compressed.vcf

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/gfa_to_fasta.py \
    ../../5-untip/unitig-normal-connected-tip.gfa \
    ../../5-untip/unitig-normal-connected-tip.fasta

python join_segments_vcf.py \
    ../../5-untip/unitig-normal-connected-tip.fasta \
    ../../5-untip/unitig-mapping-4.txt \
    ../../5-untip/unitig-normal-connected-tip.gfa \
    utig4$2compressed.fasta \
    utig3$2compressed.vcf \
    utig4$2compressed.vcf

python join_segments_vcf.py \
    ../../6-rukki/unitig-popped-unitig-normal-connected-tip.fasta \
    ../../6-rukki/unitig-popped-unitig-normal-connected-tip.paths.gaf \
    ../../6-rukki/unitig-popped-unitig-normal-connected-tip.noseq.gfa \
    rukki$2compressed.fasta \
    utig4$2compressed.vcf \
    rukki$2compressed.vcf

python rename_vcf_with_layout.py \
    ../../6-layoutContigs/unitig-popped.layout.scfmap \
    rukki$2compressed.vcf > final$2compressed.vcf

# Build chain files for lifting to chm13
minimap2 \
    -t $SLURM_CPUS_PER_TASK \
    -c \
    -x asm20 \
    ../../assembly.fasta \
    /data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.fasta > chm13_to_asm.paf

paf2chain -i chm13_to_asm.paf > chm13_to_asm.chain

minimap2 \
    -t $SLURM_CPUS_PER_TASK \
    -c \
    -x asm20 \
    /data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.fasta \
    ../../assembly.fasta > asm_to_chm13.paf
    
paf2chain -i asm_to_chm13.paf > asm_to_chm13.chain

bash /data/solarsj/hp_compression_mapping/build-chains.sh \
    ../../assembly.fasta \
    assembly.compressed.fasta \
    assembly

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/vcf_to_bed.py \
    final$2compressed.vcf > final$2compressed.bed

bash /data/solarsj/hp_compression_mapping/lift.sh \
    final$2compressed.bed \
    assembly_uncompression.chn \
    ../../assembly.fasta \
    final$2 \
    cu

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/replace_vcf_coords.py \
    final$2compressed.vcf \
    final$2lifted.bed.final > final$2vcf

java -jar ~/picard.jar LiftoverVcf \
    I=final$2vcf \
    O=final$2chm13_to_asm.lifted.vcf \
    CHAIN=chm13_to_asm.chain \
    REJECT=final$2chm13_to_asm.unlifted.vcf \
    R=/data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.fasta

java -jar ~/picard.jar LiftoverVcf \
    I=final$2vcf \
    O=final$2asm_to_chm13.lifted.vcf \
    CHAIN=asm_to_chm13.chain \
    REJECT=final$2asm_to_chm13.unlifted.vcf \
    R=/data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.fasta
