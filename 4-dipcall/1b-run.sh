#!/bin/bash

prefix=calls/$1_$2/$1_$2

touch $1.lst
echo $1 > $1.lst
seqtk subseq \
    ../2-paths-consensus/7-unitigConsensus/assembly.fasta \
    $1.lst \
    > fastas/$1.fasta
samtools faidx fastas/$1.fasta
rm $1.lst
    
touch $2.lst
echo $2 > $2.lst
seqtk subseq \
    ../2-paths-consensus/7-unitigConsensus/assembly.fasta \
    $2.lst \
    > fastas/$2.fasta
samtools faidx fastas/$2.fasta
rm $2.lst

mkdir -p calls/$1_$2

/home/solarsj/dipcall.kit.copy/run-dipcall \
    -t $SLURM_CPUS_PER_TASK \
    $1_$2 \
    fastas/$2.fasta \
    fastas/$1.fasta \
    fastas/$2.fasta > $prefix.mak

make -j2 -f $prefix.mak
mv $1_$2.* calls/$1_$2/

gunzip < $prefix.pair.vcf.gz > $prefix.pair.vcf
bgzip < $prefix.pair.vcf > $prefix.pair.vcf.bgz
tabix -p vcf $prefix.pair.vcf.bgz

# window the vcf, or just copy if no window exists
if [ -s ../3-generate-windows/lifts/chains/$2.lifted.bed.final ]; then
    bcftools view \
        -R ../3-generate-windows/lifts/chains/$2.lifted.bed.final \
        $prefix.pair.vcf.bgz > $prefix.pair.windowed.vcf
    bgzip < $prefix.pair.windowed.vcf > $prefix.pair.windowed.vcf.bgz
else
    cp $prefix.pair.vcf $prefix.pair.windowed.vcf
    cp $prefix.pair.vcf.bgz $prefix.pair.windowed.vcf.bgz
fi

tabix -p vcf $prefix.pair.windowed.vcf.bgz

# compress the windowed vcf so we can track through nodes and later lift to ref
python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/vcf_to_bed.py \
    $prefix.pair.windowed.vcf > $prefix.pair.windowed.bed

bash /data/solarsj/hp_compression_mapping/lift.sh \
    $prefix.pair.windowed.bed \
    ../3-generate-windows/lifts/chains/$2_compression.chn \
    fastas/$2.fasta \
    $prefix.pair.windowed \
    uc

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/replace_vcf_coords.py \
    $prefix.pair.windowed.vcf \
    $prefix.pair.windowed.lifted.bed.final > $prefix.pair.windowed.compressed.vcf 

# compress the unwindowed vcf so we can track through nodes and later lift to ref
python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/vcf_to_bed.py \
    $prefix.pair.vcf > $prefix.pair.bed

bash /data/solarsj/hp_compression_mapping/lift.sh \
    $prefix.pair.bed \
    ../3-generate-windows/lifts/chains/$2_compression.chn \
    fastas/$2.fasta \
    $prefix.pair \
    uc

python /data/Phillippy/projects/genochondromatosis/steven/python_scripts/replace_vcf_coords.py \
    $prefix.pair.vcf \
    $prefix.pair.lifted.bed.final > $prefix.pair.compressed.vcf 

# Call SVs
minimap2 \
    -t $SLURM_CPUS_PER_TASK \
    -c \
    -a \
    -x asm20 \
    fastas/$2.fasta fastas/$1.fasta > calls/$1_$2/$1_to_$2.sam

samtools view -b calls/$1_$2/$1_to_$2.sam | samtools sort -o calls/$1_$2/$1_to_$2.bam
samtools index calls/$1_$2/$1_to_$2.bam
svim-asm haploid calls/$1_$2/ calls/$1_$2/$1_to_$2.bam fastas/$2.fasta
mv calls/$1_$2/variants.vcf $prefix.SVs.vcf

# Lift uncompressed vcfs to ref
minimap2 \
    -t $SLURM_CPUS_PER_TASK \
    -c \
    -x asm20 \
    fastas/$2.fasta \
    /data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.fasta > $prefix.paf
    
paf2chain -i $prefix.paf > $prefix.chain

java -jar ~/picard.jar LiftoverVcf \
    I=$prefix.pair.vcf \
    O=$prefix.pair.lifted.vcf \
    CHAIN=$prefix.chain \
    REJECT=$prefix.pair.unlifted.vcf \
    R=/data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.fasta

java -jar ~/picard.jar LiftoverVcf \
    I=$prefix.pair.windowed.vcf \
    O=$prefix.pair.windowed.lifted.vcf \
    CHAIN=$prefix.chain \
    REJECT=$prefix.pair.windowed.unlifted.vcf \
    R=/data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.fasta

java -jar ~/picard.jar LiftoverVcf \
    I=$prefix.SVs.vcf \
    O=$prefix.SVs.lifted.vcf \
    CHAIN=$prefix.chain \
    REJECT=$prefix.SVs.unlifted.vcf \
    R=/data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.fasta
