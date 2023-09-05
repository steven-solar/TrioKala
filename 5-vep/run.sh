#!/bin/bash

ml load VEP python samtools

vep \
	--fork ${SLURM_CPUS_PER_TASK} \
	--force_overwrite \
    --symbol \
	--vcf \
	--gff /data/Phillippy/projects/genochondromatosis/steven/data/chm13.draft_v2.0.gene_annotation.gff3.gz \
    --fasta /data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.fasta \
   	-i ../4-dipcall/all-nodes.lifted.sorted.vcf \
    -o all-nodes.vep.vcf

vep \
	--fork ${SLURM_CPUS_PER_TASK} \
	--force_overwrite \
    --symbol \
	--vcf \
	--gff /data/Phillippy/projects/genochondromatosis/steven/data/chm13.draft_v2.0.gene_annotation.gff3.gz \
    --fasta /data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.fasta \
   	-i ../4-dipcall/all-nodes.lifted.windowed.sorted.vcf \
    -o all-nodes.windowed.vep.vcf

vep \
	--fork ${SLURM_CPUS_PER_TASK} \
	--force_overwrite \
    --symbol \
	--vcf \
	--gff /data/Phillippy/projects/genochondromatosis/steven/data/chm13.draft_v2.0.gene_annotation.gff3.gz \
    --fasta /data/Phillippy/projects/genochondromatosis/steven/data/chm13v2.0.fasta \
   	-i ../4-dipcall/all-nodes.lifted.SVs.sorted.vcf \
    -o all-nodes.SVs.vep.vcf
