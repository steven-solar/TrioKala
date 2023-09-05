### Current manual pipeline:

* Error correction on each member of the trio
    * ie. verkko/verkko_747.sh
* Combine error corrected reads
    * cat together each assembly74*/hifi-corrected.fasta.gz into one verkko/trio-hifi-corrected.fasta.gz
* Generate parental hapmers
    * meryl/run_meryl.sh -> meryl/run_hapmers.sh
* Run trio assembly
    * verkko/verkko_trio.sh
* Color the graph with the reads
    * 2-processGraph/unitig-unrolled-hifi-resolved.gfa (string graph file we want)
    * 2-processGraph/unitig-mapping-1.txt  (map of utigs to node ids in 1/*)
    * 1-buildGraph/paths.gaf (node id paths we want for coloring)
    * Convert paths of node ids to paths of utigs with scripts/replacepaths_nodeID_utigID.py -> utigs.gaf
    * run genochondromatosis/scripts/run_coloring.sh with reads and utigs.gaf
* Get list of denovo utigs (all the ones labeled ‘hap1’ in colors.csv)
* Filter tips and singletons
    * scripts/filter_singletons_tips.py
* Filter low coverage nodes also (either now or at the very end)
* Go through de novos, manually extract paths for two haplotypes (one through denovo and one not, both should be same inherited from parent)
* Run these paths through consensus to get true sequences
* Variant call denovo hap against ref hap to generate vcf file of variants (reduce min aln length)
    * denovo_comparisons/223595_pramef/dipcall/hap_to_hap/run_hap_to_hap_dipcall.sh
* Generate chain file lifting from ref hap to chm13 
    * (denovo_comparisons/223595_pramef/flo/run_flo.sh, edit flo_opts.yaml first) -> run/liftover.chn
* Lift vcf to chm13
    * https://crossmap.readthedocs.io/en/latest/
    * CrossMap.py vcf --chromid a ../flo/hap_to_ref_run/liftover.chn ../dipcall/hap_to_hap/hap_to_hap_dipcall.pair.vcf ../../../data/chm13.fa ./lifted.vcf
* Combine VCFs
* run VEP, summarize output

### Goal Automated Pipeline
* Error correction on each member of the trio
    * ie. verkko/verkko_747.sh
* Combine error corrected reads
    * cat together each assembly74*/hifi-corrected.fasta.gz into one verkko/trio-hifi-corrected.fasta.gz
* Generate parental hapmers
    * meryl/run_meryl.sh -> meryl/run_hapmers.sh
* Run trio assembly
    * verkko/verkko_trio.sh
* Color the graph with the reads
    * 2-processGraph/unitig-unrolled-hifi-resolved.gfa (string graph file we want)
    * 2-processGraph/unitig-mapping-1.txt  (map of utigs to node ids in 1/*)
    * 1-buildGraph/paths.gaf (node id paths we want for coloring)
    * Convert paths of node ids to paths of utigs with scripts/replacepaths_nodeID_utigID.py -> utigs.gaf
    * run genochondromatosis/scripts/run_coloring.sh with reads and utigs.gaf
* Get list of denovo utigs (all the ones labeled ‘hap1’ in colors.csv)
* Filter tips and singletons
    * scripts/filter_singletons_tips.py
* Filter low coverage nodes (<7.5x)
* Hop 1 nodes in each direction from de novo, extract sequence
* Temporarily cut denovo out of the graph
* GraphAlign sequence to "neighborhood" of de novo (not sure how to define this exactly)
* Consensus on top alignment nodes + de novo nodes
* Dipcall for variants, lift to chm13, VEP