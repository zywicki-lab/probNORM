# BAM example file

The example input files are prepared from the Mod-seq experiment downloaded from the SRA database (Accessions: SRR955862, SRR955864). Data preparation was performed as follows. The adapter handling step was performed with cutadapt. First, 3’ end adapters were trimmed, reads with 5’ end adapters were removed, and reads shorter than 20nt were discarded from further analysis. The remaining reads were filtered with FASTX-Toolkit to remove low-quality reads, where <80% of bases have quality >20. Next, obtained reads were mapped with bowtie2 to canonical transcripts of Saccharomyces cerevisiae retrieved with BioMart API v84. Sequences with similarity >0.95 were filtered out with CD-HIT. Unique mappings were used for reactivity profile calculation.

From the obtained dataset the 3 RNAs were filtered: RDN25-1, RDN18-1, SCR1 into example files: example/control.sorted.bam and example/treated.sorted.bam

<br>

# Counts example file

File example/counts.txt contains the three RNAs from Mod-seq experiment above, and two transcripts: 5S in vivo and 5S in vitro from SHAPE-seq 2.0 experiment derived from RNA Mapping DataBase.