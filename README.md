# pyrimid
Julia (https://julialang.org) code to count nucleotides in RNA-Seq reads aligned to a reference genome.
The software requires the ArgParse.jl package to be installed as well as the BioSequences, BioAlignments, XAM and FASTX packages from BioJulia (https://biojulia.net).

Pyrimid is designed for analyses of RNA editing and uses a window to scan the aligned reads, counting the aligned nucleotides only if all nucleotides in the window match the reference (allowing for potential RNA editing).
In this way Pyrimid avoids counting mismatches due to mis-mapped reads or sequence errors (e.g. indels in homopolymer runs) as putative editing events.

To run the software, install julia (instructions here: https://julialang.org/downloads/) and the required packages (instructions here: https://datatofish.com/install-package-julia/)

Pyrimid.jl requires two arguments:

1. A reference genome sequence in FASTA format
2. RNA-Seq reads aligned to this reference in a .bam file. The .bam file should be sorted by name (e.g using samtools sort -n)

usage: pyrimid.jl [--mapQ MAPQ] [--baseQ BASEQ] [-w CONTEXTWINDOW]
                  [--single-end] [--paired-end] [-f] [-u] [-o OUTFILE]
                  reference bam

optional arguments: -w, width of window to count in (default is 6 nt); -f, flag indicating reads (or first read in a pair) are the forward RNA strand;
-u, flag indicating to count U to C mismatches as editing events; -o, output file to write results to; --single-end, flag indicating BAM file contains unpaired reads

The default settings assume reads are paired-end, with the reverse read corresponding to the forward RNA strand (this is the case for standard Illumina paired-end RNA-Seq libraries.
These assumptions can be over-ridden with the -f and --single-end flags.
