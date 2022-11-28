# BAM format file

A BAM file (.bam) is the binary version of a SAM file, which is a tab-delimited text file that contains sequence alignment data. These formats are described in details on the SAM Tools web site: [http://samtools.github.io/hts-specs/](http://samtools.github.io/hts-specs/).

The example input BAM files is provided in **example/treated.sorted.bam** and **example/control.sorted.bam**

## Options
To show full list of available options with their description type:

    probnorm bam -h

### Required

!!! tip ""
    -t, --treated

    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Indexed BAM](#prepare-input-bam-files) file for treated sample. If provided file is not indexed, probNORM will prepare it with [SAMTools](http://www.htslib.org/)<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(requires SAMtools installation and increases the execution time)

    -c, --control

    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[Indexed BAM](#prepare-input-bam-files) file for control sample. If provided file is not indexed, probNORM will indexed it with [SAMTools](http://www.htslib.org/)<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;(requires SAMtools installation and increases the execution time)

    -o, --output

    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The name for the probNORM output file.

### Optional

probNORM provides many additional parameters thet allow the detailed adjustment of the algorithm to the analyzed data.

!!! tip ""

    -id

    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;probNORM can normalize only selected transcript. Quote transcripts ids in the comma separated list e.g. "RDN18,RDN25".<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Default: all transcripts are normalized.

    -v, --coverage

    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Coverage filtering. Only positions with coverage higher than specified are normalized.

    -p, --pvalue

    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;P-value is the probability that a nucleotide belongs to the background distribution and is not statistically significant.<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Range [0-1]. All positions with a p-value higher than the provided one are rejected from the result.<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Default: 1 -> showing all positions.

    -s, --transcript-size

    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Set the percentage of transcript covered with reactive nucleotides. Those positions are further use to set normalization<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;parameters. Default: 20 percent of provided transcript.

    -f, --constrain-files

    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The constrain files with normalized reactivities for RNAfold [ViennaRNA] and Fold [RNAStructure] will be prepared<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;for each normalized transcript. Files will be saved in <output-file-name\>_constrains.


## Prepare input BAM files

By default the input BAM files should be sorted and indexed for the max optimization of the algorithm. User can do it by running following commands ([SAMtools](http://www.htslib.org/) software is required for this steps).

    samtools sort input.bam -o input_sorted.bam
    samtools index input_sorted.bam

probNORM is able to sort and index provided files, however this step still requires the [SAMtools](http://www.htslib.org/) and significatly increases the execution time.