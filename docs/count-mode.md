# Counts format file



probNORM input file contains RT polymerase stops counts for both treated and control samples. The file consists of four tab-delimited columns:

1.  Transcript ID

2.  Position

3.  Stops count in the control sample

4.  Stops count in the treated sample


It is possible to upload multiple transcripts in one file. probNORM requires a minimum of 20 positions for each transcript to start the normalization process. At least 20% of those should have both control and treated counts higher than zero.

The example input counts file is provided in **example/counts-input.txt**

## Options
To show full list of available options with their description type:

    probnorm counts -h

### Required

!!! tip ""
    -i, --input

    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The input file listing RT polymerase stops counts for both treated and control samples.<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The file structure is as follows: 1. Transcript ID; 2. Position; 3. Stops count in the control sample;<br>
    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;4. Stops count in the treated sample;

    -o, --output

    &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The name for the probNORM output file.

### Optional

probNORM provides many additional parameters thet allow the detailed adjustment of the algorithm to the analyzed data.

!!! tip ""

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