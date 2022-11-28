## Format

<br>
The file contains full information about the normalized transcript/s. It consists of nine tab separated columns:

| Column name | Description |
|-------------|-------------|
| transcript_id | ID of normalized transcript, the same as in the input file |
| position | Position in transcript |
| stops_treated | Stops count in the treated sample: from input counts file or calculated from BAM file |
| stops_control | Stops count in the control sample: from input counts file or calculated from BAM file |
| stops_norm_control | Normalized stops count in the control sample. Stops are normalized by incorporating the normalization factor (nf). Check [The overview of probNORM method](index.md#the-overview-of-probnorm-method) point B for more information |
| reactivity | Reactivity, calculated based on the normalized control stops. Check [The overview of probNORM method](index.md#the-overview-of-probnorm-method) point D for more information |
| fold_change | The ratio between stops counts in control and treated sample |
| p_value | P-value indicates the probability of nucleotide at a given position being a part of the background, not statistically significant. Check [The overview of probNORM method](index.md#the-overview-of-probnorm-method) point C for more information |
| passed_quality_filter | Quality filter (Y - yes / N - no). Transcript positions that exceed the filtering step are those with stops count higher than zero (both control and treated samples), without any missing parameters, and with proper coverage value (when a local script is determining the stops counts |


=== "probnorm counts -i example/counts-input.txt -o output.txt"

    ```
    transcript_id	position	stops_treated	stops_control	stops_norm_control	reactivity	fold_change	p_value	passed_quality_filter
    RDN18-1	1	3095.0	3472.0	2669.1000000000004	1.0632124544542494	0.2135860512052699	0.37310634695017253	Y
    RDN18-1	2	2029.0	1148.0	882.5250000000001	2.5274855472882036	1.2010598126290937	0.03438625350046609	Y
    RDN18-1	3	315.0	360.0	276.75	0.09548691331973486	0.18676851160572655	0.38858771448505425	Y
    RDN18-1	4	264.0	405.0	311.34375	0.0	-0.23797038886541122	0.6407954148840493	Y
    RDN18-1	5	139.0	171.0	131.45625	0.018832141238058788	0.08050214738573189	0.45145693582080115	Y
    ...
    RDN18-1	1776	0	0	0.0	0.0	0	0.5	N
    RDN18-1	1777	0	0	0.0	0.0	0	0.5	N
    RDN18-1	1778	0	0	0.0	0.0	0	0.5	N
    RDN18-1	1779	0	0	0.0	0.0	0	0.5	N
    RDN18-1	1780	25.0	9.0	6.91875	0.04513784971143676	1.8533447778805348	0.002490274610317811	Y
    ```


## Summary information

After each use of probNORM the summary of run will be shown. It contains such informations as:

- input file type
- input and output file names
- parameters thresholds: coverage, p-value, reactive positions
- statictics about normalized transcripts
  
See the example below.

=== "SUMMARY"
    === "BAM input"

        ```
        ***** SUMMARY *****

            input mode: BAM
            input file/s: control: data/SRR955862_sort.sorted.bam treated:data/SRR955864_sort.sorted.bam
            output file: test.output
            min coverage: 0
            max p-value: 1.0
            min reactive positions per transcript: 20%
            selected transcripts:  all
            total number of input transcripts: 7127
            transcripts omitted due to low reactivity: 6887
            transcripts normalized: 240

        *******************
        ```
    === "COUNTS input"

        ```
        ***** SUMMARY *****

            input mode: COUNTS
            input file/s: data/counts-input.txt
            output file: test.output
            max p-value: 1.0
            min reactive positions per transcript: 20%
            total number of input transcripts: 2
            transcripts omitted due to low reactivity: 0
            transcripts normalized: 2

        *******************
        ```