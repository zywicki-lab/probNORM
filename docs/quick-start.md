# Quick start

The main file of probNORM program is **probnorm**. To quickly run probNORM on provided example files type:

    probnorm bam -t example/treated.sorted.bam -c example/control.sorted.bam -o output.txt

for BAM format input, and:

    probnorm counts -i example/counts-input.txt -o output.txt

for count format input.

!!! tip "Important"

    This command will run probNORM with the default parameters.

probNORM allows for two format of input data: [BAM file](bam-mode.md) or [custom made counts file](). Depending on the input type, the additional options may vary.
