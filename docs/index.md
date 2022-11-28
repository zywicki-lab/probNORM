# Welcome to probNORM

A new structural probing signal calculation method that eliminates read distribution bias and prevents reactivity underestimation. It is based on the analysis of background RT stops in treated and control samples of a single replicate and enables statistical discrimination of the probing-sensitive nucleotides. The reactivities obtained by probNORM are highly consistent with the structural models allowing the separation of single- and double-stranded nucleotides.

## The overview of probNORM method

<p align="center">
<img src="Fig_1_pipeline.png" width="800"/>
</p>
**A.** Counting RT polymerase stops from BAM files (marked with dashed line). 

**B.** The estimation of the distribution of log2 fold changes enables the calculation of the normalization factor (nf). 

**C.** In the background signal after correction, normalized control counts are an approximate value of counts in the treated sample. Only positions meeting the quality filtering criteria are used for further calculations. 

**D.** P-value calculation. 

&nbsp;&nbsp;&nbsp;&nbsp;*1)* The negative part of the log2 fold change distribution is extracted and mirrored on the right side; 

&nbsp;&nbsp;&nbsp;&nbsp;*2)* Next, to the newly formed distribution, probNORM is fitting the gaussian distribution - creating the background distribution (BD); 

&nbsp;&nbsp;&nbsp;&nbsp;*3)* probNORM is calculating the cumulative distribution function to estimate the p-value for each position of the transcript. 

**E.** Transcript reactivity profile. 

&nbsp;&nbsp;&nbsp;&nbsp;*1)* Normalized counts are used for reactivity calculation for each position by subtraction of normalized control signal from treated; 

&nbsp;&nbsp;&nbsp;&nbsp;*2)* The final profile is scaled with 2-8% normalization.