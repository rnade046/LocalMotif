# LocalMotif

## Overview
#### Rationale
Protein localization is regulated through various mechanisms. Among these, sequence motifs in mRNA 3’ Untranslated Regions (3’ UTRs) are known to be involved in the transport of transcripts to their subcellular compartment, where protein translation will then often occur. However, most regulatory sequence elements remain uncharacterized. Proximity labelling coupled with tandem mass spectrometry has generated comprehensive protein co-localization networks that can address this knowledge gap. Indeed, proteins that are densely connected in such networks are likely to be in similar cell compartments. If such sets of proteins share a common mRNA sequence motif, this motif could be directly or indirectly related to their localization.

#### Description
LocalMotif is a tool for motif discovery in the context of protein co-localization networks (e.g. Human Cell Map, Go, Knight *et al.* Nature 2021). This algorithm measures the clustering of proteins with a shared 3'UTR sequence motif within the network. The significance of this clustering measure is then estimated from a normal distribution approximated from a Monte Carlo Sampling Distribution. To correct for multiple hypothesis testing, we assess the false discovery rate at various thresholds against a null model of randomized sequences.

LocalMotif is a series of Java applications that can be run from the command line. However, given the time complexity, a server with job parallelization is recommended. The final step of this pipeline is an R script provided as a Jupyter Notebook for easy user interfacing.

## Prerequisites
* Java Version 8+
* Required library: The Apache Commons Mathematics Library\
(`commons-math3-3.6.1.jar` included in project)
<br><br>
* Jupiter notebook (latest)
* R

## Compiling instructions
Each step of the pipeline is required to be compile prior to execution. Unless otherwise specified, seperate project folder (e.g. `enumerateMotifs/`) should be compiled as such:
```bash
javac enumerateMotifs/*.java
```

## Input files
1. Sequences Fasta file of 3'UTR sequences corresponding to proteins in network. Additional sequences can be included but will be ignored.

2. Protein Ids - .tsv file mapping protein names (as found in the network file) to their RefSeq identifyers (as found in the .fasta file).

   Protein Name | RefSeqIds
   --- | ---
   Protein1 | XMXXX
   Protein2 | XMXXX\|XMXXXX

3. Network - user provided .tsv file of a protein-protein interaction network. 

   Protein1 | Protein2 | Score
   --- | --- | ---
   ABC | ACC | 0.75

4. Parameters file - supplied *params.txt* file. It is passed to the program as command line argument. It is important to specify the working directory, file paths and the proper parameters. A detailed explanation of these parameters can be found here. 

# Pipeline
![graphical representation of pipeline](http://url/to/img.png) 

## 1. Motif Enumeration
The motif enumeration step identifies motifs from the input sequences and generates an annotation file that maps the proteins associated to each motif.

#### Input files

1. Sequences (e.g. `3utrSequences.fasta`): user-provided 3'UTR sequences (.fasta) corresponding to proteins in network
2. Protein Ids (e.g. `proteinIds.tsv`): user-provided (.tsv) mapping protein names (as found in the network file)
3. Precompiled motifs (`motifs.tsv`): this file is included in the default project resources and lists all motifs.

#### Implementation
```bash
java enumerateMotifs motifs.tsv 3utrSequences.fasta proteinIds.tsv
```

#### Output
* Annotation file (.tsv)

## 2. Randomized sequences
This step enables the shuffling of sequences within non-overlapping windows. Randomized sequences are used to control for false discoveries.

#### Input files
1. Sequences (e.g. `3utrSequences.fasta`): user-provided 3'UTR sequences (.fasta) corresponding to proteins in network
2. Protein Ids (e.g. `proteinIds.tsv`): user-provided (.tsv) mapping protein names (as found in the network file)

#### Implementation
```bash
java randomizeSeq.RandomMain 3utrSequences.fasta proteinIds.tsv
```

#### Output
* Randomized sequences in a Fasta file (e.g.`3utrSequences_winShuffled.fasta`)

## 3. Network enrichment analysis
#### Input files
* Parameters file (`params.txt`), supplied parameters file should be updated with file paths for: 

    1. Protein Ids (e.g. `proteinIds.tsv`): user-provided (.tsv) mapping protein names (as found in the network file)
    2. Sequences (e.g. `3utrSequences.fasta`): user-provided 3'UTR sequences (.fasta) corresponding to proteins in network
    3. Annotation file (.tsv)
    4. Network - user provided .tsv file of a protein-protein interaction network. 
    
#### Compiling
```bash
javac -cp commons-math3-3.6.1.jar localEnrich/*.java
```

#### Implementation
Note : given the size of motifs to test and time complexity of this algorithm, it is recommended to execute this part of the program on a server with job parallelization capabilities.  

```bash
java localEnrich.Main params.txt
```
## 4. FDR estimation 
#### Input files
1. FWD p-values - output from local enrichment
2. Null p-values - output from local enrichment

#### Implementation
```bash
java estimateFDRs fwdPvalues nullPvalues
```
# 5. Summarize motifs
#### Input files
1. Motif details - output from local enrichment
2. p-value - value provided by user based on estimated FDR threshold 

#### Implementation
```bash
java summarizeMotifs pValue
```
## 6. Motif families
#### Input files
1. Significant motifs - 
2. Annotation subset - 

#### Implementation
R notebook - requires user input

