# LESMoNlocal

## Overview


## Dependencies
* Java Version 8+
* Required library: The Apache Commons Mathematics Library (commons-math3-3.6.1.jar)

## Implementations
![alt text](http://url/to/img.png) 

### Motif enumeration
#### Input files
1. Alphabet - user created .tsv file containing all characters to consider in the analysis. Example alphabet file is provided. 

   Character | Nucleotides
   --- | ---
   Y | C,T
   X | A,U,C,G
 
2. Sequences - user provided .fasta file containing 3'UTR sequences corresponding to proteins in network. Additional sequences can be included but will be ignored. Can be obtained from the UCSD genome browser. 
3. Protein Ids - .tsv file mapping protein names (as found in the network file) to their RefSeq identifyers (as found in the .fasta file).

   Protein Name | RefSeqIds
   --- | ---
   Protein1 | XMXXX
   Protein2 | XMXXX\|XMXXXX

#### Implementation
```bash
java enumerateMotifs alphabet.tsv 3utrSequences.fasta proteinIds.tsv outputWorkingDirectory
```
### Randomized sequences
#### Input files
1. Sequences - .fasta file containing 3'UTR sequences corresponding to proteins in network. Additional sequences can be included but will be ignored. Can be obtained from the UCSD genome browser. 
2. Protein Ids - .tsv file mapping protein names (as found in the network file) to their RefSeq identifyers (as found in the .fasta file).

   Protein Name | RefSeqIds
   --- | ---
   Protein1 | XMXXX
   Protein2 | XMXXX\|XMXXXX

#### Implementation
```bash
java randomizeSequences 3utrSequences.fasta
```
### Local enrichment
#### Input files
1. Network - user provided .tsv file of a protein-protein interaction network. 

   Protein1 | Protein2 | Score
   --- | --- | ---
   ABC | ACC | 0.75
   
2. Annotation - generated at the motif enumeration step

#### Implementation
```bash
java motifEnrichment network annotation file
```
### FDR estimation 
#### Input files
1. FWD p-values - output from local enrichment
2. Null p-values - output from local enrichment

#### Implementation
```bash
java estimateFDRs fwdPvalues nullPvalues
```
### Summarize motifs
#### Input files
1. Motif details - output from local enrichment
2. p-value - value provided by user based on estimated FDR threshold 

#### Implementation
```bash
java summarizeMotifs pValue
```
### Motif families
#### Input files
1. Significant motifs - 
2. Annotation subset - 

#### Implementation
R notebook - requires user input

