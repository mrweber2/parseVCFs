# parseVCFs
This is the repository for the development of an R package created for a final project in CPSC 499. The purpose is to parse and process VCF files to produce useful output for downstream analysis.

## convertVCF
This function takes in a raw VCF file and converts to either bed or CSV format. The user specifies format as an option.

usage:
```
convert(vcf, format="bed")
```
Returns a vector in bed or CSV format.

## overallMutRates
This function takes in a VCF file and computes the frequency of each specific single-point mutation out of all point mutations in the data set.

Possible Mutations:
```
A->C
A->G
A->T
C->A
C->G
C->T
G->A
G->C
G->T
T->A
T->C
T->G
```

usage:
```
overallMutRates(vcf)
```
Returns a list where key names equal the allele change and values equal mutation frequency.

## plotMutFreqs
This function takes in the output list of point mutation frequencies from overallMutRates and plots a histogram of frequencies.

usage:
```
plotMutFreqs(mutList)
```