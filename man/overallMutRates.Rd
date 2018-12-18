\name{overallMutRates}
\alias{overallMutRates}
\title{
Compute the Overall Mutation Rates for Point Mutations
}
\description{
This function takes in a VCF file and computes the frequency of each specific single-point mutation out of all point mutations in the data set.
}
\usage{
overallMutRates(vcf)
}
\arguments{
  \item{vcf}{
Input raw VCF file. Standard tab-separated variant file.
}
}
\value{
\code{overallMutRates} returns a list where key names equal the allele change and values equal mutation frequency.
}
\examples{
## Make sample vcf
test_vcf = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSRR871957.bam
Chr01\t10084232\t.\tC\tT\t99.2499\t.\tDP=41\tGT:PL\t1/1:0,0,0\t1/1:0,0,0\t1/1:23,3,0\t1/1:111,117,0
Chr01	10084236\t.\tC\tT\t99.2499\t.\tDP=41\tGT:PL\t1/1:0,0,0\t1/1:0,0,0\t1/1:23,3,0\t1/1:111,117,0
Chr01\t10084276\t.\tG\tC\t99.2499\t.\tDP=41\tGT:PL\t1/1:0,0,0\t1/1:0,0,0\t1/1:23,3,0\t1/1:111,120,0"

write(test_vcf, file="sample_parseVCFs.vcf")

## Execute function
MutRates <- overallMutRates("sample_parseVCFs.vcf")
MutRates
}