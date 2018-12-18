\name{convertVCF}
\alias{convertVCF}
\title{
Convert VCF Format
}
\description{
This function takes in a raw VCF file and converts to either bed or CSV format. The user specifies format as an option.
}
\usage{
convertVCF(vcf, format = "bed")
}
\arguments{
  \item{vcf}{
Input raw VCF file. Standard tab-separated variant file.
}
  \item{format}{
Format to convert to. Choices are "bed|BED" or "csv|CSV".
}
}
\value{
\code{convertVCF} returns a vector in either bed or CSV format.
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
test <- convertVCF("sample_parseVCFs.vcf", format="csv")
test
}