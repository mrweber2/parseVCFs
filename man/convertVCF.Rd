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
}
\examples{
## Make sample vcf
test_vcf = "##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSRR871957.bam
Chr01\t10084232\t.\tC\tT\t99.2499\t.\tDP=41;VDB=2.76634e-22;SGB=32.1809;MQ0F=0;AF1=1;AC1=8;DP4=0,0,0,40;MQ=23;FQ=-29.6018\tGT:PL\t1/1:0,0,0\t1/1:0,0,0\t1/1:23,3,0\t1/1:111,117,0
Chr01	10084236\t.\tC\tT\t99.2499\t.\tDP=41;VDB=2.76634e-22;SGB=32.1809;MQ0F=0;AF1=1;AC1=8;DP4=0,0,0,40;MQ=23;FQ=-29.6018\tGT:PL\t1/1:0,0,0\t1/1:0,0,0\t1/1:23,3,0\t1/1:111,117,0
Chr01\t10084276\t.\tG\tC\t99.2499\t.\tDP=41;VDB=7.12062e-27;SGB=33.124;MQ0F=0;AF1=1;AC1=8;DP4=0,0,0,41;MQ=23;FQ=-29.6018\tGT:PL\t1/1:0,0,0\t1/1:0,0,0\t1/1:23,3,0\t1/1:111,120,0"

write(test_vcf, file="sample_parseVCFs.vcf")

## Execute function
test <- convertVCF("sample_parseVCFs.vcf", format="csv")
test
}