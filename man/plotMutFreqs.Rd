\name{plotMutFreqs}
\alias{plotMutFreqs}
\title{
Plot Point Mutation Frequencies
}
\description{
This function takes in the output list of point mutation frequencies from overallMutRates and plots a histogram of frequencies.
}
\usage{
plotMutFreqs(mutList)
}
\arguments{
  \item{mutList}{
Input list of mutational frequencies from 0 to 1.00. List names are possible point mutation allele changes.}
}

\value{
\code{multList} returns a histogram of mutation frequencies for each point mutation.
}
\examples{
mutList <- as.list(c(0.05, 0.2, 0.15, 0.05, 0.1, 0.15, 0.05, 0.05, 0.05, 0.05, 0.01, 0.09))
names(mutList) <- c("A2C", "A2G", "A2T", "C2A", "C2G", "C2T", "G2A", "G2C", "G2T", "T2A", "T2C", "T2G")
plotMutFreqs(mutList)
}