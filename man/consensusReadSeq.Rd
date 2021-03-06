\name{consensusReadSeq}
\alias{consensusReadSeq}

\title{Consensus read sequence}
\description{Generate a consensus sequence from all reads in the same UMI group, based on their multiple sequence alignment.}

\usage{
consensusReadSeq(alignments, pseudo.count=1, min.coverage=0.6, 
    BPPARAM=SerialParam())
}

\arguments{
    \item{alignments}{A \linkS4class{DataFrame} with a \linkS4class{List} of alignment strings and (optionally) a List of quality scores, as produced by \code{\link{multiReadAlign}}.}
    \item{pseudo.count}{A numeric scalar specifying the pseudo-count that should be added when computing quality scores for each base.}
    \item{min.coverage}{A numeric scalar specifying the minimum proportion of reads with a base at a particular position in order to define a base position in the consensus.}
    \item{BPPARAM}{A BiocParallelParam object controlling how paralellization should be performed.}
}

\details{
This function generates a consensus sequence from a set of alignment strings, obtained from a multiple sequence alignment (MSA). 
Positions in the MSA are only retained in the consensus if the proportion of reads that have a base at that position exceeds \code{min.coverage}.
of these retained positions, the base in the consensus is chosen as the most abundant nucleotide in the MSA.

A Phred quality string is also calculated for each position in each read.
The probability of a base position being incorrect in the consensus is defined as 1 minus the (adjusted) proportion of reads that have the chosen base.
This uses the proportion as an estimate of the likelihood, and it applies an equal prior probability to each base.

For each position, the above proportion is calculated using only the set of reads that contain a base at that position.
To avoid error probabilities of zero when all reads contain the chosen base, \code{pseudo.count} is added to the total number of reads.
This effectively shrinks the probabilities towards a prior value of 0.25.
}

\value{
A \linkS4class{QualityScaledDNAStringSet} containing the consensus sequence for each UMI group, along with the associated Phred quality scores.
Each entry is named according to \code{rownames(alignments)}.
}

\author{
Florian Bieberich and Aaron Lun
}

\seealso{
    \code{\link{umiGroup}},
    \code{\link{multiReadAlign}}
}

\examples{
reads <- DNAStringSet(c("ACACTGGTTCAGGT", 
    "ACACGGTTCAGGT",
    "CGGACTGACACGGT",
    "CGGGCTGACACGGT"))
aln <- multiReadAlign(reads, c(1,1,2,2))
consensusReadSeq(aln)

qreads <- QualityScaledDNAStringSet(reads,
    PhredQuality(c("23849723948733", 
    "*(&^&23498342",
    ")(*!@!343928&^",
    "2($*&$*A323128")))
aln0 <- multiReadAlign(qreads, c(1,1,2,2))
consensusReadSeq(aln0)
}


