\name{ampliconReport}
\alias{ampliconReport}
\title{Generate report files}
\description{
Specified sequence files are read and reports generated.

Sequence files must be plaintext - no header.  Whitespace is allowed.

Four files are generated:

<SEQUENCE NAME> T/C Report.csv - analysis of the CpG containing fragments

<SEQUENCE NAME> T/C Spectra.pdf - illustration of the predicted spectra

}
\usage{
ampliconReport(fnames = NA, minMass = 1500, maxMass = 7000)
}
\arguments{
  \item{fnames}{A character vector of file names of sequences to be analyzed. If
   left to default, files are chosen interactively.}
  \item{minMass}{The minimal mass of a peak to be included in the analysis.}
  \item{maxMass}{The maximal mass of a peak to be included in the analysis.}     
}
\details{
Fragmentation is analyzed for both the C and T Cleavage reactions.
}
\value{
No value returned - reports generated are written to the directory the
sequence files are read from.
}
\references{Coolen MW, Statham AL, Gardiner-Garden M, Clark SJ. (2007) Nucleic Acids Research 35(18):e119 Genomic profiling of CpG methylation and allelic specificity using quantitative high-throughput mass spectrometry: critical evaluation and improvements.}
\author{Aaron Statham}
\seealso{\code{\link{analyzeSequenom}}}
\examples{
## choose sequence files interactively
#ampliconReport()

## generate report files from EN1.txt and SCTR.txt in the current working directory
#ampliconReport(c("EN1.txt", "SCTR.txt"))
}
\keyword{manip}
