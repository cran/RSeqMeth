\name{analyzeSequenom}
\alias{analyzeSequenom}
\title{Calculates methylation ratios and proportions}
\description{
The user is interactively asked to specify 
}
\usage{
analyzeSequenom(matchedName = NA, methName = NA, seqNames = NA, cull=TRUE, 
  writeOut = TRUE, quiet = FALSE, minMass = 1500, maxMass = 7000)
}
\arguments{
  \item{matchedName}{Name of the file to read the matched peak data from - if
  left to default then the user is asked to interactively choose the matched
  peaks table.}
  \item{methName}{If \code{cull} is \code{TRUE} then the name of the file to
  read methylation grid data from - if left to default then the user is asked to 
  interactively choose the methylation grid data. If \code{cull} is \code{FALSE} this value is
  ignored.}
  \item{seqNames}{A character vector of file names to read sequences from. The
  sequences are specified in the same order that their amplicons appear in the
  matched peaks table. If left to default then the user is asked to 
  interactively choose the sequence files used in analysis.}
  \item{cull}{logical: if \code{TRUE} (default) then methylation grid data will
  be read from \code{methName} and missing values from this grid will be used to
  cull methylation ratios that did not meet the EpiTyper Probability Threshold
  setting at the time the methylation grid was exported. if \code{FALSE} no culling
  will occur and the value of methName is ignored.}
  \item{writeOut}{logical: if \code{TRUE} (default) then the calculated 
  methylation values are writted to tab separated text files in the same 
  directory as \code{matchedName}. Two files are written:
  
  <matchedName> Weighted - contains for each informative CpGUnit (excluding
  overlaps) the methylation ratio calculated using the weighted formula.
  
  <matchedName> Proportion - contains for each informative CpGUnits (including
  overlaps) the proportion of each peak eg (0/2 CpGs methylated, 1/2 CpGs 
  methylated, 2/2 CpGs methylated)}
  \item{quiet}{logical: if \code{FALSE} no output is written to the screen}
  \item{minMass}{The minimal mass of a peak to be included in the analysis.}
  \item{maxMass}{The maximal mass of a peak to be included in the analysis.}  
}
  
\details{
Currently, only T cleavage data is interpreted.
Matched peak files containing C cleavage data have not been tested.
}
\value{
A list containing three elements:


Sequenom - Reproduces the calculations made by EpiTyper.

Weighted - Ratio for each CpGUnits (excluding overlaps) calculated using the
weighted formula.

Proportion - Proportion for each peak in each informative CpGUnit (including 
overlaps.


Each of these is a list which contains an element for each amplicon in the input
matched peaks table. Each of these amplicon elements contains two elements:


Name - The name of the amplicon as specified in the matched peaks table.

Ratios - A data frame containing for each sample in the matched peaks table the
methylation ratio or proportion as given above.
}
\references{Coolen MW, Statham AL, Gardiner-Garden M, Clark SJ. (2007) Nucleic Acids Research 35(18):e119 Genomic profiling of CpG methylation and allelic specificity using quantitative high-throughput mass spectrometry: critical evaluation and improvements.}
\author{Aaron Statham}
\seealso{\code{\link{ampliconReport}}}
\examples{
## run completely interactively, everything obtained from the user
#meth <- analyzeSequenom()

## specified matched peaks table, do not cull, sequences will be specified
## interactively by the user
#meth <- analyzeSequenom(matchedName = "Plate1 matched.txt", cull=FALSE)

## runs completely without user intervention - output is not written to files
#meth <- analyzeSequenom(matchedName = "Plate1 matched.txt", cull=FALSE, seqNames 
#= c("EN1.txt", "EN1.txt", "SCTR.txt", "SCTR.txt"), quiet=TRUE, writeOut=FALSE)

}
\keyword{manip}
