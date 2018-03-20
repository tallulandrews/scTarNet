\name{calculateConditionalCors}
\alias{calculateConditionalCors}
\title{Calculate gene-gene correlations with dCor}
\description{
	Calculates gene-gene correlations for a set of TFs to all genes in the genome using the dCor statistic from the "energy" package.
}
\usage{
	calculateConditionalCors(Mat, TFs, Dep, n.cores=1, threshold.interaction=0.01, bidirectional=TRUE, threshold.indirect=0.5, exclude.indirect=TRUE)
}
\arguments{
  \item{Mat}{gene by sample expression matrix after normalization and/or batch correction.}
  \item{TFs}{vector of gene names of the TFs of interest.}
  \item{Dep}{output from calculateTFstoTargets}
  \item{n.cores}{number of cores to run in parallel.}
  \item{threshold.interaction}{threshold for calling pair of correlations as a TF-TF interaction.}
  \item{bidirectional}{whether TF-TF interactions must be bidirection. (see Details)}
  \item{threshold.indirect}{threshold for calling a correlation as indirect.}
  \item{exclude.indirect}{whether to remove indirect interactions from "Adj"}
}
\details{
For each trio of two TFs (A, B) and a common target gene (T), calculates the
conditional distance correlations:

pdCor(A, T | B) and pdCor(B, T | A)

if either of these are approximately zero (pdCor < threshold.indirect*original dCor) then infers an indirect interaction:
pdCor(A, T | B) < T_indir * dCor(A, T) => A -> B -> T
pdCor(B, T | A) < T_indir * dCor(B, T) => B -> A -> T

Indirect interactions are optionally removed from the dCor adjacency matrix.

if the conditional correlations are more than threshold.interaction bigger than original, then infers a TF-TF interaction:
pdCor(A, T | B) > (1 + T_inter)*dCor(A, T) => A + B -> T
pdCor(B, T | A) > (1 + T_inter)*dCor(B, T) => A + B -> T

if `bidirectional` is TRUE then both of the above must be true to call a TF-TF interaction.
TF-TF interaction trios are saved to a dataframe and returned.

Runs in parallel as specified by n.cores.
}
\value{
A list of two dataframe:
	Dep : the input TF to Target dataframe with indirect interaction removed or as provided (as specified with exclude.indirect)
	Int : a dataframe of TF-TF interactions with columns: TF1, TF2, Target
}
\examples{
	Mat <- matrix(rnorm(1000), nrow=10);
	rownames(Mat) <- c("A", "B", "C", "D", "E",
			   "F", "G", "H", "I", "J")
	Dep <- calculateTFstoTargets(Mat, c("A", "B", "C", "D"), n.cores=1, mt_correction="fdr=2")
	out <- calculateConditionalCors(Mat, c("A", "B", "C", "D"), Dep, n.cores=1, threshold.indirect=0.2, threshold.interaction=0.1, exclude.indirect=FALSE)
}
\seealso{calculateTFstoTargets}
