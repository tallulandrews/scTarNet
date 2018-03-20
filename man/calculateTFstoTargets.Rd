\name{calculateTFstoTargets}
\alias{calculateTFstoTargets}
\title{Calculate gene-gene correlations with dCor}
\description{
	Calculates gene-gene correlations for a set of TFs to all genes in the genome using the dCor statistic from the "energy" package.
}
\usage{
	calculateTFstoTargets(Mat, TFs, n.cores=1, mt_correction="bon=0.05")
}
\arguments{
  \item{Mat}{gene by sample expression matrix after normalization and/or batch correction.}
  \item{TFs}{vector of gene names of the TFs of interest.}
  \item{n.cores}{number of cores to run in parallel.}
  \item{mt_correction}{string of the multiple-testing method to apply and the significance threshold, separated by an equals sign.}
}
\details{
Calculates all pairwise correlations between a provided set of TFs using the 
distance correlation (dCor) statistic from the "energy" package. Runs in 
parallel on the specified number of cores. Applies the specified multiple 
testing correction. 

The direction of the correlation is inferred using a t-test between the the 
values of the target for the top 25\% and bottom 25\% of the expression of 
the TF:
-1 = negative
0 = unclear/other
1 = positive
}
\value{
	Data frame with columns : Gene, Target, pval (dcor p.value), strength (dcor value), direction (-1, 0, 1)
}
\examples{
	Mat <- matrix(rnorm(1000), nrow=10);
	rownames(Mat) <- c("A", "B", "C", "D", "E",
			   "F", "G", "H", "I", "J")
	out <- calculateTFstoTargets(Mat, c("A", "B", "C", "D"), n.cores=1, mt_correction="fdr=0.2")
}
\seealso{getDirection}
