\name{plotDependencies}
\alias{plotDependencies}
\title{Plot dCor gene-gene correlations}
\description{
	Plot the gene-gene correlations inferred by calculateTFstoTargets.
}
\usage{
	plotDependencies(combined_dependencies, recurr_threshold=max(combined_dependencies$recurr)/2, suppress.plot=FALSE, seed.col="dodgerblue", target.col="goldenrod1", edge.col=c("blue", "black", "red"), rnd.seed=2948)
}
\arguments{
  \item{combined_dependencies}{Output from combineDependencies}
  \item{recurr_threshold}{number of each particular a correlation must be seen.}
  \item{suppress.plot}{whether to plot the network}
  \item{seed.col}{colour of seed TF nodes}
  \item{target.col}{colour of target nodes}
  \item{edge.col}{colour of edges for negative, other, and positive correlations}
  \item{rnd.seed}{random number generator seed, to ensure reproducible plots.}
}
\details{
	Filters TF-Target relationships based on the frequency it was observed across replicates. 
	Plots the resulting network using igraph with colour coded nodes and edges. 
}
\value{
	list of the edgelist of the network as a data.frame and the igraph object of the network.
}
\examples{
	set.seed(101)
	dependencies_out <- data.frame(Gene=sample(c("A", "B", "C"), 20, replace=TRUE), Targets=sample(c("D", "E", "F", "G", "H"), 20, replace=TRUE))
	dependencies_out <- unique(dependencies_out);
	dependencies_out$recurr <- rpois(nrow(dependencies_out), lambda=1)+1
	dependencies_out$direction <- sample(c(-1, 0, 1), nrow(dependencies_out), replace=TRUE)

	plotDependencies(dependencies_out, recurr_threshold=2, suppress.plot=TRUE)
}
\seealso{igraph.plotting}
