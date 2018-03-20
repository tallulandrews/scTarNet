\name{plotInteractionsWithTargets}
\alias{plotInteractionsWithTargets}
\title{Plot pdCor TF-TF interactions with their targets}
\description{
	Plot the gene-gene correlations inferred by calculateTFstoTargets, which are targets of the TF-TF interactions inferred by calculateConditionalCors.
}
\usage{
	plotInteractionsWithTargets(combined_dependencies, combined_interactions, recurr_threshold=max(combined_interactions$rucurr)/2, suppress.plot=FALSE, seed.col="dodgerblue", target.col="goldenrod1", interaction.col="black", edge.col=c("blue", "grey45", "red"))
}
\arguments{
  \item{combined_dependencies}{Output from combineDependencies}
  \item{combined_interactions}{Output from combineInteractions}
  \item{recurr_threshold}{number of each particular a correlation must be seen.}
  \item{suppress.plot}{whether to plot the network}
  \item{seed.col}{colour of seed TF nodes}
  \item{target.col}{colour of target nodes}
  \item{interaction.col}{colour of edges for TF-TF interactions}
  \item{edge.col}{colour of edges for negative, other, and positive correlations}
}
\details{
	Filters TF-TF interactions and TF-Target correlations based on the frequency it was observed across replicates. 
	Plots the resulting network using igraph with colour coded nodes and edges. 
}
\value{
	list of a dataframe of the edgelist of the network as a data.frame and the igraph object of the network.
}
\examples{
	set.seed(101)
	dependencies_out <- data.frame(Gene=rep(c("A", "B", "C"), times=3), Targets=rep(c("D", "E", "F"), each=3))
        dependencies_out$recurr <- rpois(nrow(dependencies_out), lambda=1)+1
        dependencies_out$direction <- sample(c(-1, 0, 1), nrow(dependencies_out), replace=TRUE)

	pair2targets <- data.frame(TF1=c("A", "A", "B", "B"), TF2=sample(c("B", "C", "C", "C"), 20, replace=TRUE), Targets=c("D", "E", "E", "F")) # Nope this will still create novel TF-Target pairs
	pair2targets$recurr <- rpois(nrow(pair2targets), lambda=1)+1

	pair2consistent <- data.frame(TF1=c("A","A","B"), TF2=c("B", "C", "C"))
	pair2consistent$recurr <- rpois(nrow(pair2consistent), lambda=1)+1

	inteactions_out <- list(pair2targets=pair2targets, pair2consistent=pair2consistent);

	plotInteractionsWithTargets(dependencies_out, interactions_out, recurr_threshold=2, suppress.plot=TRUE)
}
\seealso{igraph.plotting}
