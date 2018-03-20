\name{plotInteractions}
\alias{plotInteractions}
\title{Plot pdCor TF-TF interactions}
\description{
	Plot the TF-TF interactions calculated by calculateConditionalCors
}
\usage{
	plotInteractions(combined_interactions, recurr_threshold=max(combined_interactions$rucurr)/2, suppress.plot=FALSE, seed.col="dodgerblue", interaction.col="black")
}
\arguments{
  \item{combined_interactions}{Output from combineInteractions}
  \item{recurr_threshold}{number of each particular a correlation must be seen.}
  \item{suppress.plot}{whether to plot the network}
  \item{seed.col}{colour of seed TF nodes}
  \item{interaction.col}{colour of edges for TF-TF interactions}
}
\details{
	Filters TF-TF interactions based on the frequency it was observed across replicates. 
	Plots the resulting network using igraph with colour coded nodes and edges. 
}
\value{
	list of the edgelist of the network as a data.frame and the igraph object of the network.
}
\examples{
	set.seed(101)
	pair2targets <- data.frame(TF1=sample(c("A", "B", "C"), 20, replace=TRUE), TF2=sample(c("A", "B", "C"), 20, replace=TRUE), Targets=sample(c("D", "E", "F", "G", "H"), 20, replace=TRUE))
	pair2targets <- unique(pair2targets);
	pair2targets$recurr <- rpois(nrow(pair2targets), lambda=1)+1

	pair2consistent <- pair2targets[,1:2]
	pair2consistent <- unique(pair2consistent);
	pair2consistent$recurr <- rpois(nrow(pair2consistent), lambda=1)+1

	inteactions_out <- list(pair2targets=pair2targets, pair2consistent=pair2consistent);

	plotInteractions(inteactions_out, recurr_threshold=2, suppress.plot=TRUE)
}
\seealso{igraph.plotting}
