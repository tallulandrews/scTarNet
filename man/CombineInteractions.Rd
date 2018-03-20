\name{combineInteractions}
\alias{combineInteractions}
\title{Combine pdCor output}
\description{
	Combines output from calculateConditionalCors for multiple batches/normalizations.
}
\usage{
	conbineInteractions(list_of_correlation_tables)
}
\arguments{
  \item{list_of_correlation_tables}
}
\details{
	Counts the number of times TF-TF interactions were inferred across multiple batches or normalizations. And counts the number of times each TF-TF-Target trio was inferred across multiple batches.
}
\value{
	list containing two tables :
	pair2consists with columns "TF1", "TF2", "recurr"
	pair2targets with columns "TF1", "TF2", "Target", "recurr"

	"recurr" is the number of batches the relationship was obsered in.
}
\examples{
	set.seed(13)
	batchA_output <- data.frame(
		TF1=sample(c("A", "B", "C", "D"), 50, replace=TRUE),
		TF2=sample(c("A", "B", "C", "D"), 50, replace=TRUE),
		Target=sample(c("E", "F", "G", "E", "F"), 50, replace=TRUE)
		)
	batchA_output <- list(direct_int=unique(batchA_output))


	batchB_output <- data.frame(
		TF1=sample(c("A", "B", "C", "D"), 50, replace=TRUE),
		TF2=sample(c("A", "B", "C", "D"), 50, replace=TRUE),
		Target=sample(c("E", "F", "G", "E", "F"), 50, replace=TRUE)
		)
	batchB_output <- list(direct_int=unique(batchB_output))

	consensus <- combineInteractions(list(batchA_output, batchB_output));
}

