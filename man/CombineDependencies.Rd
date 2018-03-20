\name{combineDependencies}
\alias{combineDependencies}
\title{Combine dCor output}
\description{
	Combines output from calculateTFstoTargets for multiple batches/normalizations.
}
\usage{
	conbineDependencies(list_of_correlation_tables)
}
\arguments{
  \item{list_of_correlation_tables}
}
\details{
	Counts the number of times TF-Target inferred across multiple batches or normalizations. And infers the consensus direction of the relationship as follows:
	if always "other" (0) then "other" (0)
	if always "other" (0) or "positive" (1) then positive (1)
	if always "other" (0) or "negative" (1) then negative (1)
	if mix of "positive" and "negative" then "other" (0)
}
\value{
	A table with columns "Gene" (TF), "Target", "recurr" (number of batches), and "direction"
}
\examples{
	set.seed(13)
	batchA_output <- data.frame(
		Gene=sample(c("A", "B", "C", "D"), 50, replace=TRUE),
		Target=sample(c("E", "F", "G", "E", "F"), 50, replace=TRUE)
		)
	batchA_output <- unique(batchA_output)
	batchA_output$direction <- sample(c(1,0,-1), nrow(batchA_output), replace=TRUE)


	batchB_output <- data.frame(
		Gene=sample(c("A", "B", "C", "D"), 50, replace=TRUE),
		Target=sample(c("E", "F", "G", "E", "F"), 50, replace=TRUE)
		)
	batchB_output <- unique(batchB_output)
	batchB_output$direction <- sample(c(1,0,-1), nrow(batchB_output), replace=TRUE)

	consensus <- combineDependencies(list(batchA_output, batchB_output));
}

