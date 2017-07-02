\name{factor_counts}
\alias{factor_counts}
\title{Fequency of factor levels}
\description{
	Counts the frequency of each factor level.
}
\usage{
	factor_counts(vec)
}
\arguments{
  \item{vec}{vector of discrete values to be counted.}
}
\details{
	Efficiently count the frequency of each level of a factor. The provided vector is re-factored.
}
\value{
	Named vector of counts for each level of the factor.
}
\examples{
	factor_counts(c("A", "B", "B", "C","C","C"))
}

