\name{getDirection}
\alias{getDirection}
\title{get direction for dCor}
\description{
	Uses a t-test to infer the direction of the distance correlation.
}
\usage{
	getDirection(...)
}
\arguments{
  \item{...}{arguments passed to t.test()}
}
\details{
	Silently performs a t-test. Returns the sign of the test statistic if the result is significant at p-value < 0.05). If not significant or any error occurs returns a zero.
}
\value{
	1 (second group > first group), 0 (not significant), or -1 (second group < first group)
}
\examples{
	a <- rnorm(50)
	b <- rnorm(50, mean=10)
	dir <- getDirection(a, b) # 1
	dir <- getDirection(b, a) # -1
	dir <- getDirection(a, a) # 0
}

