% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{str2vec}
\alias{str2vec}
\title{array-ify a string with a possible leading path/URL/whatever}
\usage{
str2vec(x, y = "_")
}
\arguments{
\item{x}{a string}

\item{y}{a split character}
}
\value{
the array produced by splitting basename(x) on y
}
\description{
array-ify a string with a possible leading path/URL/whatever
}

