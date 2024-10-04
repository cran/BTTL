\name{bttl}
\alias{bttl}

\title{Bradley-Terry Transfer Learning}

\description{
Contains the main function of this package that is used to estimate the parameter of interest \eqn{\alpha}. The inner workings of the function are described in Hermes et al., (2024).
}

\usage{
bttl(y, M, c_s)
}

\arguments{
  \item{y}{A list consisting of \eqn{S} matrices with each matrix containing pairwise comparisons across \eqn{M} objects. The first matrix in the list should contain the pairwise comparisons on the primary attribute. The other matrices should contain pairwise comparisons on the secondary attributes. The first column of each matrix should contain the preferred object (winner) and the second column should contain the non-preferred object (loser).}
  \item{M}{Total number of objects encountered across the pairwise comparisons.}
  \item{c_s}{Constant determining the strictness of the inclusion of attributes in the informative set. Smaller values lead to smaller \eqn{\mathcal{S}}. We typically set this parameter equal to 1.}
}

\value{
\item{alpha_est}{A vector of length \eqn{M} that contains the parameter estimates \eqn{\hat{\alpha}} for all objects in the data on the primary attribute.}
\item{S_hat}{A vector containing the indices of attributes that are included in \eqn{\hat{\mathcal{S}}}.}
}

\references{
1. Hermes, S., van Heerwaarden, J., and Behrouzi, P. (2024). Multi-Attribute Preferences: A Transfer Learning Approach. arXiv preprint, arXiv:2408.10558 \cr
}

\author{
  Sjoerd Hermes\cr
  Maintainer: Sjoerd Hermes \email{sjoerd.hermes@wur.nl}
}

\examples{
\donttest{

# we first obtain the pairwise comparison data
y <- data(eba)

# now we can fit our model
result <- bttl(y, 13, 1)
}
}

