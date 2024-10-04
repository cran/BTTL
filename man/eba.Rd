\name{eba}
\alias{eba}
\docType{data}
\title{
Pairwise comparison data
}
\description{
This is a real dataset containing pairwise preferences on 13 different eba varieties across eight attributes.
}
\usage{data("eba")}
\format{
  A list with eight matrices. Each matrix consists of pairwise comparisons across one of the eight attributes. The first column contains the preferred object, whilst the second column contains the non-preferred object. In order of appearance, the matrices contain pairwise comparisons on the following attributes: consumer overall preference, color, odour, firmness, stretchability, taste, smoothness and mouldability.
}
\details{
Contains the data used in the Hermes et al. (2024) paper.
}
\source{
Data from the Hermes et al. (2024) paper is based on Olaosebikan et al. (2023).
}
\references{
1. Hermes, S., van Heerwaarden, J., and Behrouzi, P. (2024). Multi-Attribute Preferences: A Transfer Learning Approach. arXiv preprint, arXiv:2408.10558 \cr
2. Olaosebikan, O., A. Bello, K. De Sousa, R. Ndjouenkeu, M. Adesokan, E. Alamu, A. Agbona, J. Van Etten,
F. N. K´egah, D. Dufour, et al. (2023). Drivers of consumer acceptability of cassava gari-eba food products
across cultural and environmental settings using the triadic comparison of technologies approach (tricot).
Journal of the Science of Food and Agriculture. \cr
}
\examples{
data(eba)
}
\keyword{datasets}
