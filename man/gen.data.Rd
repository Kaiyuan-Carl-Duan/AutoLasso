\name{gen.data}
\alias{gen.data}
\title{
Data Generation
}
\description{
Generate a random n by n*ratio matrix with the scale of beta in (-1,1). The number of active variables is p.true. User can input the correlation of data rho and signal to noise ratio. Also change the two different correlation generate methods.
}
\usage{
gen.data(n=100,ratio=0.5,p.true=5,rho=0,cor.type=1,snr=10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of rows of generated data
}
  \item{ratio}{
Ratio is the proportion of the number of columns over that of the rows.
}
  \item{p.true}{
Number of active variables
}
  \item{rho}{
Correlation generate methods
}
  \item{cor.type}{
Different type of Sigma: Auto regressive(1), Compound symmetry(2)
}
  \item{snr}{
Signal to noise ratio
}
}
\details{
n=100,ratio=0.5,p.true=5,rho=0,cor.type=1,snr=10
}
\value{
A list contains three matrices: response variables, predicted variables, beta true and one scalar sigma true.
}
\references{
A Fast Automatic Variable Selection Method for LinearRegression Model
}
\author{
Sujit K. Ghosh, Kaiyuan Duan, Guangjie Yu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
mydata=gen.data(n=100,ratio=0.1,rho=0.25)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
