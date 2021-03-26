\name{sim.study}
\alias{sim.study}
\title{
Simulation Study
}
\description{
This program returns a graphical summary that contains Biases of all coefficients, Proportion of correct selection, Biases of active variables, and Biases of zero coefficient; a numerical summary of percentage of times each variables selected, TPR&FPR of different methods, Distribution of number of selected variables, Comparison between methods of number of active variables found, Biases of active variables, Relative Bias of sigma2, summary of overall miss classification, Real time used (subtract start time from completed) of Lasso Method and Auto Lasso Method directly.
}
\usage{
sim.study(n,N.sim,ratio,p.true,rho,cor.type,snr)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of rows of generated data
}
  \item{N.sim}{
Times of simulation
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
Different type of Sigma: Auto regressive, Compound symmetry
}
  \item{snr}{
Signal to noise ratio
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
sim.study(n=100,N.sim=100,ratio=0.5,p.true=5,rho=0,cor.type=1,snr=10)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
