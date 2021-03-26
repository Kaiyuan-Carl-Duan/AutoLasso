\name{model.fitting}
\alias{model.fitting}
\title{
Simulation Program
}
\description{
The function returns a list of simulations that summary the biases, TPR, FPR, FDR, selected index, sigma2 bias of auto lasso method and lasso1 method.
}
\usage{
model.fitting(x,N.sim=100,p.true=5,sigma.true,beta.true)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Predicted variable
}
  \item{N.sim}{
Times of simulation
}
  \item{p.true}{
Number of active variables
}
  \item{sigma.true}{
True value of sigma
}
  \item{beta.true}{
True value of beta
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
model.fitting=function(x=x,N.sim=100,p.true=5,sigma.true=sigma.true,beta.true=beta.true)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
