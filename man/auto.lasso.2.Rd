\name{auto.lasso}
\alias{auto.lasso}
\title{
Auto Lasso Method
}
\description{
Auto lasso function returns estimate beta and estimate sigma by finding minimum risk (m) largest beta's in absolute values, and obtaining least square estimates with active minimum risk (m) variables.
}
\usage{
auto.lasso.2(x,y,ratio=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{x}{
Predicted variable
}
  \item{y}{
Response variable
}
  \item{ratio}{
Ratio is the proportion of the number of columns over that of the rows.

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
auto.lasso=function(x,y,ratio=0.5)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
