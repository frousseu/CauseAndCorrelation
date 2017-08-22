\name{dsep.test}
\alias{dsep.test}

\title{
dsep.test
}
\description{
This function impliments the dsep test described in Shipley(2000), given a DAG and a covariance matrix of the variables in the DAG.  It is identical to the shipley.test function in the ggm library, except that, as a default, it only outputs the null probability of the data given the DAG.  This function is meant to be called within other functions.  You should use the shipley.test2 function if you want to get complete output.
}
\usage{
dsep.test(amat, S, n, only.null = F)
}

\arguments{
  \item{amat}{a binary matrix holding the directed acyclic graph; usually produced via the DAG function in the ggm library}
  \item{s}{the covariance matrix of the variables in the DAG}
  \item{n}{the number of independent multivariate observations involved in the covariance matrix}
  \item{only.null}{a logical argument.  If TRUE, only the null probability of the data given the DAG is output; otherwise the output mimics the shipley2.test function}
}
\details{

}
\value{
the null probability
}
\references{
Shipley, B. (2000). A new inferential test for path models based on directed acyclic graphs. Structural Equation Modeling, 7(2), 206–218.
}
\author{
  Bill Shipley
}
\note{

}



\seealso{

}
\examples{
dsep.test(amat=test,S=cov(Fig3.1.data),n=100,only.null=T)
[1] 0.8438591

> dsep.test(amat=test,S=cov(Fig3.1.data),n=100,only.null=F)
Individual d-sep claims in basis set 
A _||_ D |{ B } r= 0.051  p= 0.618 
A _||_ C |{ B } r= -0.118  p= 0.245 
A _||_ E |{ D C } r= 0.066  p= 0.521 
B _||_ E |{ A D C } r= 0.019  p= 0.855 
D _||_ C |{ B } r= -0.015  p= 0.88 
$ctest
[1] 5.648765

$df
[1] 10
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }
