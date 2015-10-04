\name{cmi}
\alias{cmi}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate mutual information between two continuous variables
%%  ~~function to do ... ~~
}
\description{This function calculates MI between two continuous variables.
    %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cmi(v1,v2,method="Kernel-Smooth",h,level=3L,k=3,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v1}{ The first continuous variable. Should be a vector.}
  \item{v2}{The second continuous variable. Should be a vector.}
  \item{method}{ The method used to calculate mutual informion. Either "Kernel-Smooth" (default) or "k-Neighbor".
%%     ~~Describe \code{cts} here~~
}
  \item{level}{ The number of levels used for plug-in bandwidth estimation in the Kernel-Smooth method (see
      the documentation for the KernSmooth package).
%%     ~~Describe \code{level} here~~
}
  \item{na.rm}{ Remove missing values if TRUE. This is required for the
      bandwidth calculation.
%%     ~~Describe \code{na.rm} here~~
}
  \item{h}{ A (double) vector of smoothing bandwidths, one for each variable. If
      missing this will be calculated using the dpik() function from the
      KernSmooth package.
%%     ~~Describe \code{h} here~~
}
  \item{k}{The number of neighbors used in the k-Neighbor method.}
  \item{\dots}{ Additional options passed to dpik() if necessary.
%%     ~~Describe \code{\dots} here~~
} 
}
\details{ It calculates mutual information between two continuous variables using either the Kernel-Smooth method (default) or k-Neighbor method. 
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns the mutual information.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...

%% \references{
%% %% ~put references to the literature/web site here ~
%% }
%% \author{
%% %%  ~~who you are~~
%% }
%% \note{
%% %%  ~~further notes~~
%% }
%% 
%% %% ~Make other sections like Warning with \section{Warning }{....} ~
%% 
%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
##################################################
# The USArrests dataset

# Kernel Smooth method
m1 <- cmi(USArrests[,1],USArrests[,2])

# k-neighbor method
m2 <- cmi(USArrests[,1],USArrests[,2],method="k-Neighbor")

# More examples
m3 <- cmi(state.x77[,1], state.x77[,2])

m4 <- cmi(state.x77[,1], state.x77[,2],method="k-Neighbor")}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
