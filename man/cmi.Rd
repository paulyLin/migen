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
function(v1,v2,method="Kernel-Smooth",h,level=3L,k=3,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cts}{ The data matrix. Each row is an observation and each column is a
      variable of interest. Should be numerical data.
%%     ~~Describe \code{cts} here~~
}
  \item{level}{ The number of levels used for plug-in bandwidth estimation (see
      the documentation for the KernSmooth package.)
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
  \item{\dots}{ Additional options passed to dpik() if necessary.
%%     ~~Describe \code{\dots} here~~
} 
    \item{v1}{ A vector for the pairwise version}
    \item{v2}{ A vector for the pairwise version}
}
\details{ The results of cmi() are in many ways similar to a correlation matrix,
    with each row and column index corresponding to a given variable.
    cminjk() and cminjk.pw() just returns the MI values without performing the
    jackknife. cmi.pw() and cminjk.pw() each only require two bandwidths, one
    for each variable. The number of processor cores used can be changed by
    setting the environment variable "OMP_NUM_THREADS" \emph{before} starting R.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns a list of 3 matrices each of size ncol(cts) by
    ncol(cts)
    \item{mi}{ The raw MI estimates.}
    \item{bcmi}{Jackknife bias corrected MI estimates (BCMI). These are each MI value
    minus the corresponding jackknife estimate of bias.}
    \item{zvalues}{Z-scores for each hypothesis that the corresponding 
    BCMI value is zero. These have poor statistical properties but can be useful
    as a rough measure of the strength of association.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
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
