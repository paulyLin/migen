\name{mmi}
\alias{mmi}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate mixed-pair mutual information between a continuous variable and a discrete variable.
%%  ~~function to do ... ~~
}
\description{This function calculates MI between a continuous variable
    and a discrete variable.
    %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mmi(cts, disc, method="Kernel-Smooth",
                level=3L, na.rm=FALSE, h, bc=TRUE,
                ki=floor(table(factor(disc))/2),...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cts}{ Continuous variable. Should be a vector.
%%     ~~Describe \code{cts} here~~
}
    \item{disc}{ Discrete variable. Should be a vector.)}
%%
  \item{level}{ The number of levels used for plug-in bandwidth estimation (see
      the documentation for the KernSmooth package.)
}
%%
  \item{method}{ The method used to calculate mutual information. Either "Kernel-Smooth" or "k-Neighbor".
}
%%     ~~Describe \code{level} here~~
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
}
\details{mminjk() and mminjk.pw() return just the MI values without performing the
    jackknife. mmi.pw() and mminjk.pw() only require one bandwidth for the
    continuous variable. The number of processor cores used can be changed by
    setting the environment variable "OMP_NUM_THREADS" \emph{before} starting R.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns a list of 3 matrices each of size ncol(cts) by
    ncol(disc). Each row index represents a continuous variable and each column
    index a discrete variable. 
    \item{mi}{The raw MI estimates.}
   

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
# A dataset with discrete and continuous variables

cts <- c(rnorm(50,-5,1), rnorm(50,0,1), rnorm(50,5,1))
disc <- c(rep("A",50), rep("B",50), rep("C",50))

# Use Kernel-Smooth method.
m1 <- mmi(cts, disc)

# Use k-Neighbor method.
m2 <- mmi(cts, disc, method="k-Neighbor")

# More examples. 

m3 <- mmi(state.x77[,1],state.region)

m4 <- mmi(state.x77[,1],state.region,method="k-Neighbor")}
