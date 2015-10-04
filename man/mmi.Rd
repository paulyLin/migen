\name{mmi}
\alias{mmi}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Calculate mixed-pair mutual information between a continuous variable and a discrete variable.
%%  ~~function to do ... ~~
}
\description{This function calculates mutual information between a continuous variable
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
    \item{disc}{ Discrete variable. Should be a vector.}

%%
  \item{method}{ The method used to calculate mutual information. Either "Kernel-Smooth" (default) or "k-Neighbor".
}

%%
  \item{level}{ The number of levels used for plug-in bandwidth estimation in the Kernel-Smooth method (see
      the documentation for the KernSmooth package.)
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
 \item{ki}{A vector of length n, where n is the number of categories 
 in the discrete variable. Each entry of ki is the number of neighbors 
 used in the k-Neighbor method.}
  \item{\dots}{ Additional options passed to dpik() if necessary.
%%     ~~Describe \code{\dots} here~~
} 
}
\details{It calculates mutual information between a 
continuous variable and a discrete variable using either the Kernel-Smooth method (default) or k-Neighbor method.
%%  ~~ If necessary, more details than the description above ~~
}
\value{ Returns the mutual information.
   

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
cts <- c(rnorm(50,-5,1), rnorm(50,0,1), rnorm(50,5,1))
disc <- c(rep("A",50), rep("B",50), rep("C",50))

# Use Kernel-Smooth method.
m1 <- mmi(cts, disc)

# Use k-Neighbor method.
m2 <- mmi(cts, disc, method="k-Neighbor")

# More examples. 

m3 <- mmi(state.x77[,1],state.region)

m4 <- mmi(state.x77[,1],state.region,method="k-Neighbor")}
