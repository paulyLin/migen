#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List cmikCpp(NumericVector bws, 
        NumericVector N, NumericVector X, NumericVector Y) 
{
    int n = N[0];

    // X differences
    NumericMatrix xdiffs(n);
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            if (j == i) continue; // Keep 0s on diagonal.
            
            xdiffs(i, j) = std::abs(X[i] - X[j]);
            xdiffs(j, i) = xdiffs(i, j); // Symmetry of abs().
        }
    }

    // Y differences
    NumericMatrix ydiffs(n);
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            if (j == i) continue; // Keep 0s on diagonal.
            
            ydiffs(i, j) = std::abs(Y[i] - Y[j]);
            ydiffs(j, i) = ydiffs(i, j); // Symmetry of abs().
        }
    }

    // Max norm distances
    NumericMatrix distmat(n);
    for (int i = 0; i < n; i++)
    {
        for (int j = i; j < n; j++)
        {
            if (j == i) continue; // Keep 0s on diagonal.
            
            distmat(i, j) = std::max(xdiffs(i, j), ydiffs(i, j));
            distmat(j, i) = distmat(i, j); // Symmetry.
        }
    }

    List ret;
    ret["xdiffs"] = xdiffs;
    ret["ydiffs"] = ydiffs;
    ret["distmat"] = distmat;
    return(ret);
} 
