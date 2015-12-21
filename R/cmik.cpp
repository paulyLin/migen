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

    // Don't need these to be vectors, keep for the moment
    // to check intermediate results.
    NumericVector s(n);
    NumericVector sN(n);
    NumericVector t(n);
    NumericVector tN(n);
    NumericVector dists(2);
    for (int i = 0; i < n; i++)
    {
        // N.b. C++ std:sort is faster than R
        //
        // R matrices are row-major. Need to check
        // if RCpp NumericMatrix is too.
        NumericVector xdi = xdiffs(i, _);
        NumericVector ydi = ydiffs(i, _);
        NumericVector d = distmat(i, _);

        // Sorting breaks connection between X & Y.
        // May need to fins another approach.
        std::sort(xdi.begin(), xdi.end());
        std::sort(ydi.begin(), ydi.end());
        std::sort(d.begin(), d.end());
 
        // get s(i)
        for (int j = 0; j < n; j++)
        {
            if (xdi[j] < bws[0]) 
            {
                s[i] += 1;
            } 
            else
            {
                break;
            }
        }
        sN[i] = xdi[s[i] + 1]; // wrong
        double xeD = d[s[i] + 1]; // wrong
        if (i == 0) dists[0] = xeD; //wrong

        // get t(i)
        for (int j = 0; j < n; j++)
        {
            if (ydi[j] < bws[1]) 
            {
                t[i] += 1;
            } 
            else
            {
                break;
            }
        }
        tN[i] = ydi[t[i] + 1]; // wrong
        double yeD = d[t[i] + 1]; // wrong
        if (i == 0) dists[1] = yeD; //wrong

    }

    List ret;
    ret["xdiffs"] = xdiffs;
    ret["ydiffs"] = ydiffs;
    ret["distmat"] = distmat;
    ret["s"] = s;
    ret["t"] = t;
    ret["sN"] = sN; // values rather than indices.
    ret["tN"] = tN;
    ret["dists"] = dists;
    return(ret);
} 
