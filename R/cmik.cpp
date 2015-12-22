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
    int sN_index = 0;
    double sN_value = 0.0;
    bool sN_candidate_found = false;
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
        // std::sort(xdi.begin(), xdi.end());
        // std::sort(ydi.begin(), ydi.end());
        // std::sort(d.begin(), d.end());
 
        // get s(i)
        sN_index = 0;
        sN_value = 0;
        sN_candidate_found = false;
        for (int j = 0; j < n; j++)
        {
            if (xdi[j] < bws[0])
            {  
                s[i] += 1;
                // Need some min s[i] == 2 stuff.
                if (s[i] > 1)
                {
                    sN[i] = j;
                } 
            }
            else // Get closest point outside the bandwidth
            {
                if (!sN_candidate_found)
                {
                    sN_index = j;
                    sN_value = xdi[j];
                    sN_candidate_found = true;
                } 
                else if (xdi[j] < sN_value)
                {
                    sN_index = j;
                    sN_value = xdi[j];
                }
            }
        }
        if (s[i] < 2)
        {
            s[i] = 2;
            sN[i] = sN_index;
        }   
        // sN[i] = xdi[s[i] + 1]; // wrong
        // double xeD = d[s[i] + 1]; // wrong
        // if (i == 0) dists[0] = xeD; //wrong

        // get t(i)
        for (int j = 0; j < n; j++)
        {
            if (ydi[j] < bws[1]) t[i] += 1; 
        }
        // tN[i] = ydi[t[i] + 1]; // wrong
        // double yeD = d[t[i] + 1]; // wrong
        // if (i == 0) dists[1] = yeD; //wrong

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
