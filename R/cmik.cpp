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

    // Finding furthest point within bw.
    int sN_inbw_index = 0;
    double sN_inbw_value = 0.0;
    bool sN_inbw_candidate_found = false;
    int tN_inbw_index = 0;
    double tN_inbw_value = 0.0;
    bool tN_inbw_candidate_found = false;

    // Finding closest point outside bw.
    int sN_outbw_index = 0;
    double sN_outbw_value = 0.0;
    bool sN_outbw_candidate_found = false;
    int tN_outbw_index = 0;
    double tN_outbw_value = 0.0;
    bool tN_outbw_candidate_found = false;

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
 
        // get s(i) and sN(i)
        sN_inbw_index = 0;
        sN_inbw_value = 0;
        sN_inbw_candidate_found = false;
        sN_outbw_index = 0;
        sN_outbw_value = 0;
        sN_outbw_candidate_found = false;
        for (int j = 0; j < n; j++)
        {
            if (xdi[j] < bws[0])
            {  
                s[i] += 1;
                
                // Find furthest point inside bw.
                if(!sN_inbw_candidate_found)
                {
                    sN_inbw_index = j;
                    sN_inbw_value = xdi[j];
                    sN_inbw_candidate_found = true;
                }
                else if (xdi[j] > sN_inbw_value) // N.b. greater than
                {
                    sN_inbw_index = j;
                    sN_inbw_value = xdi[j];
                }
            }
            else
            {
                // Find closest point outside bw.
                // Used in case where s[i] == 1, and
                // we need to set it to s[i] == 2.
                if(!sN_outbw_candidate_found)
                {
                    sN_outbw_index = j;
                    sN_outbw_value = xdi[j];
                    sN_outbw_candidate_found = true;
                }
                else if (xdi[j] < sN_outbw_value) // N.b. less than
                {
                    sN_outbw_index = j;
                    sN_outbw_value = xdi[j];
                }
            }
         }
        if (s[i] > 1)
        {
            sN[i] = sN_inbw_index;
        } 
        else
        {
            s[i] = 2;
            sN[i] = sN_outbw_index;
        }

        // get t(i)
//        for (int j = 0; j < n; j++)
//        {
//            if (ydi[j] < bws[1]) t[i] += 1; 
//        }

        // get t(i) snd tN(i)
        tN_inbw_index = 0;
        tN_inbw_value = 0;
        tN_inbw_candidate_found = false;
        tN_outbw_index = 0;
        tN_outbw_value = 0;
        tN_outbw_candidate_found = false;
        for (int j = 0; j < n; j++)
        {
            if (ydi[j] < bws[1])
            {  
                t[i] += 1;
                
                // Find furthest point inside bw.
                if(!tN_inbw_candidate_found)
                {
                    tN_inbw_index = j;
                    tN_inbw_value = ydi[j];
                    tN_inbw_candidate_found = true;
                }
                else if (ydi[j] > tN_inbw_value) // N.b. greater than
                {
                    tN_inbw_index = j;
                    tN_inbw_value = ydi[j];
                }
            }
            else
            {
                // Find closest point outside bw.
                // Used in case where s[i] == 1, and
                // we need to set it to s[i] == 2.
                if(!tN_outbw_candidate_found)
                {
                    tN_outbw_index = j;
                    tN_outbw_value = ydi[j];
                    tN_outbw_candidate_found = true;
                }
                else if (ydi[j] < tN_outbw_value) // N.b. less than
                {
                    tN_outbw_index = j;
                    tN_outbw_value = ydi[j];
                }
            }
         }
        if (t[i] > 1)
        {
            tN[i] = tN_inbw_index;
        } 
        else
        {
            t[i] = 2;
            tN[i] = tN_outbw_index;
        }

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
