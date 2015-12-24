#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector cmikCpp(NumericVector bws, 
        NumericVector N, NumericVector X, NumericVector Y, NumericVector kmax) 
{
    int n = N[0];
    int km = kmax[0];


    // X differences
    NumericMatrix xdiffs(n);
    for (int j = 0; j < n; j++)
    {
        for (int i = j; i < n; i++)
        {
            if (j == i) continue; // Keep 0s on diagonal.

            xdiffs(i, j) = std::abs(X[i] - X[j]);
            xdiffs(j, i) = xdiffs(i, j); // Symmetry of abs().
        }
    }

    // Y differences
    NumericMatrix ydiffs(n);
    for (int j = 0; j < n; j++)
    {
        for (int i = j; i < n; i++)
        {
            if (j == i) continue; // Keep 0s on diagonal.

            ydiffs(i, j) = std::abs(Y[i] - Y[j]);
            ydiffs(j, i) = ydiffs(i, j); // Symmetry of abs().
        }
    }

    // Max norm distances
    NumericMatrix distmat(n);
    for (int j = 0; j < n; j++)
    {
        for (int i = j; i < n; i++)
        {
            if (j == i) continue; // Keep 0s on diagonal.

            distmat(i, j) = std::max(xdiffs(i, j), ydiffs(i, j));
            distmat(j, i) = distmat(i, j); // Symmetry.
        }
    }

    int s = 0;
    int sN = 0;
    int t = 0;
    int tN = 0;
    double eD = 0; 

    // Need to keep k, l and m for all i.
    NumericVector k(n); 
    NumericVector l(n);
    NumericVector m(n);

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

    double xeD = 0;
    double yeD = 0;

    int kx = 0;
    int ky = 0;

    // Main loop
    for (int i = 0; i < n; i++)
    {
        s = 0;
        sN = 0;
        t = 0;
        tN = 0;
        eD = 0; 

        // N.b. R matrices are column-major.
        NumericVector xdi = xdiffs(_, i);
        NumericVector ydi = ydiffs(_, i);
        NumericVector d = distmat(_, i);
 
        // get s(i) and sN(i)
        sN_inbw_index = 0;
        sN_inbw_value = 0;
        sN_inbw_candidate_found = false;
        sN_outbw_index = 0;
        sN_outbw_value = 0;
        sN_outbw_candidate_found = false;
        xeD = 0; 
        yeD = 0;
        kx = 0;
        ky = 0;
        for (int j = 0; j < n; j++)
        {
            if (xdi[j] < bws[0])
            {  
                s += 1;
                
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
                // Used in case where s == 1, and
                // we need to set it to s == 2.
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
        // Must have s > 1
        if (s > 1)
        {
            sN = sN_inbw_index;
        } 
        else
        {
            s = 2;
            sN = sN_outbw_index;
        }

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
                t += 1;
                
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
                // Used in case where t == 1, and
                // we need to set it to t == 2.
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
        // Must have t > 1 
        if (t > 1)
        {
            tN = tN_inbw_index;
        } 
        else
        {
            t = 2;
            tN = tN_outbw_index;
        }

        xeD = d(sN);
        yeD = d(tN);

        // Find kx and ky
        for (int j = 0; j < n; j++)
        {
            if (d[j] < xeD)
            {
                kx += 1;
            }

            if (d[j] < yeD)
            {
                ky += 1;
            }
        }

        // Limit k to kmax. If this happens we then need
        // to find the corresponding distance from d.
        if (km < kx && km < ky)
        {
            k[i] = km;
            std::nth_element(d.begin(), d.begin() + km, d.end());
            // Gets (km + 1)th smallest distance (C++ indexes from 0, R from 1).
            eD = d[km];
        }
        else
        {
            // If k <= kmax we already have the correct distance.
            if (kx <= ky)
            {
                k[i] = kx;
                eD = xeD;
            }
            else
            {
                k[i] = ky;
                eD = yeD;
            }
        }

        // Get l and m
        for (int j = 0; j < n; j++)
        {
            if (xdi[j] < eD)
            {
                l[i] += 1;
            }

            if (ydi[j] < eD)
            {
                m[i] += 1;
            }
        }
    } // Main loop

    NumericVector MI(1);

    MI = digamma(N) + mean(digamma(k)) - mean(digamma(l) + digamma(m));

    return(MI);
} 
