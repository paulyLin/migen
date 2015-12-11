#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector mmikCpp(List ctsLevels, NumericVector bws, 
        NumericVector N, NumericVector cts, NumericVector Ni) 
{

    int levels = ctsLevels.size();
    int n = N[0];

    // Accumulators
    IntegerVector k(n); //make first two integers?
    IntegerVector ka(n);
    NumericVector d(n);

    int kk = 0;
    for (int i = 0; i < levels; i++)
    {
        SEXP ll = ctsLevels[i];
        NumericVector y(ll);

        for (int j = 0; j < y.size(); j++)
        {
            for (int l = 0; l < y.size(); l++)
            {
                if (std::abs(y[l] - y[j]) < bws[i])
                {
                    k[kk] += 1; // Number of points within BW

                    // Max distance to point within BW
                    if (d[kk] < std::abs(y[l] - y[j]))
                    {
                        d[kk] = std::abs(y[l] - y[j]);
                    }
                }

            }

            k[kk] = std::max(k[kk] - 1, 1);
            // This is why I should make the vectors integer

            if (d[kk] == 0 | NumericVector::is_na(d[kk]))
            {
                if (j == 0) d[kk] = std::abs(y[j] - y[j + 1]);
                else if (j == y.size() - 1) d[kk] = std::abs(y[j] - y[j - 1]); 
                else
                {
                    d[kk] = std::max(std::abs(y[j] - y[j + 1]), 
                            std::abs(y[j] - y[j - 1]));
                }

            }

            // Count points in all groups closer than d[kk]
            // Some redundancy here
            for (int l = 0; l < cts.size(); l++)
            {
                if (std::abs(cts[l] - y[j]) < d[kk])
                {
                    ka[kk] += 1;
                }
            }

            kk += 1;
        }
    }

    NumericVector MI(1);
    MI = mean(digamma(k) - digamma(ka)) 
        + sum(Ni * (-digamma(Ni))) / N
        + digamma(N);

    return MI;
}	 
