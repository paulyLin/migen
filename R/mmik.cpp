#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List mmikCpp(List ctsLevels, NumericVector bws, 
		NumericVector N, NumericVector cts) 
{

	int levels = ctsLevels.size();
	int n = N[0];

	NumericVector outvec(levels);

	// Accumulators
	NumericVector k(n); //make first two integers?
	NumericVector ka(n);
	NumericVector d(n);

	int kk = 0;
	for (int i = 0; i < levels; i++)
	{
		SEXP ll = ctsLevels[i];
		NumericVector y(ll);
		outvec[i] = sum(y);
		
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

			k[kk] = std::max(k[kk] - 1.0, 1.0);
			// This is why I should make the vectors integer

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

	//outvec[0] = sum(bws);
	return List::create(Named("outvec") = outvec, Named("k") = k, Named("d") = d, Named("ka") = ka);
}	 
