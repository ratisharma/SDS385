// Spencer Woody, 05 Oct 2016

#include <RcppEigen.h>
#include <algorithm>    // std::max

using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::LLT;
using Eigen::Lower;
using Eigen::MatrixXi;
using Eigen::Upper;
using Eigen::VectorXd;
using Eigen::VectorXi;
using Eigen::SparseVector;
typedef Eigen::MappedSparseMatrix<double>  MapMatd;
typedef Map<MatrixXi>  MapMati;
typedef Map<VectorXd>  MapVecd;
typedef Map<VectorXi>  MapVeci;

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
SEXP sgdcppr(MapMatd X, VectorXd y, VectorXd M, double stepsize0, int numpasses,
			VectorXd beta0, double lambda = 0.1, double edecay = 0.01) {
				// X is the design matrix, with each column an observation
				// y is the vector of responses
				// M is the vector of sample sizes (just ones in case of logistic regression)
				// stepsize0 is the default (first) step size for AdaGrad
				// numpasses is the number of times to cycle through the whole dataset
				// beta0 is the initial guess for beta
				// lambda is the ell-2 penalization parameter
				// edecay is the weight assigned to the new update in EMA of likelihood
				
				///////////////
				// SECTION 1 // Declare Variables
 				///////////////
				
				int N = X.cols();
				int P = X.rows();
				
				// initialize parameters
				// x is used for each observation (column) in X
				// yhat is predicted y value
				// ydelta is difference between current y and yhat
				// dotprod is x_i^T * beta, edotprod is e^dotprod
				// gsqrt is used for sqrt of historical SS of gradients
				SparseVector<double> x(P);
				int j, k; // iteration counters
				double dotprod, edotprod, yhat, ydelta, gsqrt;
				
				
				// w_hat is used for an initial guess for te intercept, alpha
				// alpha is estimated intercept term
				// grad_j is the gradient at the jth feature in X at current observation
				// g0sq is the square of the gradient for the intercept term
				// gam is AdaGrad step from previous time beta(j) is updated
				// sum_penalty is accumulated penalty for lazy updating
				double w_hat = y.sum() / M.sum();
				double alpha = log(w_hat);
				double grad_j = 0.0;
				double g0sq = 0.0;
				double gam;
				double sum_penalty = 0;
				
				// Init beta
				// Init Gsq, vector for AdaGrad (1e-3 added for numerical stability)
				// adj_grad is gradient divided by sqrt(Gsq)
				VectorXd beta(P);
				VectorXd Gsq(P);
				VectorXd adj_grad(P);
				for (int j = 0; j < P; j++) {
					Gsq(j) = 1e-3;
					adj_grad(j) = 1e-7;
					beta(j) = beta0(j);
				}
				
				// Vector for last time feature j is updated (for lazy updating)
				NumericVector last_update(P, 0.0);
				
				// Hold current EMA of loglikelihood avg, and vector for storing historic nll
				double loglikavg = 0.0;
				NumericVector logliktrace(numpasses*N, 0.0);
				
				// Global iteration counter
				k = 0; 
				
				///////////////
				// SECTION 2 // Perform SGD with AdaGrad
 				///////////////
				
				
				// LOOP # 1
				// Loop for number of passes to go through
				for(int pass = 0; pass < numpasses; pass++) {
					
					// LOOP # 2
					// Outer loop over all observations (i.e. columns in design matrix)
					for(int i = 0; i < N; i++) {
						
						// For yhat, predicted value of y at current estimate of beta
						x = X.innerVector(i); // grab ith observation of x (but only 
											  // the active features)
						dotprod = alpha + x.dot(beta);
						edotprod = exp(dotprod);
						yhat = M[i] * edotprod / (1 + edotprod);
						
						// Update EMA of loglik, add it to vector
						loglikavg = (1 - edecay) * loglikavg + edecay * (M[i]*log(1 + edotprod) - y[i]*dotprod);
						logliktrace[k] = loglikavg;
						
						// Update intercept alpha
						ydelta = y[i] - yhat;
						g0sq += ydelta * ydelta; 
						
						alpha += (stepsize0 / sqrt(g0sq)) * ydelta;
						
						// LOOP # 3
						// Inner loop; Parse over only active features within this observation
						for (SparseVector<double>::InnerIterator it(x); it; ++it) {
							
							// Grad index of feature
							j = it.index();
							
							// STEP 1
							// Lazy update for accumulated penalty
							
							// number of iterations since last update
							//Cap maximum number of recursive updates at 5, for numeric stability.
							double skip = k - last_update(j); 
							if (skip > 5){skip = 5;}
							last_update(j) = k; // update last_update with current global iterator
							
							// Add up penalties
							gam = stepsize0 * adj_grad(j);
							sum_penalty = beta(j) * ( (1 - pow(1 + lambda * gam, skip)) / (1 - lambda * gam));

							// Subtract penalties from beta(j)
							beta(j) -= sum_penalty;
							
							// ell-2 penalty for gradient
							double ell2penalty = 2 * lambda * beta(j);
							
							// Compute gradient for this feature at this observation
							grad_j = - (ydelta * it.value()) + ell2penalty;
							
							// Update Gsq vector for this feature
							Gsq(j) += grad_j * grad_j;
							
							// Compute adjusted gradient (i.e., AdaGrad)
							gsqrt = sqrt(Gsq(j));
							adj_grad(j) = grad_j / gsqrt;
							
							// Update beta
							beta(j) -= adj_grad(j) * stepsize0;
						}
						k++;
					}
					
					
					
				}
				// Penalize betas we haven't touched in a while
				for (int j=0; j<P; ++j){
					//Using (k-1) since last_updated indexes from 0, and n is based on counting rows from 1.
					double skip = (k-1) - last_update(j); 
					if (skip > 5){  skip = 5;}			
		
					//Calculate accum penalty.
					gam = stepsize0*adj_grad(j);	
					sum_penalty = beta(j) * ( (1 - pow(1 + lambda * gam, skip)) / (1 - lambda * gam));
		
					//Update beta_js to add accum penalty.
					beta(j) -= sum_penalty; 
				}			
				
				// Finally, return alpha, beta, and loglikelihood
				return List::create(Named("alpha") = alpha,
									Named("beta")  = beta,
									Named("logliktrace") = logliktrace);
}