/*
 * log functions like log add and sigma
 */
#include "log_functions.h"
#include <assert.h>
#include <iostream>
/*
 * Takes given log(a) and log(b) it computes log(a+b)
 */
double logAdd(double logA, double logB) {
	//code doesn't handle the case of +inf, it is ASSUMED that +inf won't happen
	bool a_inf = isinf(logA);
	bool b_inf = isinf(logB);
	bool a_nan = isnan(logA);
	bool b_nan = isnan(logB);
	
	if (a_nan || b_nan) {
		std::cout<<"a ="<<logA<<std::endl;
		std::cout<<"b ="<<logB<<std::endl;
	}
	/*
	std::cout<<"loga="<<logA<<std::endl;
	std::cout<<"logb="<<logB<<std::endl;
	*/


	if (0 != a_inf && 0 != b_inf) {
		//since we are assuming we dont deal with +inf, -inf and -inf is log (0+0) which is still -inf
		return logA;
	}
	/*
	log of 0 is -inf, so instead of addition, we return the other value
	the thinking behind this is if we are adding -inf and some other number
	then that means we have log(0) and log(some number) and want log(0 + some number)
	which is just log(some number)
	*/
	if (a_inf > 0 && b_inf == 0) {
		return logB;
	}
	if (a_inf == 0 && b_inf > 0) {
		return logA;
	}

	if (logB > logA) {
		std::swap(logA, logB);
	}

	if (logA <= -1 * std::numeric_limits<double>::max()) {
		return logA;
	}

	double negDiff = logB - logA;
	//if one is b is 20 orders of magnitude smaller just return a
	if (negDiff < -20.0) {
		return logA;
	}
	return logA + log(1.0 + exp(negDiff));
}

/*
 * computes the log of the sum of the vector given a vector of numbers in log space
 */
double logSigma(const std::vector<double>& vec) {
	//handle errors that might happen
	assert (vec.size() != 0);
	
	double log_sum = vec[0];
	for (unsigned int pos = 1; pos < vec.size(); pos++) {
		log_sum = logAdd(log_sum, vec[pos]);
	}
	return log_sum;
}

bool logGT(const double logA, const double logB) {
	//implemetation of the > operator, except handles -inf
	bool a_inf = isinf(logA);
	bool b_inf = isinf(logB);
	if (a_inf && b_inf) {
		//assumes -inf, -inf == -inf
		return false;
	}
	if (a_inf) {
		//-inf > some number ?
		return false;
	}
	if (b_inf) {
		//some number > -inf?
		return true;
	}
	else {
		return logA > logB;
	}
}
