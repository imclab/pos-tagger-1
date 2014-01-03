#ifndef _LOG_FUNCTIONS_H_
#define _LOG_FUNCTIONS_H

#include <math.h>
#include <algorithm>
#include <limits>
#include <vector>
#include <assert.h>

double logAdd(double logA, double logB);

double logSigma(const std::vector<double>& vec);

bool logGT(const double logA, const double logB);
#endif
