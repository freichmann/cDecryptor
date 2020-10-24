/*
 * GaussNorm.cpp
 *
 *  Created on: 24.10.2020
 *      Author: Fritz Reichmann
 */

#include "GaussNorm.h"

long double GaussNorm::GaussNorm::lnValue(const long double& iX) const {
	return -0.5*powl((iX-_mean)/_sigma,2)-logl(_sigma)-logl(sqrtl(2*M_PI));
}

GaussNorm::GaussNorm(const long double& iMean, const long double& iSigma) : _mean(iMean), _sigma(iSigma) {
}

GaussNorm::~GaussNorm() {
}

