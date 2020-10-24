/*
 * GaussNorm.h
 *
 *  Created on: 24.10.2020
 *      Author: Fritz Reichmann
 */

#ifndef SRC_GAUSSIANNORM_H_
#define SRC_GAUSSIANNORM_H_

#include <tgmath.h>

class GaussianNorm {
public:
	const long double _mean;
	const long double _sigma;

	GaussianNorm(const long double&, const long double&);
	virtual ~GaussianNorm();
	long double lnValue(const long double&) const;
};

#endif
