/*
 * GaussNorm.h
 *
 *  Created on: 24.10.2020
 *      Author: Fritz Reichmann
 */

#ifndef SRC_GAUSSIANNORM_H_
#define SRC_GAUSSIANNORM_H_

#include <tgmath.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

class GaussianNorm {
public:
	const long double _mean;
	const long double _sigma;

	GaussianNorm(const long double&, const long double&);
	virtual ~GaussianNorm();
	long double lnValue(const long double&) const;
};

#endif
