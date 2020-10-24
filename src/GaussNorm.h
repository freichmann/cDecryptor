/*
 * GaussNorm.h
 *
 *  Created on: 24.10.2020
 *      Author: Fritz Reichmann
 */

#ifndef SRC_GAUSSNORM_H_
#define SRC_GAUSSNORM_H_

#include <tgmath.h>

class GaussNorm {
public:
	const long double _mean;
	const long double _sigma;

	GaussNorm(const long double&, const long double&);
	virtual ~GaussNorm();
	long double lnValue(const long double&) const;
};

#endif
