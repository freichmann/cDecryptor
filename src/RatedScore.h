/*
 * RatedScore.h
 *
 *  Created on: Oct 26, 2020
 *      Author: freichmann
 */

#ifndef SRC_RATEDSCORE_H_
#define SRC_RATEDSCORE_H_

#include <limits>

#include "Score.h"

class RatedScore: public Score {
private:
	long double _rated=std::numeric_limits<long double>::quiet_NaN();
	void compareWithNormal(const std::unordered_map<unsigned long long, GaussianNorm>&);

public:
	RatedScore();
	virtual ~RatedScore();
	RatedScore(const RatedScore& iRatedScore);
	RatedScore(const Score&, const std::unordered_map<unsigned long long, GaussianNorm>&);
	bool operator>(const RatedScore&) const;
	bool operator<(const RatedScore&) const;
	friend std::ostream& operator<<(std::ostream&, const RatedScore&);
	RatedScore& operator=(const RatedScore&);
	const long double value() const;
};

#endif /* SRC_RATEDSCORE_H_ */
