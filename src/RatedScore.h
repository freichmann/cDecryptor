/*
 * RatedScore.h
 *
 *  Created on: Oct 26, 2020
 *      Author: Fritz Reichmann
 */

#ifndef SRC_RATEDSCORE_H_
#define SRC_RATEDSCORE_H_

#include "Score.h"

class RatedScore: public Score {
private:
	long double _rated;
	void compareWithNormal(const std::unordered_map<unsigned long long, GaussianNorm>&);

public:
	RatedScore();
	virtual ~RatedScore();
	RatedScore(const RatedScore& iRatedScore);
	RatedScore(const Score&, const std::unordered_map<unsigned long long, GaussianNorm>&);
	bool operator>(const RatedScore&) const;
	bool operator<(const RatedScore&) const;
	bool operator==(const RatedScore&) const;
	RatedScore& operator=(const RatedScore&);
	friend std::ostream& operator<<(std::ostream&, const RatedScore&);
	const long double value() const;
};

#endif /* SRC_RATEDSCORE_H_ */
