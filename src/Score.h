/*
 * Score.h
 *
 *  Created on: Oct 23, 2020
 *      Author: freichmann
 */

#ifndef SRC_SCORE_H_
#define SRC_SCORE_H_

#include <string>
#include <unordered_map>
#include <tgmath.h>
#include <iostream>


#include "NGram.h"

class Score {
private:
	long double logFac(const unsigned long long&);
	long double lnGauss(const long double& iX, const long double& iMean, const long double& iSigma);
	std::ostringstream _log;
	long double _score;
	long double score(const std::string&, const std::unordered_map<unsigned long long, NGram*>&);

public:
	Score();
	Score(const Score&);
	Score(const std::unordered_map<unsigned long long, NGram*>&, const std::string&);
	virtual ~Score();

	bool operator>(const Score&) const;
	bool operator<(const Score&) const;
	Score& operator=(const Score&);
	friend std::ostream& operator<<(std::ostream& out, const Score&);

	long double rate() const;
};
#endif /* SRC_SCORE_H_ */
