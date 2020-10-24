/*
 * Score.h
 *
 *  Created on: Oct 23, 2020
 *      Author: Fritz Reichmann
 */

#ifndef SRC_SCORE_H_
#define SRC_SCORE_H_

#include <string>
#include <unordered_map>
#include <tgmath.h>
#include <iostream>
#include <limits>

#include "NGram.h"

class Score {
private:
	std::unordered_map<unsigned long long, long double> _metrics;
	long double _rated=std::numeric_limits<long double>::quiet_NaN();

	void computeMetrics(const std::string&, const std::unordered_map<unsigned long long, NGram*>&);
	void compareWithNormal(const std::unordered_map<unsigned long long, GaussianNorm>&);

public:
	Score();
	Score(const Score&);
	Score(const std::unordered_map<unsigned long long, NGram*>&, const std::string&);
	Score(const std::unordered_map<unsigned long long, NGram*>&, const std::unordered_map<unsigned long long, GaussianNorm>&, const std::string&);
	virtual ~Score();

	bool operator>(const Score&) const;
	bool operator<(const Score&) const;
	Score& operator=(const Score&);
	friend std::ostream& operator<<(std::ostream& out, const Score&);

	const std::unordered_map<unsigned long long, long double> values() const;
	const long double likelihood() const;
};
#endif
