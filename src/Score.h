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

#include "NGram.h"

class Score {
private:
	void computeMetrics(const std::string&, const std::unordered_map<unsigned long long, NGram*>&);

protected:
	std::unordered_map<unsigned long long, long double> _metrics;

public:
	Score();
	Score(const Score&);
	Score(const std::unordered_map<unsigned long long, NGram*>&, const std::string&);
	virtual ~Score();

	Score& operator=(const Score&);
	friend std::ostream& operator<<(std::ostream& out, const Score&);

	const std::unordered_map<unsigned long long, long double> values() const;
};
#endif
