/*
 * NGram.h
 *
 *  Created on: Jul 10, 2019
 *      Author: Fritz Reichmann
 */

#ifndef SRC_NGRAM_H_
#define SRC_NGRAM_H_

#include <unordered_map>
#include <unordered_set>
#include <string>

#include "GaussianNorm.h"

class NGram {
public:
	const unsigned int _length;
	unsigned long long _count;
	std::unordered_map<std::string, unsigned long long> _NGramMap;

	NGram(const unsigned int&);
	~NGram();
	void add(const unsigned long long&, const std::string&);
};

#endif /* SRC_NGRAM_H_ */
