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

class NGram {
public:
	unsigned int _length;
	unsigned long _count;
	long double _mean;
	long double _sigma;
	std::unordered_map<std::string, unsigned long long> _NGramMap;

	NGram(unsigned int);
	~NGram();
	void add(unsigned long long, std::string);
};

#endif /* SRC_NGRAM_H_ */
