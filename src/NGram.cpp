/*
 * NGram.cpp
 *
 *  Created on: Jul 10, 2019
 *      Author: Fritz Reichmann
 */

#include <iostream>
#include <utility>
#include <limits>
#include "NGram.h"

NGram::NGram(unsigned int iLength) {
	_length=iLength;
	NGram::_count=0;
	NGram::_mean=std::numeric_limits<double>::quiet_NaN();
	NGram::_sigma=0;
}

NGram::~NGram() {
}

void NGram::add(unsigned long long iCount, std::string iString) {
	if (iString.length()==_length) {
		std::unordered_map<std::string,unsigned long long>::iterator i=_NGramMap.find(iString);
		if (i!=_NGramMap.end())
			i->second+=iCount;
		else
			_NGramMap.insert(std::pair<std::string, unsigned long long>(iString,iCount));
		_count+=iCount;
	}
}
