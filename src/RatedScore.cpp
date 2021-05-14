/*
 * RatedScore.cpp
 *
 *  Created on: Oct 26, 2020
 *      Author: Fritz Reichmann
 */

#include <limits>
#include "RatedScore.h"

RatedScore::RatedScore() : _rated(std::numeric_limits<double>::quiet_NaN()) {
}

RatedScore::~RatedScore() {
}

RatedScore::RatedScore(const RatedScore& iRatedScore) : Score(iRatedScore), _rated(iRatedScore._rated) {
}

RatedScore::RatedScore(const Score& iScore, const std::unordered_map<unsigned long long, GaussianNorm>& iGaussNorm) {
	Score::operator=(iScore);
	compareWithNormal(iGaussNorm);
}

bool RatedScore::operator>(const RatedScore& iThat) const {
	return std::isnan(iThat.value()) || this->value()>iThat.value();
}

bool RatedScore::operator<(const RatedScore& iThat) const {
	return !(this->value()>iThat.value() || (*this)==iThat);
}

bool RatedScore::operator==(const RatedScore& iThat) const {
	return this->value()==iThat.value();
}

bool RatedScore::operator!=(const RatedScore& iThat) const {
	return !(this->value()==iThat.value());
}

const long double RatedScore::value() const {
	return _rated;
}

RatedScore& RatedScore::operator=(const RatedScore& iThat) {
	Score::operator=(iThat);
	this->_rated=iThat._rated;
	return *this;
}

std::ostream& operator<<(std::ostream& iOut, const RatedScore& iRatedScore) {
	iOut << iRatedScore.value() << " " << (Score)iRatedScore;
	return iOut;
}

void RatedScore::compareWithNormal(const std::unordered_map<unsigned long long, GaussianNorm>& iGaussNorm) {
	_rated=0;
	for (std::unordered_map<unsigned long long, long double>::const_iterator aI=_metrics.cbegin(); aI!=_metrics.cend(); ++aI)
		_rated+=iGaussNorm.at(aI->first).lnValue(aI->second);
}
