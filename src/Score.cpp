/*
 * Score.cpp
 *
 *  Created on: Oct 23, 2020
 *      Author: Fritz Reichmann
 */

#include "Score.h"

Score::Score() {
}

Score::Score(const std::unordered_map<unsigned long long, NGram*>& ipNorms, const std::string &iClear) {
	computeMetrics(iClear, ipNorms);
}

Score::Score(const Score& iScore) : _metrics(iScore._metrics) {
}

Score::~Score() {
}

const std::unordered_map<unsigned long long, long double> Score::values() const {
	return _metrics;
}

std::ostream& operator<<(std::ostream& iOut, const Score& iScore) {
	bool aFirst=true;

	for (std::unordered_map<unsigned long long, long double>::const_iterator aI=iScore._metrics.cbegin(); aI!=iScore._metrics.cend(); ++aI) {
		if (!aFirst)
			iOut << " ";
		else
			aFirst=false;
		iOut << aI->first << ":" << aI->second;
	}
	return iOut;
}

Score& Score::operator=(const Score& iThat) {
	this->_metrics=std::unordered_map<unsigned long long, long double>(iThat._metrics);

	return *this;
}

void Score::computeMetrics(const std::string& iCandidate, const std::unordered_map<unsigned long long, NGram*>& ipNormNGrams) {
	_metrics.clear();

	for (std::unordered_map<unsigned long long, NGram*>::const_iterator aNormNGram=ipNormNGrams.cbegin(); aNormNGram!=ipNormNGrams.cend(); ++aNormNGram) {
		std::unordered_map<std::string, unsigned long long> aSolutionCandidateNGram;

		for (unsigned int i=0; i+aNormNGram->second->_length<=iCandidate.length(); i++) {
			const std::string aSub=iCandidate.substr(i, aNormNGram->second->_length);
			std::unordered_map<std::string, unsigned long long>::iterator b=aSolutionCandidateNGram.find(aSub);
			if (b!=aSolutionCandidateNGram.end())
				b->second++;
			else
				aSolutionCandidateNGram.insert(std::pair<std::string, unsigned long long>(aSub, 1));
		}

		long double aScore=0;
		for (std::unordered_map<std::string, unsigned long long>::const_iterator b=aSolutionCandidateNGram.cbegin(); b!=aSolutionCandidateNGram.cend(); ++b) {
			std::unordered_map<std::string, unsigned long long>::const_iterator j=aNormNGram->second->_NGramMap.find(b->first);
			long double aLnP;
			if (j!=aNormNGram->second->_NGramMap.end())
				aLnP=logl(j->second)-logl(aNormNGram->second->_count);
			else
				aLnP=logl(logl(2))-logl(aNormNGram->second->_count);

			aScore+=b->second*aLnP-lgammal(b->second+1); // Metric: Poisson(b)-Poisson(0) (+const)
		}

		_metrics.insert(std::pair<unsigned long long, long double>(aNormNGram->second->_length, aScore));
	}
}
