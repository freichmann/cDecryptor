/*
 * Score.cpp
 *
 *  Created on: Oct 23, 2020
 *      Author: Fritz Reichmann
 */

#include "Score.h"

Score::Score::Score() {
}

Score::Score::Score(const std::unordered_map<unsigned long long, NGram*>& ipNorms, const std::string &iClear) {
	computeMetrics(iClear, ipNorms);
}

Score::Score::Score(const std::unordered_map<unsigned long long, NGram*>& ipNorms, const std::unordered_map<unsigned long long, GaussianNorm>& iGaussNorm, const std::string &iClear) {
	computeMetrics(iClear, ipNorms);
	compareWithNormal(iGaussNorm);
}

Score::Score::Score(const Score& iScore) : _metrics(iScore._metrics), _rated(iScore._rated) {
}

Score::Score::~Score() {
}

const std::unordered_map<unsigned long long, long double> Score::Score::values() const {
	return _metrics;
}

const long double Score::Score::likelihood() const {
	if (std::isnan(_rated))
		throw std::string("Score not rated!");
	else
		return _rated;
}

bool Score::operator>(const Score& iThat) const {
	return this->likelihood()>iThat.likelihood();
}

bool Score::operator<(const Score& iThat) const {
	return this->likelihood()<iThat.likelihood();
}

std::ostream& operator<<(std::ostream& iOut, const Score& iScore) {
	bool aFirst=true;

	iOut << iScore.likelihood() << " ";

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
	this->_rated=iThat._rated;

	return *this;
}

long double logFac(const unsigned long long& iK) {
	long double aLogFac=0;
	for (unsigned long long i=2; i<=iK; i++)
		aLogFac+=logl(i);
	return aLogFac;
}

void Score::Score::compareWithNormal(const std::unordered_map<unsigned long long, GaussianNorm>& iGaussNorm) {
	_rated=0;
	for (std::unordered_map<unsigned long long, long double>::const_iterator aI=_metrics.cbegin(); aI!=_metrics.cend(); ++aI)
		_rated+=iGaussNorm.at(aI->first).lnValue(aI->second);
}

void Score::Score::computeMetrics(const std::string& iCandidate, const std::unordered_map<unsigned long long, NGram*>& ipNormNGrams) {
	_metrics.clear();

	std::unordered_map<std::string, unsigned long long> aSolutionCandidateNGram;

	for (std::unordered_map<unsigned long long, NGram*>::const_iterator aNormNGram=ipNormNGrams.cbegin(); aNormNGram!=ipNormNGrams.cend(); ++aNormNGram) {
		long double aLoopScore=0;
		const unsigned int aLength=aNormNGram->second->_length;

		aSolutionCandidateNGram.clear();
		for (unsigned int i=0; i+aLength<=iCandidate.length(); i++) {
			const std::string aSub=iCandidate.substr(i, aLength);
			std::unordered_map<std::string, unsigned long long>::iterator b=aSolutionCandidateNGram.find(aSub);
			if (b!=aSolutionCandidateNGram.end())
				b->second++;
			else
				aSolutionCandidateNGram.insert(std::pair<std::string, unsigned long long>(aSub, 1));
		}

		const long double aNeverSeen=logl(2)/(long double)aNormNGram->second->_count;
		for (std::unordered_map<std::string, unsigned long long>::const_iterator b=aSolutionCandidateNGram.cbegin(); b!=aSolutionCandidateNGram.cend(); ++b) {
			long double aP;
			std::unordered_map<std::string, unsigned long long>::const_iterator j=aNormNGram->second->_NGramMap.find(b->first);
			if (j!=aNormNGram->second->_NGramMap.end())
				aP=(long double)j->second/(long double)aNormNGram->second->_count;
			else
				aP=aNeverSeen;

			aLoopScore+=b->second*logl(aP)-logFac(b->second); // Metric: Poisson(b)-Poisson(0) (+const)
		}

		_metrics.insert(std::pair<unsigned long long, long double>(aLength,aLoopScore));
	}
}
