/*
 * Score.cpp
 *
 *  Created on: Oct 23, 2020
 *      Author: freichmann
 */

#include "Score.h"

Score::Score::Score() {
	_score=0.0;
}

Score::Score::Score(const std::unordered_map<unsigned long long, NGram*>& ipNorms, const std::string &iClear) : _score(score(iClear, ipNorms)) {
}

Score::Score::Score(const Score& iScore) : _score(iScore._score) {
	_log.str(iScore._log.str());
}

Score::Score::~Score() {
}

long double Score::Score::rate() const {
	return _score;
}

bool Score::operator>(const Score& iThat) const {
	return this->_score>iThat._score;
}

bool Score::operator<(const Score& iThat) const {
	return this->_score<iThat._score;
}

std::ostream& operator<<(std::ostream& iOut, const Score& iScore) {
	iOut << iScore.rate() << " " << iScore._log.str();
	return iOut;
}

Score& Score::operator=(const Score& iThat) {
	this->_score=iThat._score;
	_log=std::ostringstream(iThat._log.str());
	return *this;
}

long double Score::Score::lnGauss(const long double& iX, const long double& iMean, const long double& iSigma) {
	return -0.5*powl((iX-iMean)/iSigma,2)-logl(iSigma)-logl(sqrtl(2*M_PI));
}

long double logFac(const unsigned long long& iK) {
	long double aLogFac=0;
	for (unsigned long long i=2; i<=iK; i++)
		aLogFac+=logl(i);
	return aLogFac;
}

long double Score::Score::score(const std::string& iCandidate, const std::unordered_map<unsigned long long, NGram*>& ipNormNGrams) {
	long double aScore=0;
	bool aFirstLog=true;

	_log.str("");
	_log.clear();

	std::unordered_map<std::string, unsigned long long> aSolutionCandidateNGram;

	for (std::unordered_map<unsigned long long, NGram*>::const_iterator aNormNGram=ipNormNGrams.cbegin(); aNormNGram!=ipNormNGrams.cend(); ++aNormNGram) {
		long double aLoopScore=0;
		const unsigned int aLength=aNormNGram->second->_length;

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

			aLoopScore+=b->second*logl(aP)-logFac(b->second);
		}

		if (aNormNGram->second->_mean<0) {
			const long double aLnGaussScore=lnGauss(aLoopScore, aNormNGram->second->_mean, aNormNGram->second->_sigma);
			aScore+=aLnGaussScore;
			if (!aFirstLog)
				_log << " ";
			else
				aFirstLog=false;
			_log << aLength << ":" << aLnGaussScore;
		} else
			aScore+=aLoopScore;

		aSolutionCandidateNGram.clear();
	}

	return aScore;
}
