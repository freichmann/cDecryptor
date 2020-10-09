#include <getopt.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <list>
#include <vector>
#include <random>
#include <iomanip>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <tgmath.h>

#include "NGram.h"

//Globals
std::mutex aBestScoreMux;
std::mutex aOutputMux;
long double aGlobalBestScore;
std::string aGlobalBestSolution;

long double logFac(unsigned long long iK) {
	long double aLogFac=0;
	for (unsigned long long i=2; i<=iK; i++)
		aLogFac+=logl(i);
	return aLogFac;
}

long double lnGauss(const long double iX, const long double iMean, const long double iSigma) {
	return -0.5*powl((iX-iMean)/iSigma,2)-logl(iSigma)-logl(sqrtl(2*M_PI));
}

long double score(const std::string *ipCandidate, const std::unordered_map<unsigned long long, NGram*> *ipNorms, std::ostringstream *iLog) {
	long double aScore=0;

	if (iLog!=NULL)
		iLog->str("");

	std::unordered_map<std::string, unsigned long long> *aObserved=new std::unordered_map<std::string, unsigned long long>();

	for (std::unordered_map<unsigned long long, NGram*>::const_iterator aNorm=ipNorms->begin(); aNorm!=ipNorms->end(); ++aNorm) {
		long double aLoopScore=0;
		unsigned int aLength=aNorm->second->_length;

		for (unsigned int i=0; i+aLength<=ipCandidate->length(); i++) {
			std::string aSub=ipCandidate->substr(i , aLength);
			std::unordered_map<std::string, unsigned long long>::iterator b=aObserved->find(aSub);
			if (b!=aObserved->end())
				b->second++;
			else
				aObserved->insert(std::pair<std::string, unsigned long long>(aSub, 1));
		}

		const long double aNeverSeen=logl(2)/(long double)aNorm->second->_count;
		for (std::unordered_map<std::string, unsigned long long>::iterator b=aObserved->begin(); b!=aObserved->end(); ++b) {
			long double aP;
			std::unordered_map<std::string, unsigned long long>::iterator j=aNorm->second->_NGramMap->find(b->first);
			if (j!=aNorm->second->_NGramMap->end())
				aP=(long double)j->second/(long double)aNorm->second->_count;
			else
				aP=aNeverSeen;

			aLoopScore+=b->second*logl(aP)-logFac(b->second);
		}

		if (aNorm->second->_mean<0) {
			long double aLnGaussScore=lnGauss(aLoopScore, aNorm->second->_mean, aNorm->second->_sigma);
			aScore+=aLnGaussScore;
			if (iLog!=NULL)
				(*iLog) << aLength << ":" << aLnGaussScore << " ";
		} else
			aScore+=aLoopScore;

		aObserved->clear();
	}

	delete aObserved;
	return aScore;
}

void readCipher(const std::string& aCipherfile, std::string *ipCipher) {
	std::ifstream myfile(aCipherfile);
	if (myfile.is_open()) {
		while (!myfile.eof()) {
			std::string aLine;
			myfile >> aLine;
			*ipCipher += aLine;
		}
		myfile.close();
	}
}

void readNorms(std::list<std::string> afiles, std::unordered_map<unsigned long long, NGram*>* aNorms) {
	for (std::list<std::string>::iterator i = afiles.begin(); i != afiles.end(); ++i) {
		std::cout << "Reading norm file " << *i << std::endl;
		std::ifstream myfile(*i);
		if (myfile.is_open()) {
			while (!myfile.eof()) {
				std::string aString;
				myfile >> aString;
				transform(aString.begin(), aString.end(), aString.begin(), ::tolower);
				unsigned long long aCount;
				myfile >> aCount;
				if (aString.length() > 0 && aCount > 0) {
					std::unordered_map<unsigned long long, NGram*>::iterator i=aNorms->find(aString.length());
					NGram* aNorm;
					if (i == aNorms->end()) {
						aNorm = new NGram(aString.length());
						aNorms->insert(std::pair<unsigned long long, NGram*>(aString.length(), aNorm));
					} else
						aNorm = (*i).second;

					aNorm->add(aCount, aString);
				}
			}
			myfile.close();
		} else
			std::cout << "Failed to open file" << (*i) << std::endl;
	}
}

void randomizeMap(std::unordered_map<char, char> *ipUniqueSymbols, const long double &iRandom) {
	std::default_random_engine aGenerator;
	aGenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<unsigned int> aCharDistribution((unsigned int)'a', (unsigned int)'z');
	std::uniform_real_distribution<long double> aDoubleDistribution(0,1);

	for (std::unordered_map<char, char>::iterator i=ipUniqueSymbols->begin(); i!=ipUniqueSymbols->end(); ++i)
		if (aDoubleDistribution(aGenerator)<iRandom)
			i->second=(char)(aCharDistribution(aGenerator));
}

void initMap(const std::string& iCipher, std::unordered_map<char, char> *ipUniqueSymbols) {
	std::default_random_engine aGenerator;
	aGenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<unsigned int> aCharDistribution((unsigned int)'a', (unsigned int)'z');

	ipUniqueSymbols->clear();

	std::unordered_set<char>* apSymbols = new std::unordered_set<char>();
	for (unsigned long long i = 0; i < iCipher.length(); i++)
		apSymbols->insert(iCipher[i]);

	for (std::unordered_set<char>::iterator i=apSymbols->begin(); i!=apSymbols->end(); ++i)
		ipUniqueSymbols->insert(std::pair<char, char>(*i, (char)(aCharDistribution(aGenerator))));

	delete apSymbols;
}

void buildClear(const std::string& iCipherString, std::unordered_map<char, char> *ipSymbolMap, std::string* iClear) {
	iClear->clear();
	for (unsigned int i=0; i<iCipherString.length(); i++)
		(*iClear)+=ipSymbolMap->find(iCipherString[i])->second;
}

void insertSymbols(const std::string& iCipherString, std::unordered_map<char, char> *ipSymbolMap, const std::string *ipLetters, const unsigned int iPos) {
	for (unsigned int i=0; i<ipLetters->length(); i++)
		ipSymbolMap->find(iCipherString[iPos+i])->second=ipLetters->at(i);
}

bool checkBest(long double aLoopBestScore, const std::string *ipClear) {
	aBestScoreMux.lock();
	if (aLoopBestScore > aGlobalBestScore) {
		aGlobalBestScore = aLoopBestScore;
		aGlobalBestSolution = std::string(*ipClear);
		aBestScoreMux.unlock();
		return true;
	}
	aBestScoreMux.unlock();
	return false;
}

void computeScoreStatistics(const std::string& iTextFile, std::unordered_map<unsigned long long, NGram*>* iopNorms, const std::string* ipCipherString) {
	std::string* apString = new std::string();
	std::ifstream aFile(iTextFile);
	if (aFile.is_open())
		while (!aFile.eof())
			aFile >> *apString;
	else
		throw "Can not open file";

	aFile.close();
	std::vector<long double> aScores;
	for (std::unordered_map<unsigned long long, NGram*>::iterator i=iopNorms->begin(); i != iopNorms->end(); ++i) {
		aScores.clear();
		std::unordered_map<unsigned long long, NGram*>* apSubNorm = new std::unordered_map<unsigned long long, NGram*>();
		apSubNorm->insert(std::pair<unsigned long long, NGram*>(i->first, i->second));
		for (unsigned long long aPos = 0; aPos < apString->length() - ipCipherString->length()+1; aPos++) {
			std::string aSub = apString->substr(aPos, ipCipherString->length());
			aScores.push_back(score(&aSub, apSubNorm, nullptr));
		}
		delete apSubNorm;
		long double aDouble = 0;
		for (std::vector<long double>::iterator j = aScores.begin(); j != aScores.end(); ++j)
			aDouble += (*j);
		i->second->_mean=aDouble/aScores.size();
		aDouble = 0;
		for (std::vector<long double>::iterator j = aScores.begin(); j != aScores.end(); ++j)
			aDouble += powl((*j) - i->second->_mean, 2);
		i->second->_sigma=sqrtl(aDouble/(aScores.size() - 1));
	}
	delete apString;
}

void log() {
	std::cout << std::endl;
};

template<typename T, typename ... Args>
void log(T first, Args ... args) {
	std::cout << first << " ";
	log(args ...);
}

template<typename ... Args>
void logTime(Args ... args) {
	time_t now_c = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	aOutputMux.lock();
	std::cout << std::put_time(std::localtime(&now_c), "%c") << " ";
	log(args ...);
	aOutputMux.unlock();
}

void hillclimber(const unsigned long long iThread, const std::unordered_map<unsigned long long, NGram*> *ipNorms, const std::string& iCipherString, const long double *ipWorstScore, const std::string &iSeed, const long double iRandom, const long double iMaxIter) {
	std::unordered_map<char, char> *apSymbolMap=new std::unordered_map<char, char>();

	std::default_random_engine aGenerator;
	aGenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<unsigned int> aIntDistribution(0, iCipherString.length());
	std::uniform_real_distribution<long double> aDoubleDistribution(0,1);

	unsigned long long aFails=0;
	long double aTolInit=0.02;
	long double aTolFac=0.95;
	long double aMaxTol=0.005;
	long double aCurTol=aTolInit;

	initMap(iCipherString, apSymbolMap);

	if (iThread==0)
		insertSymbols(iCipherString, apSymbolMap, &iSeed, 0);

	std::ostringstream *iLog=new std::ostringstream();
	std::string* apClear=new std::string();

	long double aLoopBestScore;
	std::string aLoopBestSolution;

	while (true) {
		long double aLastScore=*ipWorstScore;
		bool aLoopImproved;

		buildClear(iCipherString, apSymbolMap, apClear);
		aLoopBestScore=score(apClear, ipNorms, iLog);

		logTime("DEBUG Thread:", iThread, "Restart", "Score:", aLoopBestScore, *apClear);

		do {
			aLoopImproved=false;

			std::unordered_map<std::string, unsigned long long>* apTestMap=ipNorms->find(1)->second->_NGramMap;
			unsigned int aOffset=aIntDistribution(aGenerator);
			for (unsigned int aCounter=0; aCounter<iCipherString.length(); aCounter++) {
				unsigned int aPos=(aCounter+aOffset)%iCipherString.length();
				std::string aBestChoiceSoFar(1, apSymbolMap->find(iCipherString[aPos])->second);

				for (std::unordered_map<std::string, unsigned long long>::const_iterator aTestNGram=apTestMap->begin(); aTestNGram!=apTestMap->end(); ++aTestNGram) {
					insertSymbols(iCipherString, apSymbolMap, &aTestNGram->first, aPos);
					buildClear(iCipherString, apSymbolMap, apClear);
					long double aCurrentScore=score(apClear, ipNorms, iLog);

					long double aTol=aCurTol*aDoubleDistribution(aGenerator);

					if (aCurrentScore*(1-aTol)>aLastScore) {
						aLastScore=aCurrentScore;
						aBestChoiceSoFar=aTestNGram->first;

						if (aCurrentScore>aLoopBestScore) {
							aLoopBestScore=aCurrentScore;
							aLoopBestSolution=*apClear;
							aLoopImproved=true;
							aFails=0;

							if (checkBest(aCurrentScore, apClear))
								logTime("Thread:", iThread, "Score:", aLoopBestScore, iLog->str(), aCurTol, *apClear);
						}
					}
				}
				insertSymbols(iCipherString, apSymbolMap, &aBestChoiceSoFar, aPos);
			}
		} while (aLoopImproved);

		buildClear(iCipherString, apSymbolMap, apClear);
		aLoopBestScore=score(apClear, ipNorms, iLog);
		logTime("DEBUG Thread:", iThread, "Give Up", "Score:", aLoopBestScore, *apClear);

		aFails++;

		aCurTol*=aTolFac;
		if (aCurTol<aMaxTol)
			aCurTol=aTolInit;

		if (aFails<iMaxIter && iRandom>0) {
			aBestScoreMux.lock();
			insertSymbols(iCipherString, apSymbolMap, &aLoopBestSolution, 0);
			aBestScoreMux.unlock();
			randomizeMap(apSymbolMap, iRandom);
		} else {
			initMap(iCipherString, apSymbolMap);
			aCurTol=aTolInit;
			aFails=0;
		}
	}

	delete iLog;
	delete apSymbolMap;
	delete apClear;
}

int main( int argc, char* argv[] ) {
	std::list<std::string> aNGramsFiles;
	std::string aCipherFile;
	std::string aTextFile;
	std::string aSeed="";
	unsigned int aMaxIter=0;
	unsigned int aThreadsCount=1;
	long double aRandom=0.0;

	{
		int c;
		while( ( c = getopt(argc, argv, "l:c:t:s:w:r:x:") ) != -1 ) {
			switch(c) {
			case 'l':
				if(optarg)
					aNGramsFiles.push_back(optarg);
				break;
			case 'c':
				if(optarg)
					aCipherFile=optarg;
				break;
			case 'w':
				if(optarg)
					aTextFile=optarg;
				break;
			case 's':
				if(optarg)
					aSeed=optarg;
				break;
			case 't':
				if(optarg)
					aThreadsCount=atoi(optarg);
				break;
			case 'r':
				if(optarg)
					aRandom=std::stold(optarg);
				break;
			case 'x':
				if(optarg)
					aMaxIter=atoi(optarg);
				break;
			}
		}
	}

	std::string *apCipherString=new std::string();
	readCipher(aCipherFile, apCipherString);

	std::cout << "Cipher length " << apCipherString->length() << std::endl;
	std::cout << *apCipherString << std::endl;

	std::cout << "Randomize fraction: " << aRandom << std::endl;
	std::cout << "Random re-initialization after " << aMaxIter << " iterations" << std::endl;

	std::unordered_map<unsigned long long, NGram*> *apNorms=new std::unordered_map<unsigned long long, NGram*>();
	readNorms(aNGramsFiles, apNorms);

	computeScoreStatistics(aTextFile, apNorms, apCipherString);

	{
		long double aLnPerfect=0.0;
		for (std::unordered_map<unsigned long long, NGram*>::iterator i=apNorms->begin(); i!=apNorms->end(); ++i) {
			long double aLnNGramPerfect=-logl(sqrtl(2.0*M_PI) * i->second->_sigma);
			std::cout << std::setprecision(6) << "NGram length:" << i->second->_length << " NGrams:" << i->second->_NGramMap->size() << " Samples:" << i->second->_count << " Mean:" << i->second->_mean << " StdDev:" << i->second->_sigma << " Perfect: " << aLnNGramPerfect << std::endl;
			aLnPerfect += aLnNGramPerfect;
		}
		std::cout << "Maximum reachable score " << aLnPerfect << std::endl;
	}

	if (aSeed.length()>0) {
		std::ostringstream *iLog=new std::ostringstream();
		const long double aSeedScore=score(&aSeed, apNorms, iLog);
		std::cout << "Seed: " << aSeedScore << " " << iLog->str() << " " << aSeed << std::endl;
		delete iLog;
	}

	std::string *apString=new std::string(apCipherString->length(), '.');
	const long double *apWorstScore=new long double(score(apString, apNorms, NULL));
	aGlobalBestScore=*apWorstScore;
	delete apString;

	std::vector<std::thread> aThreads[aThreadsCount];
	for (unsigned long long aThread=0; aThread<aThreadsCount; aThread++)
		aThreads->push_back(std::thread(&hillclimber, aThread, apNorms, *apCipherString, apWorstScore, aSeed, aRandom, aMaxIter));

	logTime(aThreadsCount, "threads started.");

	for (std::vector<std::thread>::iterator i=aThreads->begin(); i!=aThreads->end(); ++i)
		(*i).join();

	delete apNorms;
	delete apCipherString;
	delete apWorstScore;

	return EXIT_SUCCESS;
}
