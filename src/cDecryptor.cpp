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
#include "Lock.h"

//Globals
std::mutex aBestScoreMutex;
std::mutex aOutputMutex;
long double aGlobalBestScore;
std::string aGlobalBestSolution;
bool aVerbose=false;

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
	bool aFirstLog=true;

	if (iLog!=NULL)
		iLog->str("");

	std::unordered_map<std::string, unsigned long long> *aObserved=new std::unordered_map<std::string, unsigned long long>();

	for (std::unordered_map<unsigned long long, NGram*>::const_iterator aNorm=ipNorms->cbegin(); aNorm!=ipNorms->cend(); ++aNorm) {
		long double aLoopScore=0;
		const unsigned int aLength=aNorm->second->_length;

		for (unsigned int i=0; i+aLength<=ipCandidate->length(); i++) {
			const std::string aSub=ipCandidate->substr(i, aLength);
			std::unordered_map<std::string, unsigned long long>::iterator b=aObserved->find(aSub);
			if (b!=aObserved->end())
				b->second++;
			else
				aObserved->insert(std::pair<std::string, unsigned long long>(aSub, 1));
		}

		const long double aNeverSeen=logl(2)/(long double)aNorm->second->_count;
		for (std::unordered_map<std::string, unsigned long long>::const_iterator b=aObserved->cbegin(); b!=aObserved->cend(); ++b) {
			long double aP;
			std::unordered_map<std::string, unsigned long long>::const_iterator j=aNorm->second->_NGramMap->find(b->first);
			if (j!=aNorm->second->_NGramMap->end())
				aP=(long double)j->second/(long double)aNorm->second->_count;
			else
				aP=aNeverSeen;

			aLoopScore+=b->second*logl(aP)-logFac(b->second);
		}

		if (aNorm->second->_mean<0) {
			const long double aLnGaussScore=lnGauss(aLoopScore, aNorm->second->_mean, aNorm->second->_sigma);
			aScore+=aLnGaussScore;
			if (iLog!=NULL) {
				if (!aFirstLog) {
					(*iLog) << " ";
				} else
					aFirstLog=false;
				(*iLog) << aLength << ":" << aLnGaussScore;
			}
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
	Lock aLock(aBestScoreMutex);
	if (aLoopBestScore > aGlobalBestScore) {
		aGlobalBestScore = aLoopBestScore;
		aGlobalBestSolution = std::string(*ipClear);
		return true;
	}
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
	Lock aLock(aOutputMutex);
	std::cout << std::put_time(std::localtime(&now_c), "%c") << " ";
	log(args ...);
}

void hillclimber(const unsigned long long iThread, const std::unordered_map<unsigned long long, NGram*> *ipNorms, const std::string& iCipherString, const std::string &iSeed, const long double iRandom, const long double iMaxIter, const long double iFuzzy) {

	std::unordered_map<char, char> *apSymbolMap=new std::unordered_map<char, char>();

	std::default_random_engine aGenerator;
	aGenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<unsigned int> aIntDistribution(0, iCipherString.length());
	std::uniform_real_distribution<long double> aDoubleDistribution(0, 1.0);

	unsigned long long aFails=0;
	long double aCurrentTolerance=0.02;

	initMap(iCipherString, apSymbolMap);

	if (iThread==0)
		insertSymbols(iCipherString, apSymbolMap, &iSeed, 0);

	std::ostringstream *iLog=new std::ostringstream();
	std::string* apClear=new std::string();

	buildClear(iCipherString, apSymbolMap, apClear);
	long double aClimberBestScore=score(apClear, ipNorms, iLog);
	std::string aClimberBestSolution=*apClear;

	while (true) {
		bool aLoopImproved;

		buildClear(iCipherString, apSymbolMap, apClear);
		long double aLoopBestScore=score(apClear, ipNorms, iLog);

		if (checkBest(aLoopBestScore, apClear))
			logTime("Thread:", iThread, "Score:", aLoopBestScore, iLog->str(), "Tolerance:", aCurrentTolerance, *apClear);

		if (aVerbose)
			logTime("DEBUG Thread:", iThread, "Restart", "Tolerance:", aCurrentTolerance, "Score:", aLoopBestScore, *apClear);

		do {
			aLoopImproved=false;
			long double aLastScore=aLoopBestScore;
			unsigned int aTolerated=0;

			std::unordered_map<std::string, unsigned long long>* apAlphabet=ipNorms->find(1)->second->_NGramMap;
			for (std::unordered_map<char, char>::iterator aMappedSymbol=apSymbolMap->begin(); aMappedSymbol!=apSymbolMap->end(); ++aMappedSymbol) {
				const char aBefore=aMappedSymbol->second;
				char aBestChoiceSoFar=aBefore;
				for (std::unordered_map<std::string, unsigned long long>::const_iterator aMappedLetter=apAlphabet->begin(); aMappedLetter!=apAlphabet->end(); ++aMappedLetter) {
					if (aMappedLetter->first.at(0)!=aBefore) {
						aMappedSymbol->second=aMappedLetter->first.at(0);
						buildClear(iCipherString, apSymbolMap, apClear);

						long double aCurrentScore=score(apClear, ipNorms, iLog);
						long double aTolerance=aCurrentTolerance*aDoubleDistribution(aGenerator);

						if (aCurrentScore*(1.0-aTolerance)>aLastScore) {
							if (aCurrentScore<aLastScore)
								aTolerated++;

							aLastScore=aCurrentScore;
							aBestChoiceSoFar=aMappedLetter->first.at(0);

							if (aCurrentScore>aLoopBestScore) {
								aLoopBestScore=aCurrentScore;
								aLoopImproved=true;

								if (aCurrentScore>aClimberBestScore) {
									aClimberBestScore=aCurrentScore;
									aClimberBestSolution=*apClear;
									aFails=0;

									if (checkBest(aCurrentScore, apClear))
										logTime("Thread:", iThread, "Score:", aLoopBestScore, iLog->str(), "Tolerance:", aCurrentTolerance, *apClear);
								}
							}
						}
					}
				}
				aMappedSymbol->second=aBestChoiceSoFar;
			}

			if (aTolerated>iFuzzy*iCipherString.length())
				aCurrentTolerance*=0.95;
			else {
				aCurrentTolerance*=1.05;
				if (aCurrentTolerance>1)
					aCurrentTolerance=1;
			}
		} while (aLoopImproved);

		if (aVerbose) {
			buildClear(iCipherString, apSymbolMap, apClear);
			aLoopBestScore=score(apClear, ipNorms, iLog);

			logTime("DEBUG Thread:", iThread, "Give Up", "Tolerance:", aCurrentTolerance, "Score:", aLoopBestScore, *apClear);
		}
		aFails++;

		if (aFails<iMaxIter && iRandom>0) {
			insertSymbols(iCipherString, apSymbolMap, &aClimberBestSolution, 0);
			randomizeMap(apSymbolMap, iRandom);
		} else {
			initMap(iCipherString, apSymbolMap);
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
	unsigned int aMaxIter=250;
	unsigned int aThreadsCount=1;
	long double aRandom=0.0;
	long double aFuzzy=0.05;

	std::cout << "cDecryptor Version 15.10.2020 17:59" << std::endl;

	{
		int c;
		while( ( c = getopt(argc, argv, "l:c:t:f:s:w:r:x:v") ) != -1 ) {
			switch(c) {
			case 'c':
				if(optarg)
					aCipherFile=optarg;
				break;
			case 'f':
				if(optarg)
					aFuzzy=std::stold(optarg);
				break;
			case 'l':
				if(optarg)
					aNGramsFiles.push_back(optarg);
				break;
			case 'r':
				if(optarg)
					aRandom=std::stold(optarg);
				break;
			case 's':
				if(optarg)
					aSeed=optarg;
				break;
			case 't':
				if(optarg)
					aThreadsCount=atoi(optarg);
				break;
			case 'v':
				aVerbose=true;
				break;
			case 'w':
				if(optarg)
					aTextFile=optarg;
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
	std::cout << "Tolerance factor: " << aFuzzy << std::endl;
	std::cout << "Parallel threads: " << aThreadsCount << std::endl;

	std::unordered_map<unsigned long long, NGram*> *apNorms=new std::unordered_map<unsigned long long, NGram*>();
	readNorms(aNGramsFiles, apNorms);

	computeScoreStatistics(aTextFile, apNorms, apCipherString);

	{
		long double aLnPerfect=0.0;
		for (std::unordered_map<unsigned long long, NGram*>::iterator i=apNorms->begin(); i!=apNorms->end(); ++i) {
			long double aLnNGramPerfect=-logl(sqrtl(2.0*M_PI) * i->second->_sigma);
			std::cout << setiosflags(std::ios::fixed) << std::setprecision(6) << "NGram length:" << i->second->_length << " NGrams:" << i->second->_NGramMap->size() << " Samples:" << i->second->_count << " Mean:" << i->second->_mean << " StdDev:" << i->second->_sigma << " Perfect: " << aLnNGramPerfect << std::endl;
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
		aThreads->push_back(std::thread(&hillclimber, aThread, apNorms, *apCipherString, aSeed, aRandom, aMaxIter, aFuzzy));

	logTime(aThreadsCount, "threads started.");

	for (std::vector<std::thread>::iterator i=aThreads->begin(); i!=aThreads->end(); ++i)
		(*i).join();

	delete apNorms;
	delete apCipherString;
	delete apWorstScore;

	return EXIT_SUCCESS;
}
