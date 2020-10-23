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

#include "NGram.h"
#include "Lock.h"
#include "Score.h"

//Globals
std::mutex aBestScoreMutex;
std::mutex aOutputMutex;
Score aGlobalBestScore;
std::string aGlobalBestSolution;
bool aVerbose=false;

void readCipher(const std::string& aCipherfile, std::string& iCipher) {
	std::ifstream myfile(aCipherfile);
	if (myfile.is_open()) {
		while (!myfile.eof()) {
			std::string aLine;
			myfile >> aLine;
			iCipher += aLine;
		}
		myfile.close();
	}
}

void readNorms(std::list<std::string> afiles, std::unordered_map<unsigned long long, NGram*>& aNorms) {
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
					std::unordered_map<unsigned long long, NGram*>::iterator i=aNorms.find(aString.length());
					NGram* aNorm;
					if (i == aNorms.end()) {
						aNorm = new NGram(aString.length());
						aNorms.insert(std::pair<unsigned long long, NGram*>(aString.length(), aNorm));
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

void partiallyShuffleMap(std::unordered_map<char, char>& iUniqueSymbols, const long double &iRandom) {
	std::default_random_engine aGenerator;
	aGenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<unsigned int> aCharDistribution((unsigned int)'a', (unsigned int)'z');
	std::uniform_real_distribution<long double> aDoubleDistribution(0,1);

	for (std::unordered_map<char, char>::iterator i=iUniqueSymbols.begin(); i!=iUniqueSymbols.end(); ++i)
		if (aDoubleDistribution(aGenerator)<iRandom)
			i->second=(char)(aCharDistribution(aGenerator));
}

void randomMapInit(const std::string& iCipher, std::unordered_map<char, char>& iSymbolMap) {
	std::default_random_engine aGenerator;
	aGenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<unsigned int> aCharDistribution((unsigned int)'a', (unsigned int)'z');

	iSymbolMap.clear();

	std::unordered_set<char> aSymbols;
	for (unsigned long long i = 0; i < iCipher.length(); i++)
		aSymbols.insert(iCipher[i]);

	for (std::unordered_set<char>::iterator i=aSymbols.begin(); i!=aSymbols.end(); ++i)
		iSymbolMap.insert(std::pair<char, char>(*i, (char)(aCharDistribution(aGenerator))));
}

void buildClear(const std::string& iCipherString, std::unordered_map<char, char>& iSymbolMap, std::string& iClear) {
	iClear.clear();
	for (unsigned int i=0; i<iCipherString.length(); i++)
		iClear+=iSymbolMap.find(iCipherString[i])->second;
}

void insertSymbols(const std::string& iCipherString, std::unordered_map<char, char>& iSymbolMap, const std::string& iLetters, const unsigned int iPos) {
	for (unsigned int i=0; i<iLetters.length(); i++)
		iSymbolMap.find(iCipherString[iPos+i])->second=iLetters.at(i);
}

bool checkIfGlobalBest(Score& aLoopBestScore, const std::string& iClear) {
	Lock aLock(aBestScoreMutex);
	if (aLoopBestScore > aGlobalBestScore) {
		aGlobalBestScore = aLoopBestScore;
		aGlobalBestSolution = std::string(iClear);
		return true;
	}
	return false;
}

void computeScoreStatistics(const std::string& iTextFile, std::unordered_map<unsigned long long, NGram*>& ioNorms, const std::string& iCipherString) {
	std::string* apString = new std::string();
	std::ifstream aFile(iTextFile);
	if (aFile.is_open())
		while (!aFile.eof())
			aFile >> *apString;
	else
		throw "Can not open file";

	aFile.close();
	std::vector<Score> aScores;
	for (std::unordered_map<unsigned long long, NGram*>::iterator i=ioNorms.begin(); i != ioNorms.end(); ++i) {
		aScores.clear();
		std::unordered_map<unsigned long long, NGram*> aSubNorm;
		aSubNorm.insert(std::pair<unsigned long long, NGram*>(i->first, i->second));
		for (unsigned long long aPos = 0; aPos < apString->length() - iCipherString.length()+1; aPos++) {
			std::string aSub = apString->substr(aPos, iCipherString.length());
			aScores.push_back(Score(aSubNorm, aSub));
		}
		long double aDouble = 0;
		for (std::vector<Score>::iterator j = aScores.begin(); j != aScores.end(); ++j)
			aDouble += j->rate();
		i->second->_mean=aDouble/aScores.size();
		aDouble = 0;
		for (std::vector<Score>::iterator j = aScores.begin(); j != aScores.end(); ++j)
			aDouble += powl(j->rate() - i->second->_mean, 2);
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

void hillclimber(const unsigned long long& iThread, const std::unordered_map<unsigned long long, NGram*>& iNorms, const std::string& iCipherString, const std::string &iSeedString, const long double& iRandomFraction, const long double& iMaxIter, const long double& iFuzzy) {

	std::unordered_map<char, char> aSymbolMap;

	std::default_random_engine aGenerator;
	aGenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<unsigned int> aIntDistribution(0, iCipherString.length());
	std::uniform_real_distribution<long double> aDoubleDistribution(0, 1.0);

	unsigned long long aConsecutiveFailuresToImprove=0;
	long double aCurrentTolerance=0.02;

	randomMapInit(iCipherString, aSymbolMap);

	if (iThread==0)
		insertSymbols(iCipherString, aSymbolMap, iSeedString, 0);

	std::string aClear;

	Score aClimberBestScore(iNorms, aClear);
	std::string aClimberBestSolution=aClear;

	while (true) {
		bool aLoopImproved;

		buildClear(iCipherString, aSymbolMap, aClear);
		Score aLoopBestScore(iNorms, aClear);

		if (checkIfGlobalBest(aLoopBestScore, aClear))
			logTime("Thread:", iThread, "Score:", aLoopBestScore, "Tolerance:", aCurrentTolerance, aClear);

		if (aVerbose)
			logTime("DEBUG Thread:", iThread, "Restart", "Tolerance:", aCurrentTolerance, "Score:", aLoopBestScore, aClear);

		do {
			aLoopImproved=false;
			Score aLastScore(aLoopBestScore);
			unsigned int aTolerated=0;

			std::unordered_map<std::string, unsigned long long> aAlphabet=iNorms.find(1)->second->_NGramMap;
			for (std::unordered_map<char, char>::iterator aMappedSymbol=aSymbolMap.begin(); aMappedSymbol!=aSymbolMap.end(); ++aMappedSymbol) {
				const char aBefore=aMappedSymbol->second;
				char aBestChoiceSoFar=aBefore;
				for (std::unordered_map<std::string, unsigned long long>::const_iterator aMappedLetter=aAlphabet.begin(); aMappedLetter!=aAlphabet.end(); ++aMappedLetter) {
					if (aMappedLetter->first.at(0)!=aBefore) {
						aMappedSymbol->second=aMappedLetter->first.at(0);
						buildClear(iCipherString, aSymbolMap, aClear);

						Score aCurrentScore(iNorms, aClear);
						long double aTolerance=aCurrentTolerance*aDoubleDistribution(aGenerator);

						if (aCurrentScore.rate()*(1.0-aTolerance)>aLastScore.rate()) {
							if (aCurrentScore<aLastScore)
								aTolerated++;

							aLastScore=aCurrentScore;
							aBestChoiceSoFar=aMappedLetter->first.at(0);

							if (aCurrentScore>aLoopBestScore) {
								aLoopBestScore=aCurrentScore;
								aLoopImproved=true;

								if (aCurrentScore>aClimberBestScore) {
									aClimberBestScore=aCurrentScore;
									aClimberBestSolution=aClear;
									aConsecutiveFailuresToImprove=0;

									if (checkIfGlobalBest(aCurrentScore, aClear))
										logTime("Thread:", iThread, "Score:", aLoopBestScore, "Tolerance:", aCurrentTolerance, aClear);
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
			buildClear(iCipherString, aSymbolMap, aClear);
			aLoopBestScore=Score(iNorms, aClear);

			logTime("DEBUG Thread:", iThread, "Give Up", "Tolerance:", aCurrentTolerance, "Score:", aLoopBestScore, aClear);
		}
		aConsecutiveFailuresToImprove++;

		if (aConsecutiveFailuresToImprove<iMaxIter && iRandomFraction>0) {
			insertSymbols(iCipherString, aSymbolMap, aClimberBestSolution, 0);
			partiallyShuffleMap(aSymbolMap, iRandomFraction);
		} else {
			randomMapInit(iCipherString, aSymbolMap);
			aConsecutiveFailuresToImprove=0;
		}
	}
}

long double computePerfectScore(std::unordered_map<unsigned long long, NGram*>& aNorms) {
	long double aLnPerfect = 0.0;
	for (std::unordered_map<unsigned long long, NGram*>::iterator i =
			aNorms.begin(); i != aNorms.end(); ++i) {
		long double aLnNGramPerfect = -logl(
				sqrtl(2.0 * M_PI) * i->second->_sigma);
		std::cout << setiosflags(std::ios::fixed) << std::setprecision(6)
		<< "NGram length:" << i->second->_length << " NGrams:"
		<< i->second->_NGramMap.size() << " Samples:"
		<< i->second->_count << " Mean:" << i->second->_mean
		<< " StdDev:" << i->second->_sigma << " Perfect: "
		<< aLnNGramPerfect << std::endl;
		aLnPerfect += aLnNGramPerfect;
	}
	return aLnPerfect;
}

int main( int argc, char* argv[] ) {
	std::list<std::string> aNGramsFiles;
	std::string aCipherFile;
	std::string aTextFile;
	std::string aSeed="";
	unsigned int aMaxIter=250;
	unsigned int aThreadsCount=1;
	long double aRandom=0.2;
	long double aFuzzy=0.05;

	std::cout << "cDecryptor Version 23.10.2020 17:17" << std::endl;

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

	std::string aCipherString;
	readCipher(aCipherFile, aCipherString);

	std::cout << "Cipher: " << aCipherString << std::endl;
	std::cout << "Cipher length: " << aCipherString.length() << std::endl;
	std::cout << "Randomize fraction: " << aRandom << std::endl;
	std::cout << "Random re-initialization after " << aMaxIter << " iterations" << std::endl;
	std::cout << "Tolerance factor: " << aFuzzy << std::endl;
	std::cout << "Parallel threads: " << aThreadsCount << std::endl;

	std::unordered_map<unsigned long long, NGram*> aNorms;
	readNorms(aNGramsFiles, aNorms);

	computeScoreStatistics(aTextFile, aNorms, aCipherString);
	std::cout << computePerfectScore(aNorms) << std::endl;

	if (aSeed.length()>0)
		std::cout << "Seed: " << Score(aNorms, aSeed) << " " << aSeed << std::endl;

	aGlobalBestScore=Score(aNorms, std::string(aCipherString.length(), '.'));

	std::vector<std::thread> aThreads[aThreadsCount];
	for (unsigned long long aThread=0; aThread<aThreadsCount; aThread++)
		aThreads->push_back(std::thread(&hillclimber, aThread, aNorms, aCipherString, aSeed, aRandom, aMaxIter, aFuzzy));

	logTime(aThreadsCount, "threads started.");

	for (std::vector<std::thread>::iterator i=aThreads->begin(); i!=aThreads->end(); ++i)
		(*i).join();

	return EXIT_SUCCESS;
}
