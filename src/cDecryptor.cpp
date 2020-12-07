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
#include <csignal>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

#include "NGram.h"
#include "Lock.h"
#include "RatedScore.h"
#include "Options.h"

//Globals
std::mutex aGlobalBestScoreMutex;
std::mutex aGlobalOutputMutex;
RatedScore aGlobalBestScore;
std::string aGlobalBestSolution;
std::unordered_map<unsigned long long, GaussianNorm> aGlobalScoreStatistics;

void readCipher(const std::string& iCipherFilename, std::string& iCipher) {
	std::ifstream aFile(iCipherFilename);
	if (aFile.is_open()) {
		while (!aFile.eof()) {
			std::string aLine;
			aFile >> aLine;
			iCipher += aLine;
		}
		aFile.close();
	} else
		throw "Failed to open "+iCipherFilename;
}

void readNorms(std::list<std::string>& iFileNames, std::unordered_map<unsigned long long, NGram*>& iNorms) {
	for (std::list<std::string>::iterator i = iFileNames.begin(); i != iFileNames.end(); ++i) {
		std::cout << "Reading norm file " << *i << std::endl;
		std::ifstream aFile(*i);
		if (aFile.is_open()) {
			while (!aFile.eof()) {
				std::string aString;
				aFile >> aString;
				transform(aString.begin(), aString.end(), aString.begin(), ::tolower);
				unsigned long long aCount;
				aFile >> aCount;
				if (aString.length() > 0 && aCount > 0) {
					std::unordered_map<unsigned long long, NGram*>::iterator i=iNorms.find(aString.length());
					NGram* aNorm;
					if (i == iNorms.end()) {
						aNorm = new NGram(aString.length());
						iNorms.insert(std::pair<unsigned long long, NGram*>(aString.length(), aNorm));
					} else
						aNorm = (*i).second;

					aNorm->add(aCount, aString);
				}
			}
			aFile.close();
		} else
			std::cout << "Failed to open file" << (*i) << std::endl;
	}
}

void partiallyShuffleMap(std::unordered_map<char, unsigned int>& iSymbolMap, const long double &iRandom) {
	std::default_random_engine aGenerator;
	aGenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<unsigned int> aCharDistribution(0, (unsigned int)'z'-(unsigned int)'a');
	std::uniform_real_distribution<long double> aDoubleDistribution(0,1);

	for (std::unordered_map<char, unsigned int>::iterator i=iSymbolMap.begin(); i!=iSymbolMap.end(); ++i)
		if (aDoubleDistribution(aGenerator)<iRandom)
			i->second=(char)(aCharDistribution(aGenerator));
}

void randomMapInit(const std::string& iCipher, std::unordered_map<char, unsigned int>& iSymbolMap, std::vector<char>& iLetterVec) {
	std::default_random_engine aGenerator;
	aGenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<unsigned int> aIntDistribution(0, (unsigned int)'z'-(unsigned int)'a');

	iSymbolMap.clear();
	iLetterVec.clear();

	std::unordered_set<char> aSymbols;
	for (unsigned long long i=0; i<iCipher.length(); i++)
		aSymbols.insert(iCipher[i]);

	for (std::unordered_set<char>::iterator i=aSymbols.begin(); i!=aSymbols.end(); ++i)
		iSymbolMap.insert(std::pair<char, unsigned int>(*i, aIntDistribution(aGenerator)));

	for (char c='a'; c<='z'; c++)
		iLetterVec.push_back(c);
}

std::string buildClear(const std::string& iCipherString, std::unordered_map<char, unsigned int>& iSymbolMap, std::vector<char>& iLetterVec, const Options& iOptions) {
	std::string iClear;
	for (unsigned int i=0; i<iCipherString.length(); i++)
		iClear+=iLetterVec.at(iSymbolMap.find(iCipherString[i])->second % (iOptions._diskSize>0?iOptions._diskSize:iLetterVec.size()));
	return iClear;
}

void insertSymbols(const std::string& iCipherString, std::unordered_map<char, unsigned int>& iSymbolMap, const std::vector<char>& iLetterVec, const std::string& iLetters, const unsigned int iPos) {
	for (unsigned int i=0; i<iLetters.length(); i++) {
		int j=iLetterVec.size()-1;
		while (iLetterVec.at(j)!=iLetters.at(i) && j>=0)
			j--;
		if (j>=0)
			iSymbolMap.find(iCipherString[iPos+i])->second=(unsigned int)j;
	}
}

bool checkIfGlobalBest(const RatedScore& iRatedScore, const std::string& iClear) {
	Lock aLock(aGlobalBestScoreMutex);
	if (iRatedScore>aGlobalBestScore) {
		aGlobalBestScore=iRatedScore;
		aGlobalBestSolution=std::string(iClear);
		return true;
	}
	return false;
}

void computeScoreStatistics(const std::string& iTextFile, std::unordered_map<unsigned long long, NGram*>& ioNorms, const std::string& iCipherString) {
	std::string aString;
	std::ifstream aFile(iTextFile);

	if (aFile.is_open())
		while (!aFile.eof())
			aFile >> aString;
	else
		throw "Can not open file " + iTextFile;

	aFile.close();
	std::vector<Score> aScores;
	aGlobalScoreStatistics.clear();

	for (std::unordered_map<unsigned long long, NGram*>::iterator i=ioNorms.begin(); i != ioNorms.end(); ++i) {
		aScores.clear();
		std::unordered_map<unsigned long long, NGram*> aSubNorm;
		aSubNorm.insert(std::pair<unsigned long long, NGram*>(i->first, i->second));
		for (unsigned long long aPos = 0; aPos < aString.length() - iCipherString.length()+1; aPos++) {
			std::string aSub = aString.substr(aPos, iCipherString.length());
			aScores.push_back(Score(aSubNorm, aSub));
		}

		long double aMean = 0;
		for (std::vector<Score>::iterator j = aScores.begin(); j != aScores.end(); ++j)
			aMean += j->values().find(i->first)->second;
		aMean/=aScores.size();
		long double aSigma = 0;
		for (std::vector<Score>::iterator j = aScores.begin(); j != aScores.end(); ++j)
			aSigma += powl(j->values().find(i->first)->second - aMean, 2);
		aSigma=sqrtl(aSigma/(aScores.size() - 1));

		aGlobalScoreStatistics.insert(std::pair<unsigned long long , GaussianNorm>(i->first,GaussianNorm(aMean,aSigma)));
	}
}

void log() {
	std::cout << std::endl;
};

template<typename T, typename ... Args>
void log(T first, Args ... iArgs) {
	std::cout << first << " ";
	log(iArgs ...);
}

template<typename ... Args>
void logTime(Args ... iArgs) {
	time_t aNow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
	Lock aLock(aGlobalOutputMutex);
	std::cout << std::put_time(std::localtime(&aNow), "%c") << " ";
	log(iArgs ...);
}

bool checkIfLocalBest(const RatedScore& iCandidateScore, const std::string& iCandidateSolution,
		RatedScore& ioPreviousBestScore, std::string& oPreviousBestSolution,
		const unsigned long long& iThread, long double& iCurrentTolerance) {

	if (iCandidateScore>ioPreviousBestScore) {
		ioPreviousBestScore = iCandidateScore;
		oPreviousBestSolution = iCandidateSolution;

		if (checkIfGlobalBest(iCandidateScore, iCandidateSolution))
			logTime("Thread:", iThread, "Score:", iCandidateScore, "Tolerance:", iCurrentTolerance, iCandidateSolution);

		return true;
	} else
		return false;
}

void hillclimber(const unsigned long long& iThread, const std::unordered_map<unsigned long long, NGram*>& iNorms, const std::string& iCipherString, const Options& iOptions) {
	std::unordered_map<char, unsigned int> aSymbolMap;
	std::vector<char> aLetterVec;
	std::default_random_engine aGenerator;
	aGenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
	std::uniform_int_distribution<unsigned int> aIntDistribution(0, iCipherString.length());
	std::uniform_real_distribution<long double> aDoubleDistribution(0, 1.0);
	unsigned long long aCounterUntilReset=iOptions._maxiter;
	long double aCurrentTolerance=0.02;

	randomMapInit(iCipherString, aSymbolMap, aLetterVec);

	if (iThread==0)
		insertSymbols(iCipherString, aSymbolMap, aLetterVec, iOptions._seed, 0);

	while (true) {
		std::string aClimberBestSolution=buildClear(iCipherString, aSymbolMap, aLetterVec, iOptions);
		RatedScore aClimberBestScore(Score(iNorms, aClimberBestSolution), aGlobalScoreStatistics);

		if (checkIfGlobalBest(aClimberBestScore, aClimberBestSolution))
			logTime("Thread:", iThread, "Score:", aClimberBestScore, "Tolerance:", aCurrentTolerance, aClimberBestSolution);

		while (aCounterUntilReset) {
			RatedScore aLoopBestScore;
			{
				std::string aClear=buildClear(iCipherString, aSymbolMap, aLetterVec, iOptions);
				aLoopBestScore=RatedScore(Score(iNorms, aClear), aGlobalScoreStatistics);

				if (iOptions._verbose)
					logTime("DEBUG Thread:", iThread, "Restart", "Tolerance:", aCurrentTolerance, "Score:", aLoopBestScore, aClear);

				if (checkIfLocalBest(aLoopBestScore, aClear, aClimberBestScore, aClimberBestSolution, iThread, aCurrentTolerance))
					aCounterUntilReset=iOptions._maxiter;
			}

			bool aLoopImproved;
			do {
				aLoopImproved=false;
				RatedScore aLastScore=aLoopBestScore;
				unsigned int aTolerated=0;

				for (std::unordered_map<char, unsigned int>::iterator aMappedSymbol=aSymbolMap.begin(); aMappedSymbol!=aSymbolMap.end(); ++aMappedSymbol) {
					const unsigned int aBefore=aMappedSymbol->second;
					char aBestChoiceSoFar=aBefore;
					for (unsigned int i=0; i<aLetterVec.size(); i++) {
						if (i!=aBefore) {
							aMappedSymbol->second=i;
							std::string aCandidateSolution=buildClear(iCipherString, aSymbolMap, aLetterVec, iOptions);

							RatedScore aCandidateScore(Score(iNorms, aCandidateSolution), aGlobalScoreStatistics);
							long double aTolerance=aCurrentTolerance*aDoubleDistribution(aGenerator);

							if (aCandidateScore.value()*(1.0-aTolerance)>aLastScore.value()) {
								if (aCandidateScore<aLastScore)
									aTolerated++;

								aLastScore=aCandidateScore;
								aBestChoiceSoFar=i;

								if (aCandidateScore>aLoopBestScore) {
									aLoopBestScore=aCandidateScore;
									aLoopImproved=true;

									if (checkIfLocalBest(aLoopBestScore, aCandidateSolution, aClimberBestScore, aClimberBestSolution, iThread, aCurrentTolerance))
										aCounterUntilReset=iOptions._maxiter;
								}
							}
						}
					}
					aMappedSymbol->second=aBestChoiceSoFar;
				}

				if (aTolerated>iOptions._fuzzy*iCipherString.length())
					aCurrentTolerance*=0.95;
				else {
					aCurrentTolerance*=1.05;
					if (aCurrentTolerance>1)
						aCurrentTolerance=1;
				}
			} while (aLoopImproved);

			if (iOptions._verbose) {
				std::string aClear=buildClear(iCipherString, aSymbolMap, aLetterVec, iOptions);
				aLoopBestScore=RatedScore(Score(iNorms, aClear), aGlobalScoreStatistics);

				logTime("DEBUG Thread:", iThread, "Give Up", "Tolerance:", aCurrentTolerance, "Score:", aLoopBestScore, aClear);
			}

			if (iOptions._random>0) {
				insertSymbols(iCipherString, aSymbolMap, aLetterVec, aClimberBestSolution, 0);
				partiallyShuffleMap(aSymbolMap, iOptions._random);
			}

			aCounterUntilReset--;
		}
		randomMapInit(iCipherString, aSymbolMap, aLetterVec);
		aCounterUntilReset=iOptions._maxiter;
	}
}

void printBestPossibleScore(std::unordered_map<unsigned long long, NGram*>& iNorms) {
	long double aLnPerfect = 0.0;
	for (std::unordered_map<unsigned long long, NGram*>::iterator i=iNorms.begin(); i != iNorms.end(); ++i) {
		long double aLnNGramPerfect = -logl(sqrtl(2.0 * M_PI) * aGlobalScoreStatistics.at(i->first)._sigma);
		std::cout << setiosflags(std::ios::fixed) << std::setprecision(6)
				<< "NGram length:" << i->second->_length << " NGrams:"
				<< i->second->_NGramMap.size() << " Samples:"
				<< i->second->_count << " Mean:" << aGlobalScoreStatistics.at(i->first)._mean
				<< " StdDev:" << aGlobalScoreStatistics.at(i->first)._sigma << " Perfect: "
				<< aLnNGramPerfect << std::endl;
		aLnPerfect += aLnNGramPerfect;
	}
	std::cout << "Best possible score: " << aLnPerfect << std::endl;
}

void signalHandler(const int iSigNum) {
	std::cout << "Interrupt signal " << iSigNum << " received. Exiting." << std::endl;
	exit(iSigNum);
}

void parseOptions(const int iArgc, char* iArgv[], Options& oOptions) {
	int aInt;
	while ((aInt = getopt(iArgc, iArgv, "c:d:f:l:r:s:t:vw:x:")) != -1) {
		switch (aInt) {
		case 'c':
			if (optarg)
				oOptions._cipherfile = optarg;
			break;
		case 'd':
			if (optarg)
				oOptions._diskSize = atoi(optarg);
			break;
		case 'f':
			if (optarg)
				oOptions._fuzzy = std::stold(optarg);
			break;
		case 'l':
			if (optarg)
				oOptions._ngramsfiles.push_back(optarg);
			break;
		case 'r':
			if (optarg)
				oOptions._random = std::stold(optarg);
			break;
		case 's':
			if (optarg)
				oOptions._seed = optarg;
			break;
		case 't':
			if (optarg)
				oOptions._threadscount = atoi(optarg);
			break;
		case 'v':
			oOptions._verbose = true;
			break;
		case 'w':
			if (optarg)
				oOptions._textfile = optarg;
			break;
		case 'x':
			if (optarg)
				oOptions._maxiter = atoi(optarg);
			break;
		}
	}
}

int main(int iArgc, char* iArgv[]) {
	try {
		std::cout << "cDecryptor Version 7.12.2020 20:33" << std::endl;
		signal(SIGINT, signalHandler);

		Options aOptions;
		parseOptions(iArgc, iArgv, aOptions);

		std::string aCipherString;
		std::cout << "Reading cipher file " << aOptions._cipherfile << std::endl;
		readCipher(aOptions._cipherfile, aCipherString);

		std::cout << "Cipher: " << aCipherString << std::endl;
		std::cout << "Cipher length: " << aCipherString.length() << std::endl;
		std::cout << "Randomize fraction: " << aOptions._random << std::endl;
		std::cout << "Random re-initialization after " << aOptions._maxiter << " iterations" << std::endl;
		std::cout << "Tolerance factor: " << aOptions._fuzzy << std::endl;
		std::cout << "Parallel threads: " << aOptions._threadscount << std::endl;

		std::unordered_map<unsigned long long, NGram*> aNorms;
		readNorms(aOptions._ngramsfiles, aNorms);

		std::cout << "Analyzing sample text file: " << aOptions._textfile << std::endl;
		computeScoreStatistics(aOptions._textfile, aNorms, aCipherString);

		printBestPossibleScore(aNorms);

		if (aOptions._seed.length()>0)
			std::cout << "Seed: " << RatedScore(Score(aNorms, aOptions._seed), aGlobalScoreStatistics) << " " << aOptions._seed << std::endl;

		aGlobalBestScore=RatedScore(Score(aNorms, std::string(aCipherString.length(), '.')), aGlobalScoreStatistics);

		std::vector<std::thread> aThreads[aOptions._threadscount];

		for (unsigned long long aThread=0; aThread<aOptions._threadscount; aThread++)
			aThreads->push_back(std::thread(&hillclimber, aThread, aNorms, aCipherString, aOptions));

		logTime(aOptions._threadscount, "threads started.");

		for (std::vector<std::thread>::iterator i=aThreads->begin(); i!=aThreads->end(); ++i)
			(*i).join();

		for (std::unordered_map<unsigned long long, NGram*>::iterator aI=aNorms.begin(); aI!=aNorms.end(); ++aI)
			delete aI->second;
	}
	catch (std::string& iString) {
		std::cerr << "Error caught: " << iString << std::endl;
	}
	return EXIT_SUCCESS;
}
