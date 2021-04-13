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
#include <sstream>

#include "NGram.h"
#include "Lock.h"
#include "RatedScore.h"
#include "Options.h"

enum Decision { REJECT, TOLERATE, ACCEPT };

//Globals
std::mutex aGlobalBestScoreMutex;
std::mutex aGlobalOutputMutex;
RatedScore aGlobalBestScore;
std::unordered_map<char, unsigned int> aGlobalBestMap;
std::vector<char> aGlobalBestVector;
std::unordered_map<unsigned long long, GaussianNorm> aGlobalScoreStatistics;
std::default_random_engine aGlobalRandomEngine;
std::uniform_real_distribution<long double> aGlobalConstantDistribution(0,1);

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

bool partiallyShuffleMap(std::unordered_map<char, unsigned int>& iSymbolMap, const long double &iRandom) {
	static std::uniform_int_distribution<unsigned int> aCharDistribution(0, (unsigned int)'z'-(unsigned int)'a');
	bool aShuffled=false;

	for (std::unordered_map<char, unsigned int>::iterator i=iSymbolMap.begin(); i!=iSymbolMap.end(); ++i)
		if (aGlobalConstantDistribution(aGlobalRandomEngine)<iRandom) {
			i->second=(char)(aCharDistribution(aGlobalRandomEngine));
			aShuffled=true;
		}
	return aShuffled;
}

bool partiallyShuffleLetters(std::vector<char>& iVector, const long double &iRandom) {
	std::uniform_int_distribution<unsigned int> aCharDistribution(0, iVector.size());
	bool aShuffled=false;

	for (std::vector<char>::iterator aFrom=iVector.begin(); aFrom!=iVector.end(); ++aFrom)
		for (std::vector<char>::iterator aTo=iVector.begin(); aTo!=iVector.end(); ++aTo)
			if (aGlobalConstantDistribution(aGlobalRandomEngine)<iRandom) {
				std::iter_swap(aFrom, aTo);
				aShuffled=true;
			}
	return aShuffled;
}

void randomMapInit(const std::string& iCipher, std::unordered_map<char, unsigned int>& oSymbolMap, std::vector<char>& oLetterVec) {
	oSymbolMap.clear();
	oLetterVec.clear();

	std::unordered_set<char> aSymbols;
	for (unsigned long long i=0; i<iCipher.length(); i++)
		aSymbols.insert(iCipher.at(i));

	std::uniform_int_distribution<unsigned int> aSymbolsDistribution(0, (unsigned int)'z'-(unsigned int)'a');
	for (std::unordered_set<char>::iterator i=aSymbols.begin(); i!=aSymbols.end(); ++i)
		oSymbolMap.insert(std::pair<char, unsigned int>(*i, aSymbolsDistribution(aGlobalRandomEngine)));

	std::vector<char> aVec;
	for (char c='a'; c<='z'; c++)
		aVec.push_back(c);
	unsigned int aOriginalVectorSize=aVec.size();

	for (unsigned int aI=0; aI<aOriginalVectorSize; aI++) {
		std::uniform_int_distribution<unsigned int> aIntDistribution(0, aVec.size()-1);
		unsigned int aPos=aIntDistribution(aGlobalRandomEngine);
		oLetterVec.push_back(*(aVec.begin()+aPos));
		aVec.erase(aVec.begin()+aPos);
	}
}

std::string buildClear(const std::string& iCipherString, const std::unordered_map<char, unsigned int>& iSymbolMap, const std::vector<char>& iLetterVec, const Options& iOptions) {
	std::string iClear;
	for (unsigned int i=0; i<iCipherString.length(); i++) {
		if (iOptions._diskSize==0)
			iClear+=iLetterVec.at(iSymbolMap.find(iCipherString[i])->second);
		else
			iClear+=iLetterVec.at((iSymbolMap.find(iCipherString[i])->second + i ) % iOptions._diskSize);
	}
	return iClear;
}

void insertSymbols(std::unordered_map<char,unsigned int>& iSymbolMap, const std::string& iCipherString, std::vector<char>& iLetterVec, const Options& iOptions) {
	std::vector<char> aChars;
	for (unsigned int i=0; i<std::min(iOptions._seed.length(), iCipherString.length()); i++) {
		std::vector<char>::iterator aIndex = std::find(aChars.begin(), aChars.end(), iOptions._seed.at(i));
		if (aIndex == aChars.end()) {
			aChars.push_back(iOptions._seed.at(i));
			aIndex=aChars.end()-1;
		}
		if (iOptions._diskSize==0)
			iSymbolMap.find(iCipherString[i])->second=std::distance(aChars.begin(), aIndex);
		else
			iSymbolMap.find(iCipherString[i])->second=(std::distance(aChars.begin(), aIndex)-i) % iOptions._diskSize;
	}
	for (std::vector<char>::iterator aIt=iLetterVec.begin(); aIt!=iLetterVec.end(); ++aIt)
		if (std::find(aChars.begin(), aChars.end(), *aIt)==aChars.end())
			aChars.push_back(*aIt);
	iLetterVec=aChars;
}

bool checkIfGlobalBest(const RatedScore& iRatedScore, const std::vector<char>& iVector, const std::unordered_map<char, unsigned int>& iMap) {
	Lock aLock(aGlobalBestScoreMutex);
	if (iRatedScore>aGlobalBestScore) {
		aGlobalBestScore = iRatedScore;
		aGlobalBestMap = iMap;
		aGlobalBestVector = iVector;
		return true;
	}
	return false;
}

std::pair<const std::vector<char>, const std::unordered_map<char, unsigned int>> getGlobalBest() {
	Lock aLock(aGlobalBestScoreMutex);

	return std::pair<const std::vector<char>, const std::unordered_map<char, unsigned int>>(aGlobalBestVector, aGlobalBestMap);
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

std::string concat(const std::vector<char>& iCandidateLetterVector) {
	std::string aChiffreDisk;
	for (std::vector<char>::const_iterator i = iCandidateLetterVector.begin();
			i != iCandidateLetterVector.end(); ++i)
		aChiffreDisk += (*i);
	return aChiffreDisk;
}

bool printIfGlobalBest(const RatedScore& iScore, const std::string& iCipher, const unsigned long long& iThread, const std::vector<char>& iVector, const std::unordered_map<char, unsigned int>& iMap, const Options& iOptions) {
	if (checkIfGlobalBest(iScore, iVector, iMap)) {
		std::string iClear=buildClear(iCipher, iMap, iVector, iOptions);
		if (iOptions._diskSize==0)
			logTime("Thread:", iThread, "Score:", iScore, "Clear:", iClear);
		else
			logTime("Thread:", iThread, "Score:", iScore, "Clear:", iClear, "Chiffredisk:", concat(iVector));
		return true;
	}
	return false;
}

bool acceptCandidate(const RatedScore& iPreviousScore, const RatedScore& iCandidateScore) {
	return (std::isnan(iPreviousScore.value()) || iCandidateScore > iPreviousScore);
}

Decision tolerateCandidate(const RatedScore& iBestScore, const RatedScore& iCandidateScore, long double iTemperature) {
	//	std::cout << "BestScore: " << iBestScore.value() << std::endl;
	//	std::cout << "CandidateScore: " << iCandidateScore.value() << std::endl;
	//	std::cout << "Temperature: " << iTemperature << std::endl;
	//	std::cout << "tolerance threshold: " << std::exp(-(iBestScore.value()-iCandidateScore.value())/(-iBestScore.value()*iTemperature)) << std::endl;

	if (acceptCandidate(iBestScore, iCandidateScore))
		return ACCEPT;
	else {
		if (
				iBestScore.value() != iCandidateScore.value()
				&&
				aGlobalConstantDistribution(aGlobalRandomEngine) < std::exp(-(iBestScore.value()-iCandidateScore.value())/(-iBestScore.value()*iTemperature))
		)
			return TOLERATE;
		else
			return REJECT;
	}
}

RatedScore optimizeSymbols(const std::string& iCipherString,
		const Options& iOptions,
		const std::unordered_map<unsigned long long, NGram*>& iNorms,
		std::unordered_map<char, unsigned int>& oCandidateMap,
		std::vector<char>& iCandidateVector,
		const unsigned long long& iThread) {
	RatedScore aBestScore;
	bool aSymbolsChanged;
	do {
		aSymbolsChanged = false;
		for (std::unordered_map<char, unsigned int>::iterator aMappedSymbol =
				oCandidateMap.begin(); aMappedSymbol != oCandidateMap.end();
				++aMappedSymbol) {
			const unsigned int aBefore = aMappedSymbol->second;
			unsigned int aBestSymbolSoFar = aBefore;
			for (unsigned int i = 0; i < iCandidateVector.size(); i++)
				if (i != aBefore) {
					aMappedSymbol->second = i;
					std::string aCandidateString = buildClear(iCipherString, oCandidateMap, iCandidateVector, iOptions);
					RatedScore aCandidateScore( Score(iNorms, aCandidateString), aGlobalScoreStatistics);
					if (acceptCandidate(aBestScore, aCandidateScore)) {
						aBestScore=aCandidateScore;
						aBestSymbolSoFar = i;
						aSymbolsChanged = true;
					}
				}
			aMappedSymbol->second = aBestSymbolSoFar;
		}
	} while (aSymbolsChanged);

	return aBestScore;
}

void hillclimber(const unsigned long long& iThread,
		const std::unordered_map<unsigned long long, NGram*>& iNorms,
		const std::string& iCipherString,
		const Options& iOptions) {
	unsigned long long aCounterUntilReset=iOptions._maxiter;
	long double aTemperature=1.0;
	RatedScore aClimberBestScore;
	std::unordered_map<char, unsigned int> aCandidateMap;
	std::vector<char> aCandidateVector;

	randomMapInit(iCipherString, aCandidateMap, aCandidateVector);

	if (iThread==0) {
		insertSymbols(aCandidateMap, iCipherString, aCandidateVector, iOptions);
		if (iOptions._verbose)
			std::cout << "Thread " << iThread << " applied with seed " << buildClear(iCipherString, aCandidateMap, aCandidateVector, iOptions) << std::endl;
	}

	while (true) {
		RatedScore aLoopBestScore(Score(iNorms, buildClear(iCipherString, aCandidateMap, aCandidateVector, iOptions)), aGlobalScoreStatistics);
		std::unordered_map<char, unsigned int> aLoopBestMap;
		std::vector<char> aLoopBestVector;

		while (aCounterUntilReset>0 && iOptions._random>0) {
			std::string aCandidateString=buildClear(iCipherString, aCandidateMap, aCandidateVector, iOptions);
			RatedScore aCandidateScore=RatedScore(Score(iNorms, aCandidateString), aGlobalScoreStatistics);

			if (iOptions._verbose)
				logTime("DEBUG Thread:", iThread, "Start", "Score:", aCandidateScore, "Counter until reset:", aCounterUntilReset, "Cool down:", aTemperature, aCandidateString, concat(aCandidateVector));

			aCandidateScore=optimizeSymbols(iCipherString, iOptions, iNorms, aCandidateMap, aCandidateVector, iThread);

			if (iOptions._verbose && iOptions._diskSize>0)
				std::cout << "Symbol step: " << buildClear(iCipherString, aCandidateMap, aCandidateVector, iOptions) << std::endl;

			if (acceptCandidate(aLoopBestScore, aCandidateScore)) {
				aLoopBestScore=aCandidateScore;
				aLoopBestMap=aCandidateMap;
				aLoopBestVector=aCandidateVector;
				if (acceptCandidate(aClimberBestScore, aLoopBestScore)) {
					aClimberBestScore=aLoopBestScore;
					aCounterUntilReset=iOptions._maxiter;
				}
				printIfGlobalBest(aCandidateScore, iCipherString, iThread, aCandidateVector, aCandidateMap, iOptions);
			} else {
				if (TOLERATE==tolerateCandidate( aClimberBestScore, aCandidateScore, aTemperature)) {
					aLoopBestScore=aCandidateScore;
					aLoopBestMap=aCandidateMap;
					aLoopBestVector=aCandidateVector;
				}
			}

			if (iOptions._diskSize>0) {
				bool aVectorImprovement;
				do {
					aVectorImprovement=false;
					if (iOptions._verbose && aVectorImprovement)
						std::cout << "Begin vector iteration: " << buildClear(iCipherString, aCandidateMap, aCandidateVector, iOptions) << std::endl;
					for (std::vector<char>::iterator aFrom=aCandidateVector.begin(); aFrom!=aCandidateVector.end(); ++aFrom)
						for (std::vector<char>::iterator aTo=aFrom+1; aTo!=aCandidateVector.end(); ++aTo) {
							iter_swap(aFrom, aTo);
							if (iOptions._verbose)
								std::cout << "Testing vector: " << concat(aCandidateVector) << std::endl;
							std::unordered_map<char, unsigned int> aInnerCandidateMap(aCandidateMap);
							std::vector<char> aInnerCandidateVector(aCandidateVector);

							RatedScore aMapCandidateScore=optimizeSymbols(iCipherString, iOptions, iNorms, aInnerCandidateMap, aInnerCandidateVector, iThread);

							if (acceptCandidate(aLoopBestScore, aMapCandidateScore)) {
								aLoopBestScore=aMapCandidateScore;
								aLoopBestMap=aCandidateMap;
								aLoopBestVector=aInnerCandidateVector;
								if (acceptCandidate(aClimberBestScore, aLoopBestScore)) {
									aVectorImprovement=true;
									aClimberBestScore=aLoopBestScore;
									aCounterUntilReset=iOptions._maxiter;
								}
								printIfGlobalBest(aMapCandidateScore, iCipherString, iThread, aCandidateVector, aCandidateMap, iOptions);
							} else {
								if (REJECT!=tolerateCandidate( aClimberBestScore, aMapCandidateScore, aTemperature)) {
									aLoopBestScore=aMapCandidateScore;
									aLoopBestMap=aCandidateMap;
									aLoopBestVector=aCandidateVector;
								} else
									iter_swap(aFrom, aTo);
							}
						}
					if (iOptions._verbose && aVectorImprovement)
						std::cout << "End vector iteration: " << buildClear(iCipherString, aCandidateMap, aCandidateVector, iOptions) << std::endl;

				} while (aVectorImprovement);

				if (iOptions._verbose && iOptions._diskSize>0)
					std::cout << "Vector step: " << buildClear(iCipherString, aCandidateMap, aCandidateVector, iOptions) << std::endl;
			}

			if (iOptions._verbose) {
				std::string aClear=buildClear(iCipherString, aCandidateMap, aCandidateVector, iOptions);
				aCandidateScore=RatedScore(Score(iNorms, aClear), aGlobalScoreStatistics);

				logTime("DEBUG Thread:", iThread, "Stuck", "Score:", aCandidateScore, "Counter until reset:", aCounterUntilReset, "Cool down:", aTemperature, aClear, concat(aCandidateVector));
			}

			if (aTemperature*iOptions._random>0) {
				std::pair<std::vector<char>, std::unordered_map<char, unsigned int>> aBest=getGlobalBest();
				aCandidateMap=aBest.second;
				aCandidateVector=aBest.first;
				bool aShuffled=false;
				while (!aShuffled) {
					aShuffled |= partiallyShuffleMap(aCandidateMap, aTemperature*iOptions._random/2.0);
					aShuffled |= partiallyShuffleLetters(aCandidateVector, aTemperature*iOptions._random/4.0);
				}
			}

			aCounterUntilReset--;
			{
				long double aDouble=(long double)aCounterUntilReset/(long double)iOptions._maxiter;
				if (aDouble<aTemperature)
					aTemperature=aDouble;
			}
		}
		randomMapInit(iCipherString, aCandidateMap, aCandidateVector);
		aCounterUntilReset=iOptions._maxiter;
		aTemperature=1.0;
	}
}

void printBestPossibleScore(std::unordered_map<unsigned long long, NGram*>& iNorms) {
	long double aLnPerfect = 0.0;
	for (std::unordered_map<unsigned long long, NGram*>::iterator i=iNorms.begin(); i != iNorms.end(); ++i) {
		long double aLnNGramPerfect = -logl(sqrtl(2.0 * M_PI) * aGlobalScoreStatistics.at(i->first)._sigma);
		std::cout << setiosflags(std::ios::fixed)
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
	std::cout << "Received signal " << iSigNum << std::endl;
	exit(iSigNum);
}

void parseOptions(const int iArgc, char* iArgv[], Options& oOptions) {
	int aInt;
	while ((aInt = getopt(iArgc, iArgv, "c:d:f:l:p:r:s:t:vw:x:z:")) != -1) {
		switch (aInt) {
		case 'c':
			if (optarg)
				oOptions._cipherfile = optarg;
			break;
		case 'd':
			if (optarg)
				oOptions._diskSize = atoi(optarg);
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
		case 'z':
			if (optarg)
				oOptions._transpositionfile = optarg;
			break;
		}
	}
}

std::string transpose(const std::string& iCipherString, const std::string& iFileName) {
	std::cout << "Reading transposition file: " << iFileName
			<< std::endl;
	std::string aTransposition;
	std::ifstream aFile(iFileName);
	std::string aTransposed;

	if (aFile.is_open()) {
		while (!aFile.eof()) {
			std::string aLine;
			aFile >> aLine;
			aTransposition += aLine;
		}
		aFile.close();
		std::cout << "Applying transposition " << aTransposition << std::endl;
		std::vector<unsigned int> aVector;
		std::stringstream aStringStream(aTransposition);
		for (int i; aStringStream >> i;) {
			aVector.push_back(i);
			if (aStringStream.peek() == ',')
				aStringStream.ignore();
		}
		for (std::size_t i = 0; i < aVector.size(); i++)
			aTransposed += iCipherString.at(aVector.at(i));
	} else
		throw "Failed to open " + iFileName;
	return aTransposed;
}

void printCipherStats(std::string& aCipherString) {
	std::cout << "Cipher: " << aCipherString << std::endl;
	std::cout << "Cipher length: " << aCipherString.length() << std::endl;
	std::unordered_set<char> aSymbols;
	for (unsigned long int i = 0; i < aCipherString.length(); i++)
		aSymbols.insert(aCipherString.at(i));
	unsigned int N = aSymbols.size();
	std::cout << "Cipher Symbols count: " << N << std::endl;
	std::cout << "Unicity distance for homophonic ciphers: "
			<< 1.47*N-49.91+0.45*(0.5 + N)*std::log(N)-0.45*(N-25.5)* std::log(N-26) << std::endl;
}

int main(int iArgc, char* iArgv[]) {
	try {
		std::cout << "cDecryptor Version 12.4.2020 22:42" << std::endl;
		std::cout << std::setprecision(17);
		signal(SIGINT, signalHandler);
		aGlobalRandomEngine.seed(std::chrono::system_clock::now().time_since_epoch().count());

		Options aOptions;
		parseOptions(iArgc, iArgv, aOptions);

		std::string aCipherString;
		std::cout << "Reading cipher file " << aOptions._cipherfile << std::endl;
		readCipher(aOptions._cipherfile, aCipherString);
		printCipherStats(aCipherString);

		if (aOptions._diskSize>0)
			std::cout << "Chiffre Disk Size " << aOptions._diskSize << std::endl;

		if (aOptions._transpositionfile.length()>0) {
			aCipherString=transpose(aCipherString, aOptions._transpositionfile);
			std::cout << "After transposition: " << aCipherString <<                     std::endl;
		}

		std::cout << "Randomize fraction: " << aOptions._random << std::endl;
		std::cout << "Random re-initialization after " << aOptions._maxiter << " iterations" << std::endl;
		std::cout << "Parallel threads: " << aOptions._threadscount << std::endl;

		std::unordered_map<unsigned long long, NGram*> aNorms;
		readNorms(aOptions._ngramsfiles, aNorms);

		std::cout << "Analyzing sample text file: " << aOptions._textfile << std::endl;
		computeScoreStatistics(aOptions._textfile, aNorms, aCipherString);

		printBestPossibleScore(aNorms);

		if (aOptions._seed.length()>0) {
			std::cout << "Seed: " << RatedScore(Score(aNorms, aOptions._seed), aGlobalScoreStatistics) << " " << aOptions._seed << std::endl;
			if (aOptions._diskSize>0)
				std::cout << "Warning: Providing a seed for chiffre disk is not properly implemented." << std::endl;
		}

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
	catch (...) {
		std::cerr << "Unhandled error caught." << std::endl;
	}

	return EXIT_SUCCESS;
}
