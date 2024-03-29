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

#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif

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

bool partiallyShuffleVec(std::vector<char>& iVector, const long double &iRandom) {
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

std::unordered_set<char> distinctSymbols(const std::string &iCipherString) {
	std::unordered_set<char> aSymbols;
	for (unsigned long long i = 0; i < iCipherString.length(); i++)
		aSymbols.insert(iCipherString.at(i));
	return aSymbols;
}

std::string buildClear(const std::string& iCipherString, const std::unordered_map<char, unsigned int>& iMap, const std::vector<char>& iVec, const Options& iOptions) {
	std::string iClear;
	for (unsigned int i=0; i<iCipherString.length(); i++) {
		if (iOptions._diskSize==0)
			iClear+=iVec.at(iMap.find(iCipherString[i])->second);
		else
			iClear+=iVec.at((iMap.find(iCipherString[i])->second + i ) % (signed int)iOptions._diskSize);
	}
	return iClear;
}

void randomVecInit(std::vector<char> &oVec) {
	oVec.clear();
	std::vector<char> aVec;

	for (char c='a'; c<='z'; c++)
		aVec.push_back(c);

	unsigned int aOriginalVectorSize=aVec.size();

	for (unsigned int aI=0; aI<aOriginalVectorSize; aI++) {
		std::uniform_int_distribution<unsigned int> aIntDistribution(0, aVec.size()-1);
		unsigned int aPos = aIntDistribution(aGlobalRandomEngine);
		oVec.push_back(*(aVec.begin()+aPos));
		aVec.erase(aVec.begin()+aPos);
	}
}

void randomMapInit(std::unordered_map<char, unsigned int> &oMap, const std::string &iCipherString, const unsigned int iSize) {
	oMap.clear();
	std::unordered_set<char> aSymbols = distinctSymbols(iCipherString);
	std::uniform_int_distribution<unsigned int> aSymbolsDistribution(0, iSize-1);

	for (std::unordered_set<char>::iterator i = aSymbols.begin(); i != aSymbols.end(); ++i)
		oMap.insert(std::pair<char, unsigned int>(*i, aSymbolsDistribution(aGlobalRandomEngine)));
}

void randomMapVecInit(std::unordered_map<char, unsigned int>& oMap, std::vector<char>& oVec, const std::string& iCipherString) {
	randomVecInit(oVec);
	randomMapInit(oMap, iCipherString, oVec.size());
}

std::string concat(const std::vector<char>& iVec) {
	std::string aString;
	for (std::vector<char>::const_iterator i = iVec.begin(); i != iVec.end(); ++i)
		aString += (*i);
	return aString;
}

void insertSymbols(std::unordered_map<char, unsigned int>& oMap, std::vector<char>& oVec, const std::string& iCipherString, const Options& iOptions) {
	oVec.clear();
	if (iOptions._diskSize>0)
		for (unsigned int i=0; i<iOptions._seedmap.length(); i++)
			oVec.push_back(iOptions._seedmap.at(i));
	else
		randomVecInit(oVec);

	randomMapInit(oMap, iCipherString, oVec.size());

	for (unsigned int i=0; i<std::min(iOptions._seed.length(), iCipherString.length()); i++) {
		std::vector<char>::iterator aVecIndex=std::find(oVec.begin(), oVec.end(), iOptions._seed.at(i));

		if (aVecIndex != oVec.end()) {
			if (iOptions._diskSize==0)
				oMap.find(iCipherString[i])->second=aVecIndex-oVec.begin();
			else
				oMap.find(iCipherString[i])->second=(((aVecIndex-oVec.begin()) % iOptions._diskSize) - (i % iOptions._diskSize) + iOptions._diskSize) % iOptions._diskSize;
		}
	}
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

void printGlobalBest(const std::string &iCipher, const std::unordered_map<char, unsigned int> &iMap, const std::vector<char> &iVector, const Options &iOptions, const unsigned long long &iThread, const RatedScore &iScore) {
	std::string iClear = buildClear(iCipher, iMap, iVector, iOptions);
	if (iOptions._diskSize == 0)
		logTime("Thread:", iThread, "Score:", iScore, iClear);
	else
		logTime("Thread:", iThread, "Score:", iScore, iClear, concat(iVector));
}

// not thread safe
void setGlobalBest(const RatedScore &iScore, const std::unordered_map<char, unsigned int> &iMap, const std::vector<char> &iVector) {
	aGlobalBestScore = iScore;
	aGlobalBestMap = iMap;
	aGlobalBestVector = iVector;
}

bool acceptIfGlobalBest(const RatedScore& iScore, const std::vector<char>& iVector, const std::unordered_map<char, unsigned int>& iMap) {
	Lock aLock(aGlobalBestScoreMutex);
	if (iScore>aGlobalBestScore) {
		setGlobalBest(iScore, iMap, iVector);
		return true;
	}
	return false;
}

bool printIfGlobalBest(const RatedScore& iScore, const std::string& iCipher, const unsigned long long& iThread, const std::vector<char>& iVector, const std::unordered_map<char, unsigned int>& iMap, const Options& iOptions) {
	if (acceptIfGlobalBest(iScore, iVector, iMap)) {
		printGlobalBest(iCipher, iMap, iVector, iOptions, iThread, iScore);
		return true;
	}
	return false;
}

Decision tolerateCandidate(const RatedScore& iBestScore, const RatedScore& iCandidateScore, long double iTemperature) {
	if (iCandidateScore>iBestScore)
		return ACCEPT;
	else {
		if (
				iBestScore != iCandidateScore
				&&
				aGlobalConstantDistribution(aGlobalRandomEngine) < std::exp(-(iBestScore.value()-iCandidateScore.value())/(-iBestScore.value()*iTemperature))
		)
			return TOLERATE;
		else
			return REJECT;
	}
}

RatedScore optimizeSymbols(
		std::unordered_map<char, unsigned int>& oCandidateMap,
		const std::string& iCipherString,
		const Options& iOptions,
		const std::unordered_map<unsigned long long, NGram*>& iNorms,
		const std::vector<char>& iCandidateVector,
		const unsigned long long& iThread) {
	RatedScore aBestScore;
	bool aFirst=true;
	bool aSymbolsChanged;

	do {
		aSymbolsChanged = false;
		for (std::unordered_map<char, unsigned int>::iterator aMappedSymbol = oCandidateMap.begin(); aMappedSymbol != oCandidateMap.end(); ++aMappedSymbol) {
			const unsigned int aBefore = aMappedSymbol->second;
			unsigned int aBestSymbolSoFar = aBefore;
			for (unsigned int i = 0; i < iCandidateVector.size(); i++)
				if (i != aBefore) {
					aMappedSymbol->second = i;
					std::string aCandidateString = buildClear(iCipherString, oCandidateMap, iCandidateVector, iOptions);
					RatedScore aCandidateScore(Score(iNorms, aCandidateString), aGlobalScoreStatistics);
					if (aCandidateScore>aBestScore || aFirst) {
						aBestScore=aCandidateScore;
						aBestSymbolSoFar = i;
						aSymbolsChanged = true;
						aFirst=false;
					}
				}
			aMappedSymbol->second = aBestSymbolSoFar;
		}
	} while (aSymbolsChanged);

	return aBestScore;
}

void hillclimber(const unsigned long long& iThread,
		const std::unordered_map<unsigned long long, NGram*>& iNorms,
		const std::string& iCipher,
		const Options& iOptions) {
	unsigned long long aCounterUntilReset=iOptions._maxiter;
	long double aTemperature;
	unsigned long long aGiveUp=iOptions._giveUp;

	while (iOptions._giveUp==0 || aGiveUp<=iOptions._giveUp ) {
		std::vector<char> aClimberVec=getGlobalBest().first;
		std::unordered_map<char, unsigned int> aClimberMap=getGlobalBest().second;
		RatedScore aClimberScore(Score(iNorms, buildClear(iCipher, aClimberMap, aClimberVec, iOptions)), aGlobalScoreStatistics);

		aCounterUntilReset=iOptions._maxiter;
		aTemperature=1.0;

		while (aCounterUntilReset>0) {
			{
				std::unordered_map<char, unsigned int> aCandidateMap=aClimberMap;
				std::vector<char> aCandidateVec=aClimberVec;

				while (!partiallyShuffleVec(aCandidateVec, aTemperature*iOptions._random));

				if (iOptions._verbose) {
					std::string aString=buildClear(iCipher, aCandidateMap, aCandidateVec, iOptions);
					RatedScore aScore=RatedScore(Score(iNorms, aString), aGlobalScoreStatistics);
					logTime("DEBUG Thread:", iThread, "Start Map", "Score:", aScore, "Counter until reset:", aCounterUntilReset, "Cool down:", aTemperature, aString, concat(aCandidateVec));
				}

				RatedScore aCandidateScore=optimizeSymbols(aCandidateMap, iCipher, iOptions, iNorms, aCandidateVec, iThread);

				if (aCandidateScore>aClimberScore) {
					aClimberScore=aCandidateScore;
					aClimberMap=aCandidateMap;
					aClimberVec=aCandidateVec;
					aCounterUntilReset=iOptions._maxiter;
					printIfGlobalBest(aClimberScore, iCipher, iThread, aClimberVec, aClimberMap, iOptions);
				}

				if (iOptions._verbose) {
					std::string aString=buildClear(iCipher, aCandidateMap, aCandidateVec, iOptions);
					RatedScore aScore=RatedScore(Score(iNorms, aString), aGlobalScoreStatistics);
					logTime("DEBUG Thread:", iThread, "Stuck Map", "Score:", aScore, "Counter until reset:", aCounterUntilReset, "Cool down:", aTemperature, aString, concat(aCandidateVec));
				}
			}

			if (iOptions._diskSize>0) {
				bool aVecImprovement;
				do {
					if (iOptions._verbose) {
						std::string aString=buildClear(iCipher, aClimberMap, aClimberVec, iOptions);
						RatedScore aScore=RatedScore(Score(iNorms, aString), aGlobalScoreStatistics);
						logTime("DEBUG Thread:", iThread, "Start Vec", "Score:", aScore, "Counter until reset:", aCounterUntilReset, "Cool down:", aTemperature, aString, concat(aClimberVec));
					}

					aVecImprovement=false;
					std::vector<char> aCandidateVec=aClimberVec;
					for (std::vector<char>::iterator aFrom=aCandidateVec.begin(); aFrom!=aCandidateVec.end(); ++aFrom)
						for (std::vector<char>::iterator aTo=aFrom+1; aTo!=aCandidateVec.end(); ++aTo) {
							iter_swap(aFrom, aTo);
							if (iOptions._verbose)
								std::cout << "Testing vector: " << concat(aCandidateVec) << std::endl;

							std::unordered_map<char, unsigned int> aCandidateMap(aClimberMap);
							RatedScore aCandidateScore=optimizeSymbols(aCandidateMap, iCipher, iOptions, iNorms, aCandidateVec, iThread);

							if (aCandidateScore>aClimberScore) {
								aVecImprovement=true;
								aClimberScore=aCandidateScore;
								aClimberMap=aCandidateMap;
								aClimberVec=aCandidateVec;
								aCounterUntilReset=iOptions._maxiter;
								printIfGlobalBest(aClimberScore, iCipher, iThread, aClimberVec, aClimberMap, iOptions);
							} else
								iter_swap(aFrom, aTo);
						}

					if (iOptions._verbose) {
						std::string aString=buildClear(iCipher, aClimberMap, aClimberVec, iOptions);
						RatedScore aScore=RatedScore(Score(iNorms, aString), aGlobalScoreStatistics);
						logTime("DEBUG Thread:", iThread, "Stuck Vec", "Score:", aScore, "Counter until reset:", aCounterUntilReset, "Cool down:", aTemperature, aString, concat(aClimberVec));
					}

				} while (aVecImprovement);
			}

			aCounterUntilReset--;
			{
				long double aDouble=(long double)aCounterUntilReset/(long double)iOptions._maxiter;
				if (aDouble<aTemperature)
					aTemperature=aDouble;
			}
		}

		aGiveUp++;
		if (iOptions._verbose)
			logTime("DEBUG Thread:", iThread, "Iteration", aGiveUp, "Re-shuffling from best known");
	}
}

void printBestPossibleScore(std::unordered_map<unsigned long long, NGram*>& iNorms) {
	long double aLnPerfect = 0.0;
	for (std::unordered_map<unsigned long long, NGram*>::iterator i=iNorms.begin(); i != iNorms.end(); ++i) {
		long double aLnNGramPerfect = -logl(sqrtl(2.0 * M_PI) * aGlobalScoreStatistics.at(i->first)._sigma);
		std::cout
		<< setiosflags(std::ios::fixed)
		<< "NGram length:" << i->second->_length << " NGrams:"
		<< i->second->_NGramMap.size() << " Samples:"
		<< i->second->_count << " Mean:" << aGlobalScoreStatistics.at(i->first)._mean
		<< " StdDev:" << aGlobalScoreStatistics.at(i->first)._sigma << " Optimum: "
		<< aLnNGramPerfect << std::endl;
		aLnPerfect += aLnNGramPerfect;
	}
	std::cout << "Score optimum: " << aLnPerfect << std::endl;
}

void signalHandler(const int iSigNum) {
	std::cout << "Received signal " << iSigNum << std::endl;
	exit(iSigNum);
}

void parseOptions(const int iArgc, char* iArgv[], Options& oOptions) {
	int aInt;
	while ((aInt = getopt(iArgc, iArgv, "c:d:f:g:l:m:p:r:s:t:vw:x:z:")) != -1) {
		switch (aInt) {
		case 'c':
			if (optarg)
				oOptions._cipherfile = optarg;
			break;
		case 'd':
			if (optarg)
				oOptions._diskSize = atoi(optarg);
			break;
		case 'g':
			if (optarg)
				oOptions._giveUp = atoi(optarg);
			break;
		case 'l':
			if (optarg)
				oOptions._ngramsfiles.push_back(optarg);
			break;
		case 'm':
			if (optarg)
				oOptions._seedmap = optarg;
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
		std::vector<unsigned long long> aVector;
		std::stringstream aStringStream(aTransposition);
		for (int i; aStringStream >> i;) {
			aVector.push_back(i);
			if (aStringStream.peek() == ',')
				aStringStream.ignore();
		}
		for (std::size_t i = 0; i < std::min(aVector.size(),iCipherString.length()); i++)
			aTransposed += iCipherString.at(aVector.at(i));
	} else
		throw "Failed to open " + iFileName;
	return aTransposed;
}

void printCipherStats(std::string& aCipherString) {
	std::cout << "Cipher: " << aCipherString << std::endl;
	std::cout << "Cipher length: " << aCipherString.length() << std::endl;
	std::unordered_set<char> aSymbols;
	for (unsigned long long i = 0; i < aCipherString.length(); i++)
		aSymbols.insert(aCipherString.at(i));
	unsigned int N = aSymbols.size();
	std::cout << "Cipher Symbols count: " << N << std::endl;
	if (N>26)
		std::cout << "Unicity distance for homophonic ciphers: "
		<< 1.47*N-49.91+0.45*(0.5 + N)*std::log(N)-0.45*(N-25.5)* std::log(N-26) << std::endl;
}

int main(int iArgc, char* iArgv[]) {
	try {
		std::cout << "cDecryptor Version 12.10.2021 14:23" << std::endl;
		std::cout << std::setprecision(16);
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
			std::cout << "Cipher string after transposition: " << aCipherString << std::endl;
			std::cout << "Cipher length after transposition: " << aCipherString.length() << std::endl;
		}

		std::cout << "Randomize ratio in iteration: " << aOptions._random << std::endl;
		if (aOptions._random<=0 || aOptions._random>1.0) {
			std::cerr << "Random ratio must be smaller or equal than 1.0 and larger than 0.0." << std::endl;
			return EXIT_FAILURE;
		}

		std::cout << "Re-shuffling from best known solution after " << aOptions._maxiter << " iterations without improvement" << std::endl;

		if (aOptions._giveUp==0)
			std::cout << "Never give up with";
		else
			std::cout << "Give up after " << aOptions._giveUp;
		std::cout << " re-shufflings" << std::endl;

		if (aOptions._threadscount==0)
			aOptions._threadscount = std::thread::hardware_concurrency();
		std::cout << "Parallel threads: " << aOptions._threadscount << std::endl;

		std::unordered_map<unsigned long long, NGram*> aNorms;
		readNorms(aOptions._ngramsfiles, aNorms);

		if (aOptions._seed.length()>0) {
			std::cout << "Seed: " << aOptions._seed << std::endl;
			if (aOptions._diskSize>0) {
				if (aOptions._seedmap.length()>0)
					std::cout << "Map: " << aOptions._seedmap << std::endl;
				else {
					std::cout << "Seed for chiffre disk only with providing map." << std::endl;
					return EXIT_FAILURE;
				}
			}
		}

		std::cout << "Building statistics of scores of sample text file " << aOptions._textfile << " for snippets of length " << aCipherString.length() << std::endl;
		computeScoreStatistics(aOptions._textfile, aNorms, aCipherString);

		printBestPossibleScore(aNorms);

		{
			std::vector<char> aVector;
			std::unordered_map<char, unsigned int> aMap;

			if (aOptions._seed.length()>0)
				insertSymbols(aMap, aVector, aCipherString, aOptions);
			else
				randomMapVecInit(aMap, aVector, aCipherString);

			optimizeSymbols(aMap, aCipherString, aOptions, aNorms, aVector, 0);

			setGlobalBest(RatedScore(Score(aNorms, buildClear(aCipherString, aMap, aVector, aOptions)), aGlobalScoreStatistics), aMap, aVector);
		}

		std::vector<std::thread> aThreads[aOptions._threadscount];

		for (unsigned long long aThread=0; aThread<aOptions._threadscount; aThread++)
			aThreads->push_back(std::thread(&hillclimber, aThread, aNorms, aCipherString, aOptions));

		logTime(aOptions._threadscount, "threads started.");

		for (std::vector<std::thread>::iterator i=aThreads->begin(); i!=aThreads->end(); ++i)
			i->join();

		for (std::unordered_map<unsigned long long, NGram*>::iterator aI=aNorms.begin(); aI!=aNorms.end(); ++aI)
			delete aI->second;
	}
	catch (std::string& iString) {
		std::cerr << "Error: " << iString << std::endl;
	}
	catch (...) {
		std::cerr << "Unspecified error." << std::endl;
	}

	std::cout << "Program terminated." << std::endl;

	return EXIT_SUCCESS;
}
