/*
 * Options.h
 *
 *  Created on: 14.11.2020
 *      Author: Fritz Reichmann
 */

#ifndef SRC_OPTIONS_H_
#define SRC_OPTIONS_H_

#include <string>
#include <list>

struct Options {
	std::list<std::string> _ngramsfiles;
	std::string _cipherfile;
	std::string _textfile;
	std::string _transpositionfile="";
	std::string _seed="";
	std::string _seedmap="";
	unsigned int _maxiter=10;
	unsigned int _threadscount=0;
	unsigned int _diskSize=0;
	unsigned int _giveUp=0;
	long double _random=0.025;
	bool _verbose=false;
};

#endif /* SRC_OPTIONS_H_ */
