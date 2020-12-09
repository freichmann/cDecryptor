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
	std::string _seed="";
	unsigned int _maxiter=250;
	unsigned int _threadscount=1;
	unsigned int _diskSize=0;
	long double _random=0.005;
	long double _fuzzy=0.015;
	bool _verbose=false;
};

#endif /* SRC_OPTIONS_H_ */
