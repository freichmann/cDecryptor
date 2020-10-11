/*
 * Lock.h
 *
 *  Created on: 11.10.2020
 *      Author: fritz
 */

#include <mutex>

#ifndef SRC_LOCK_H_
#define SRC_LOCK_H_

class Lock {
	std::mutex &m_;

public:
	Lock(std::mutex &);
	~Lock();
};

#endif /* SRC_LOCK_H_ */
