/*
 * Lock.cpp
 *
 *  Created on: 11.10.2020
 *      Author: Fritz Reichmann
 */

#include "Lock.h"

Lock::Lock(std::mutex &m) : m_(m) {
	m.lock();
}

Lock::~Lock() {
	m_.unlock();
}
