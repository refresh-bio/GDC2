/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#ifndef DISK_MUTEX_H
#define DISK_MUTEX_H

#include <mutex>
using namespace std;

class DiskMutex
{
	mutable mutex mtx;

public:
	void lock() 
	{
		mtx.lock();
	}

	void unlock()
	{
		mtx.unlock();
	}
};

#endif