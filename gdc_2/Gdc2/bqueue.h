/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#ifndef _B_QUEUE_H
#define _B_QUEUE_H
// Generic multithreading queues

#include "defs.h"
#include <queue>
#include <list>
#include <thread>
#include <mutex>
#include <condition_variable>

using namespace std;

// Multithreading queue with registering mechanism:
//   * The queue can report whether it is in wainitng for new data state or there will be no new data
template<typename T> class CBufferedRegisteringQueue
{
	typedef queue<T, list<T>> queue_t;
	queue_t q;
	bool is_completed;
	int n_producers;
	uint32_t n_elements;
	uint32_t n_size; // queue size
	uint32_t n_threshold; 
	double threshold_factor;

	mutable mutex mtx;								// The mutex to synchronise on
	condition_variable cv_queue_empty;
	condition_variable cv_queue_full;
	condition_variable cv_queue_threshold;

public:
	CBufferedRegisteringQueue(int _n_size, double _threshold_factor)
	{
		Restart(1, _n_size, _threshold_factor); // only one producer
	};

	~CBufferedRegisteringQueue()
	{};

	void Restart(int _n_producers, int _n_size, double _threshold_factor)
	{
		unique_lock<mutex> lck(mtx);

		is_completed = false;
		n_producers = _n_producers;
		n_elements = 0;
		n_size = _n_size;
		n_threshold = _n_size * _threshold_factor;
		threshold_factor = _threshold_factor;
	}

	bool IsEmpty()
	{
		lock_guard<mutex> lck(mtx);
		return n_elements == 0;
	}

	bool IsCompleted()
	{
		lock_guard<mutex> lck(mtx);
		return n_elements == 0 && n_producers == 0;
	}

	void MarkCompleted()
	{
		lock_guard<mutex> lck(mtx);
		n_producers--;
		if (!n_producers)
			cv_queue_empty.notify_all();
	}

	void Push(T data)
	{
		unique_lock<mutex> lck(mtx);
		if (n_elements == n_size) 
		{
			cv_queue_threshold.wait(lck, [this]{return (this->n_elements < this->n_threshold); });
		}
		bool was_empty = n_elements == 0;
		q.push(data);
		++n_elements;

		if (was_empty)
			cv_queue_empty.notify_all();
	}
	bool Pop(T &data)
	{
		unique_lock<mutex> lck(mtx);
		if (n_elements < n_threshold)
			cv_queue_threshold.notify_one();

		cv_queue_empty.wait(lck, [this]{return !this->q.empty() || !this->n_producers; });
		if (n_elements == 0)
			return false;
		data = q.front();
		q.pop();
		--n_elements;

		return true;
	}

	uint32_t GetSize()
	{
		return n_elements;
	}
};


#endif
