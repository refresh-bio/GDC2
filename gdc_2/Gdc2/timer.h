/*
 This file is a part of GDC software distributed under GNU GPL 2 licence.
 The homepage of the GDC project is http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=gdc&subpage=about
 
 Authors: Sebastian Deorowicz, Agnieszka Danek, Marcin Niemiec
 
 Version: 2.0
 Date   : 2015-March-05
 */

#ifndef _TIMER_H
#define _TIMER_H

#ifdef WIN32
#include <windows.h>

// **********************************************************
typedef struct 
{
	LARGE_INTEGER start;
	LARGE_INTEGER stop;
} stop_watch_t;

// **********************************************************
class CStopWatch 
{
	stop_watch_t timer;
	LARGE_INTEGER frequency;
	double LIToSecs( LARGE_INTEGER & L);

public:
	CStopWatch();
	void StartTimer( );
	void StopTimer( );
	double GetElapsedTime();
};

// **********************************************************
typedef struct
{
	ULARGE_INTEGER start;
	ULARGE_INTEGER stop;
} thread_watch_t;

// **********************************************************
class CThreadWatch
{
	thread_watch_t timer_kernel, timer_user;
	LARGE_INTEGER frequency;
	double LIToSecs( LARGE_INTEGER & L);

public:
	CThreadWatch();
	void StartTimer( );
	void StopTimer( );
	double GetElapsedTime();
};

#else
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

typedef struct 
{
	timeval start;
	timeval stop;
} stop_watch_t;

class CStopWatch 
{
	stop_watch_t timer;

public:
	CStopWatch() {};
	void StartTimer( );
	void StopTimer( );
	double GetElapsedTime();
};

typedef timeval thread_watch_t;

// **********************************************************
class CThreadWatch
{
	thread_watch_t start_kernel, start_user;
	thread_watch_t stop_kernel , stop_user;

public:
	CThreadWatch();
	void StartTimer( );
	void StopTimer( );
	double GetElapsedTime();
};

#endif

#endif
