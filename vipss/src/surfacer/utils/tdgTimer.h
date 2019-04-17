#ifndef TDGTIMER_H
#define TDGTIMER_H

#include <stdio.h>
#include <sys/types.h>
#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif
//#include <time.h>
#include <utils/Pm_MeshLite.H>

class tdgTimer
{
private:
	//clock_t m_StartTime;
	//clock_t m_StopTime;
	//Array<clock_t> m_SplitTime;
#ifdef WIN32
	LARGE_INTEGER m_freq;
	LARGE_INTEGER m_StartTime;
	LARGE_INTEGER m_StopTime;
	Array<LARGE_INTEGER> m_SplitTimes;
#else
	timeval m_StartTime;
	timeval m_StopTime;
	Array<timeval> m_SplitTimes;
#endif
	int m_isRunning;

public:
    tdgTimer();
	~tdgTimer();
	//struct timeval t0;

	void Reset();
	void Start();
	void Stop();
	void Split();

	double GetStartTime();
	double GetStopTime();
	double GetElapsedTime();

	double GetSplitTime(const int i);

	double GetTotalSplitTime();
	double GetAvgSplitTime();
	double GetMinSplitTime();
	double GetMaxSplitTime();
};
#endif
