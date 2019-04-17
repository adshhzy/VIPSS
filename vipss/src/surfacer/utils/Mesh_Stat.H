// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Stat_h
#define Stat_h

#if 0
{
	NEST { STAT(Svdeg); for (i=0;i<10;i++) Svdeg+=vdeg[i]; }
	SSTAT(Svanum,va.num());
}
#endif

// GetenvValue("STATFILES") -> store data values in files.

class Stat {
  public:
	Stat(const char* pname=0, int pprint=0, int isstatic=0);
	~Stat();
	void setname(const char* pname);
	void setprint(int pprint);
	void zero();
	void terminate();
	double enter(double f);	// ret same f
	Stat& operator+=(double f);
	Stat& operator+=(const Stat& st);
	const char* name() const;
	int num() const;
	double smin() const;
	double smax() const;
	double avg() const;
	double var() const;
	double sdv() const;
	double sum() const;
	const char* string() const; // no leading name, no trailing '\n'
	const char* namestring() const; // operator<< uses namestring format
	friend ostream& operator<<(ostream& s, const Stat& st);
  private:
	const char* iname;
	int iprint;		// print statistics in destructor
	int n;
	double isum;
	double isum2;
	double imin;
	double imax;
	ofstream* fos;		// if getenv("STATFILES")
	void output(double f) const;
	DISABLECOPY(Stat);
};

#define STAT(Svar) Stat Svar(#Svar,1)
#define STATNP(Svar) Stat Svar(#Svar,0)
#define SSTAT(Svar,v) { static Stat Svar(#Svar,1,1); Svar+=v; }

//----------------------------------------------------------------------------

inline const char* Stat::name() const { return iname; }
inline int Stat::num() const { return n; }
inline double Stat::smin() const { return imin; }
inline double Stat::smax() const { return imax; }
inline double Stat::sdv() const { return mysqrt(var()); }
inline double Stat::avg() const { return MESHassertw1(n)?0:isum/n; }
inline double Stat::sum() const { return isum; }
inline void Stat::setprint(int pprint) { iprint=pprint; }

inline Stat& Stat::operator+=(double f)
{
	n++; isum+=f; isum2+=f*f;
	if (f<imin) imin=f;
	if (f>imax) imax=f;
	if (fos) output(f);
	return *this;
}

inline double Stat::enter(double f)
{
	*this+=f;
	return f;
}

#endif

