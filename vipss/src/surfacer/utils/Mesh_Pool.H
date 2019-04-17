// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1993 Hugues H. Hoppe; All rights reserved.

#ifndef Pool_h
#define Pool_h

// See notes in Pool.C

#define POOLALLOCATION(class) \
void* operator new(size_t s) \
{ return s==sizeof(class)?pool.alloc():pool.specialalloc(s); } \
void operator delete(void* p, size_t s) \
{ s==sizeof(class)?pool.free(p):pool.specialfree(p,s); } \
static Pool pool

#define ALLOCATEPOOL(class) Pool class::pool(sizeof(class),#class)

class Pool {
  public:
	Pool(unsigned pesize, const char* pname);
	~Pool();
	void* alloc();
	void free(void* p);
	void* specialalloc(size_t s);
	void specialfree(void* p, size_t s);
  private:
	struct Link { Link* next; };
	struct Chunk { Chunk* next; };
	int esize;
	const char* iname;
	Link* h;
	Chunk* chunkh;
	int nalloc;
	int offset;
	void grow();
	DISABLECOPY(Pool);
};

//----------------------------------------------------------------------

inline void* Pool::alloc()
{
	if (!h) grow();
	Link* p=h; h=p->next; return p;
}

inline void Pool::free(void* pp)
{
	Link* p=(Link*)pp;
	p->next=h; h=p;
}

#endif

