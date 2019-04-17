// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Queue_h
#define Queue_h

#if 0
{
	Queue<Vertex> queuev;
	while (!queue.empty()) delete queue.dequeue();
	ForStack(stack,Point*,p) { consider(p); } EndFor;
}
#endif

#include <utils/Mesh_Pool.H>


class BQueue {
  public:
	BQueue();
	~BQueue();
	void clear();
	void enqueue(Univ e);	// "insertlast"
	Univ dequeue();		// ret e, else die
	Univ front() const;	// ret e, else die (only look)
	void insertfirst(Univ e);
	int empty() const;
	int length() const;	// slow O(n)
	int contains(Univ e) const; // slow O(n)
	void addtoend(BQueue& q);    // quick merge, clears q
	void addtofront(BQueue& q);  // quick merge, clears q
	static const BQueue EMPTY;
  private:
    friend class BQueueIter;
	struct Node {
		POOLALLOCATION(Node);
		Univ data;
		Node* next;
	};
	Node* ifront;
	Node* irear;
	DISABLECOPY(BQueue);
};

class BQueueIter {
  public:
	BQueueIter(const BQueue& s);
	virtual ~BQueueIter();
	operator void*() const;
	void next();
	Univ operator()() const;
  private:
	const BQueue::Node* n;
	// operator= is safe
};

//----------------------------------------------------------------------------

inline BQueueIter::BQueueIter(const BQueue& s) : n(s.ifront) { }
inline BQueueIter::~BQueueIter() { }
inline void BQueueIter::next() { n=n->next; }
inline BQueueIter::operator void*() const { return n?(void*)1:0; }
inline Univ BQueueIter::operator()() const { return n->data; }


inline BQueue::BQueue() : ifront(0), irear(0) { }
inline int BQueue::empty() const { return !ifront; }
inline Univ BQueue::front() const { return ifront->data; }

inline void BQueue::enqueue(Univ e)
{
	Node* n=new Node;
	n->data=e;
	n->next=0;
	if (irear) irear->next=n;
	else ifront=n;
	irear=n;
}

inline void BQueue::insertfirst(Univ e)
{
	Node* n=new Node;
	n->data=e;
	n->next=ifront;
	ifront=n;
	if (!irear) irear=n;
}

inline Univ BQueue::dequeue()
{
	Node* n=MESHassertv(ifront);
	Univ e=n->data;
	ifront=n->next;
	delete n;
	if (!ifront) irear=0;
	return e;
}

//----------------------------------------------------------------------------

template<class T> class QueueIter;

template<class T>
class Queue : public BQueue {
  public:
	Queue() { }
	~Queue() { }
	inline void enqueue(T e) { BQueue::enqueue(Conv<T>::e(e)); }
	inline T dequeue() { return Conv<T>::d(BQueue::dequeue()); }
	inline T front() const { return Conv<T>::d(BQueue::front()); }
	inline void insertfirst(T e) { BQueue::insertfirst(Conv<T>::e(e)); }
	inline int contains(T e) const
	{ return BQueue::contains(Conv<T>::e(e)); }
	inline void addtoend(Queue<T>& q) { BQueue::addtoend(q); }
	inline void addtofront(Queue<T>& q) { BQueue::addtofront(q); }
// 	typedef QueueIter<T> Iter;
};

template<class T>
class QueueIter : public BQueueIter {
  public:
	inline QueueIter(const Queue<T>& q) : BQueueIter(q) { }
	inline ~QueueIter() { }
	inline T operator()() const
	{return Conv<T>::d(BQueueIter::operator()());}
};

#define ForQueue(S,T,V) { for (QueueIter<T> zz(S);zz;zz.next()) { T V=zz();
								  
#endif

