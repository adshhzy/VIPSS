// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Stack_h
#define Stack_h

#if 0
{
	Stack<Edge> stack;
	while (!stack.empty()) delete stack.pop();
	ForStack(stack,Point*,p) { consider(p); } EndFor;
}
#endif

#include <utils/Mesh_Pool.H>


class BStack {
  public:
	BStack();
	~BStack();
	void clear();
	void push(Univ e);
	Univ pop();		// ret e, else die
	Univ top() const;	// ret e, else die
	int empty() const;
	int height() const;	// slow O(n)
	int contains(Univ e) const; // slow O(n)
	int remove(Univ e);	    // slow O(n), ret: was_there
	static const BStack EMPTY;
	POOLALLOCATION(BStack);
  private:
    friend class BStackIter;
	struct Node {
		POOLALLOCATION(Node);
		Univ data;
		Node* next;
	};
	Node* itop;
	DISABLECOPY(BStack);
};

class BStackIter {
  public:
	BStackIter(const BStack& s);
	~BStackIter();
	operator void*() const;
	Univ next();
	Univ operator()() const;
	void reinit(const BStack& s);
  private:
	const BStack::Node* n;
	// operator= is safe
};

//----------------------------------------------------------------------------

inline BStackIter::BStackIter(const BStack& s) : n(s.itop) { }
inline BStackIter::~BStackIter() { }
inline Univ BStackIter::next() { Univ d=n->data; n=n->next; return d; }
inline BStackIter::operator void*() const { return n?(void*)1:0; }
inline Univ BStackIter::operator()() const { return n->data; }
inline void BStackIter::reinit(const BStack& s) { n=s.itop; }


inline BStack::BStack() : itop(0) { }
inline int BStack::empty() const { return !itop; }
inline Univ BStack::top() const { return itop->data; }

inline void BStack::push(Univ e)
{
	Node* n=new Node;
	n->data=e;
	n->next=itop;
	itop=n;
}

inline Univ BStack::pop()
{
	Node* n=MESHassertv(itop);
	Univ e=n->data;
	itop=n->next;
	delete n;
	return e;
}

//----------------------------------------------------------------------------

template<class T> class StackIter;

template<class T>
class Stack : public BStack {
  public:
	inline Stack() { }
	inline ~Stack() { }
	inline void push(T e) { BStack::push(Conv<T>::e(e)); }
	inline T pop() { return Conv<T>::d(BStack::pop()); }
	inline T top() { return Conv<T>::d(BStack::top()); }
	inline int contains(T e) const
	{ return BStack::contains(Conv<T>::e(e)); }
	inline int remove(T e)
	{ return BStack::remove(Conv<T>::e(e)); }
// 	typedef StackIter<T> Iter;
};

template<class T>
class StackIter : public BStackIter {
  public:
	inline StackIter(const Stack<T>& s) : BStackIter(s) { }
	inline ~StackIter() { }
	inline T next() { return Conv<T>::d(BStackIter::next()); }
	inline T operator()() const
	{return Conv<T>::d(BStackIter::operator()());}
	inline void reinit(const Stack<T>& s) { BStackIter::reinit(s); }
};

#define ForStack(S,T,V) { for (StackIter<T> zz(S);zz;zz.next()) { T V=zz();

#endif

