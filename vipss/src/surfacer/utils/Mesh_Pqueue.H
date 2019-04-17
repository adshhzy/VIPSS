// This may look like C code, but it is really -*- C++ -*-
// Copyright (c) 1992 Hugues H. Hoppe; All rights reserved.

#ifndef Pqueue_h
#define Pqueue_h

#include <utils/Mesh_Map.H>

#include <utils/Mesh_Pool.H>


/** \class BPqueue Mesh_PQueue.H utils/Mesh_PQueue.H

  \ingroup MeshesQueue

  \brief Priority queue which can be resized

  Base class for all the other queue types

  This should really be a templated class, but instead it uses a void
  pointer (Univ) to allow you to attach data to the queue items.

  Two methods of use: either put the elements all in (enterUnsorted) then sort,
  or put each element in, sorting as you go (slower).
  */
class BPqueue {
 public:
     /// Initialize the queue; default size is zero
    BPqueue(int size=0);
    ///
    ~BPqueue();
    /// Empty the queue
    void clear();
    /** Add an element to the queue.  
      @param e is a void pointer - can be a pointer to a class, or an integer
      @param pri Priority. Must be greater than or equal to zero
      */
    void enter(Univ e, double pri); 
    /// Clear all elements out of the queue
    int empty() const;
    /// Number of elements in the queue
    int num() const;
    /// Get the minimum priority element without popping it off (ret e, else die)
    Univ bmin() const;           
    /// Return the current minimum priority (die if no elems in queue)
    double minpriority() const;  
    /// Pop the minimum priority element off (die if no elements)
    Univ removemin();           
    /// Just put it in anywhere
    void enterUnsorted(Univ e, double pri);
    /// Sort
    void sort();
 protected:
    struct Node { Univ e; double pri; };
    int isize;
    int inum;
    Node* ar;
    // {ar,isize,inum} could be made Array<Node> but for efficiency
    void resize(int newsize);
    void nswitch(int n1, int n2);
    void adjust(int n, int up, int down);
 private:
    DISABLECOPY(BPqueue);
};

/** \class BHPqueue Mesh_PQueue.H utils/Mesh_PQueue.H

  \ingroup MeshesQueue

  \brief Priority queue with a hash table

  This should really be a templated class, but instead it uses a void
  pointer (Univ) to allow you to attach data to the queue items. Pqueue,
  below, wraps a template around this class.

  Essentially adds a hash table to the priority queue so it is
  easy to find elements.
  */
class BHPqueue : private BPqueue {
 public:
     /// default size is zero
    BHPqueue(int size=0);
    ///
    ~BHPqueue();
    /// Enter element in the queue, pri>=0 always
    void enter(Univ e, double pri); // pri>=0 always
    int empty() const { return BPqueue::empty(); }
    int num() const { return BPqueue::num(); }
    Univ bmin() const { return BPqueue::bmin(); }
    double minpriority() const { return BPqueue::minpriority(); }
    /// Pop the top one off
    Univ removemin();
    /// Completely empty queue
    void clear();
    /// Is this element in the queue? Uses pointer comparison
    int contains(Univ e) const;
    /// Get the specific element's priorty, returns < 0 if not in queue
    double retrieve(Univ e) const; 
    /// Remove the specific element, returns < 0 if not in queue
    double remove(Univ e);         // ret pri or <0
    /// Change the element's priority, returns < 0 if not in queue. pri should be > 0
    double update(Univ e, double pri); // ret prevpri or <0
    /// Enter the element or change the element's priority if in queue already, returns < 0 if not in queue. pri should be > 0
    double enterupdate(Univ e, double pri); // ret prevpri or <0
    /// Just put it in anywhere
    void enterUnsorted(Univ e, double pri);
    ///
    void sort();
 private:
    MMap<Univ,int> m;            // element -> index in array
    void nswitch(int n1, int n2);
    void adjust(int n, int up, int down);
    int find(Univ e) const;     // ret -1 if not there
};

//----------------------------------------------------------------------------

inline int BPqueue::empty() const { return !inum; }
inline int BPqueue::num() const { return inum; }
inline Univ BPqueue::bmin() const { MESHassertx(inum); return ar[0].e; }
inline double BPqueue::minpriority() const { MESHassertx(inum); return ar[0].pri; }

inline void BPqueue::enterUnsorted(Univ e, double pri)
{
    MESHassertx(pri>=0);
    if (inum==isize) resize(isize*2);
    ar[inum].e=e;
    ar[inum].pri=pri;
    inum++;
}

inline void BHPqueue::enterUnsorted(Univ e, double pri)
{
    m.enter(e,inum);
    BPqueue::enterUnsorted(e,pri);
}

//----------------------------------------------------------------------------

/** \class Pqueue Mesh_PQueue.H utils/Mesh_PQueue.H

  \ingroup MeshesQueue

  \brief Priority queue which can be resized

  Basically a template wrapper around BPqueue

  Should not use this if template type is not a pointer or basic
  type because it will do a lot of copy constructor calls.

  Two methods of use: either put the elements all in (enterUnsorted) then sort,
  or put each element in, sorting as you go (slower).
  */
template<class T>
class Pqueue : public BPqueue {
 public:
     ///
    inline Pqueue() { }
    ///
    inline virtual ~Pqueue() { }
    /// Templated enter function
    inline virtual void enter(T e, double pri) { BPqueue::enter(Conv<T>::e(e),pri); }
    /// Get the element on the top of the queue (does not pop)
    inline T bmin() const { return Conv<T>::d(BPqueue::bmin()); }
    /// Pop the element off the queue
    inline virtual T removemin() { return Conv<T>::d(BPqueue::removemin()); }
    /// Put on queue without sorting
    inline virtual void enterUnsorted(T e, double pri)
    { BPqueue::enterUnsorted(Conv<T>::e(e),pri); }
};

/** \class HPqueue Mesh_PQueue.H utils/Mesh_PQueue.H

  \ingroup MeshesQueue

  \brief Priority queue with a hash table which can be resized

  Basically a template wrapper around BHPqueue.

  Should not use this if template type is not a pointer or basic
  type because it will do a lot of copy constructor calls.
  */
template<class T>
class HPqueue : public BHPqueue {
 public:
     ///
    inline HPqueue() { }
    ///
    inline ~HPqueue() { }
    /// Enter element in the queue, pri>=0 always
    inline virtual void enter(T e, double pri) { BHPqueue::enter(Conv<T>::e(e),pri); }
    /// Get the smallest element without popping Fail if no elements.
    inline T bmin() const { return Conv<T>::d(BHPqueue::bmin()); }
    /// Pop the smallest element. Fail if no elements.
    inline virtual T removemin() { return Conv<T>::d(BHPqueue::removemin()); }
    /// Put in the queue without sorting
    inline virtual void enterUnsorted(T e, double pri)
    { BHPqueue::enterUnsorted(Conv<T>::e(e),pri); }
    /// Does this queue contain this element? Return non-zero if true
    inline int contains(T e) const
    { return BHPqueue::contains(Conv<T>::e(e)); }
    /// Get this element's priority value. Returns < 0 if not in queue
    inline double retrieve(T e) const
    { return BHPqueue::retrieve(Conv<T>::e(e)); }
    /// Remove element from queue. Returns < 0 if not in queue
    inline double remove(T e)
    { return BHPqueue::remove(Conv<T>::e(e)); }
    /// Change the priority value. Returns old priority value, or < 0 if not in queue
    inline double update(T e, double pri)
    { return BHPqueue::update(Conv<T>::e(e),pri); }
    /// Enter into the queue; returns old priority value if in queue already, < 0 if not in queue
    inline virtual double enterupdate(T e, double pri)
    { return BHPqueue::enterupdate(Conv<T>::e(e),pri); }
};

#endif

