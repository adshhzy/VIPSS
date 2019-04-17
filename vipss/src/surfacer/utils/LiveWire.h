#pragma once
#include <math.h>
#include <stack>
#include <queue> 
#include <../src/utils/LiveWire/MyStack.h>
#include <../src/utils/LiveWire/MyQueue.h>
#include <../src/utils/LiveWire/MyList.h>
#include <../src/utils/LiveWire/MyCPoint.h>
#include <utils/Rn_Defs.H>
#include <utils/Mesh_Array.H>
#include <iostream>
#include <fstream>

class DistMap;
class MyImage;
class MyCPoint;

using namespace std;

// Object "LiveWire" is created by the GUI for each image opened by the user.
// A pointer to the corresponding image is passed along by the constructor (see below).
// "Segmentation" object allows efficient object segmentation via extraction of contours 
// and regions based on two very basic algorithms for image analysis.
//
// Contour extraction is activated when a user enters a sequence of seeds on the 
// desired object boundary (left mouse clicks). The GUI calls 'addSeed' 
// function after each mouse click. The seeds are automatically connected via "live-wire" 
// which is a short continuous "path" of pixels with high intensity contrast (note that
// high contrast normally occurs on object boundaries).
// A stack of points on a high-contrast path to the last entered seed is returned by 
// 'findPath(MyCPoint p)' function called by the GUI each time a mouse moves 
// to a new pixel 'p'. As a user clicks/adds a new seed (addSeed), the current live-wire is 
// permanently added to the list of contour points (m_contour). As the mouse moves on, live-wire
// returns high-contrast paths to the newly added seed.
// A right mouse click completes the contour: GUI calls 'addLastSeed(p)' function
// that adds current live-wire to the contour and also closes a contour by 
// computing a live-wire path between the contour's first seed and the last one (point p).
//
// Region extraction is activated by "growRegion(MyCPoint seed)" function called by the GUI when 
// a user clicks on a point in the desired region while holding button 'r' down on the keyboard.
// This function returns a mask (2D binary array of values 0 and 1) for the region
// that was grown from the specified seed by adding neighboring points with small intensity 
// difference. This mask (2D array) has the same size as the image.

class LiveWire
{
public:
	LiveWire(void);
	LiveWire(int rows, int cols, float** imgData); //Initializes a contour, links with an image
	~LiveWire(void);    
	bool addSeed(R2Pt pt); //Called by the GUI when a user makes a left or right mouse click
	R2Pt removeLastSeed();
	bool addLastSeed(R2Pt pt); // Called by the GUI when a user makes a right mouse click
	                            // This function 'loops' the contour
	void clear(); // Clears all contour points and seeds (called when a user presses "DEL")

	// Allocates, initializes, and returns a path (stack) of points from p to 
	// the last entereted seed, if any. GUI calls this function as a mouse moves
	// to a new point p.
	MyStack<MyCPoint>* findPath(MyCPoint p);	// for internal calls
	Array<R2Pt> findPath(R2Pt p);				// for external calls

	// This function is activated by a left (reset=false) or right (reset=true) mouse clicks 
	// when "r" key is down. It "grows" a region startind from the specified seed and returns
	// it via a binary "mask" of the same size as the image. The region grows by adding
	// to each pixel p already in the region all neighbors q such that |Ip-Iq|<T where
	// T is some fixed threshold. If 'reset' flag is false than the algorithm should
	// add a new region to existing set of masked pixels. Otherwise, the mask should
	// be cleared-up (all values reset to zero) first.
	void growRegion(char ** mask, R2Pt region_seed, bool reset = false);
	
	// Returns all points in the contour. GUI draws these pixels in green color
	// on the top of the image each time a display is refreshed
	inline const MyList<MyCPoint> * getContour() {return &m_contour;}
	// These allow for external access of contour since MyList and MyCPoint are made just for LiveWire class
	inline const R2Pt getContourPt(int i) {return R2Pt(m_contour.retrieve(i).x, m_contour.retrieve(i).y);}
	inline const int getContourSize() {return m_contour.getLength();}
	inline const int getLastSeedPos() {return m_iLastSeedPos;}
	//inline const R2Pt getLastSeed() { return R2Pt(m_seeds.retrieveLast().x,m_seeds.retrieveLast().y); }

	// Returns the status of the contour (Is it closed already or not?)
	inline bool isClosed() {return m_loop;}
	
	// Just to check if there are any points in the contour already
	//inline bool empty() {return m_contour.isEmpty();}
	inline bool empty() {return m_seeds.isEmpty();}

private:

	// Internal method that allows to "precompute" live-wire paths to the last seed
	// from any image pixel. Is called inside "addSeed" for each new mouse click.
	void computePaths(MyCPoint seed);

	MyList<MyCPoint> m_contour; // stores permanent contour points
	MyList<MyCPoint> m_seeds; // stored all enetered seed points
	MyImage * m_im; // stores a pointer to image data
	bool m_loop; // a closed contour flag
	int m_iLastSeedPos; //posistion of the last seed in the contour

	double ** m_penalty;   // A 2D table of individual pixel penalties 
	                       // Low penalty m_penalty[x][y] means high image contrast
	                       // at image pixel (x,y).

	DistMap * m_dist_map;  // Distance Map object for computing and
	                       // storing "live-wire" paths to the last seed
						   // In particular, for each pixel/node p it should store 
	                       // a "parent" neighboring node connecting p to the last seed 
	                       // along a "good" (high contrast) path. m_dist_path should 
	                       // also store a distance from 'p' to the sees along that path 
	                       // (distance = accumulated sum of penalties of nodes along the path).
};
