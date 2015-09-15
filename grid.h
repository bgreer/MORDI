#ifndef GRID_H
#define GRID_H
// the purpose of a grid is to make it easier to parallelize things
// since the code splits up problems by wavenumber, 
// it's useful to hold an entire grid of data in a distributed way

// The grid itself is a construct that keeps information about where the grid
// data is held. It can be held on a single processor (parallel=false), or
// it can be held uniformly across every processor (parallel=true). 

// There are some helper functions in parallel.h that operate like MPI calls
// on the grids. You can send/recv grids, and distribute/collect grids.

// If you are unfamiliar with C++, the template is so that you can make grids
// of any data type, like a grid of floats or a grid of ints. The compiler 
// creates a specific implementation for float or int at compile time if
// needed elsewhere in the code.
#include <cstring>
#include "settings.h"

using namespace std;

extern settings glob_set; // contains parallel params

template<class Type> class grid
{
public:
	Type **data;
	int nx_local, ny_local; // my local size
	int nx, ny; // total, distributed size
	int x0, x1, y0, y1; // inclusive bounds in each direction
	bool alloc, parallel;

	// CONSTRUCTORS //
	grid ();
	grid (int x, int y, bool doparallel);
	bool contains (int x, int y);
	Type get (int x, int y);
	void set (int x, int y, Type val);
	Type* operator [] (int y);
	void init (int x, int y, bool doparallel);
	double norm(double val);
	void print(string fname, int mode);
	void dealloc();
	~grid ();
};
#include "grid.tcc"
#endif
