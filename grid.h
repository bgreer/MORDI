#ifndef GRID_H
#define GRID_H
// the purpose of a grid is to make it easier to parallelize things
// since the code splits up problems by wavenumber, 
// it's useful to hold an entire grid of data in a distributed way
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
