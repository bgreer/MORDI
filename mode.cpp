
#include "mode.h"

mode::mode (int k, int n)
	{
		kval = k;
		nval = n;
		has_sker = false;
		proc_owner = -1;
		err = 0.0;
	}


void mode::allocateRealKernel (int x, int y, int z)
	{
		int iz;
		sker = new grid<double> [z];
		for (iz=0; iz<z; iz++)
			sker[iz].init(x, y, false);
	}

void mode::allocateFTKernel (int x, int y, int z, bool parallel)
	{
		int iz, iy;
		sker_ft = new grid<complex<double> > [z];
		for (iz=0; iz<z; iz++)
			sker_ft[iz].init(x,y,parallel);
	}

void mode::deallocRealKernel (int z)
	{
		int iz;
		for (iz=0; iz<z; iz++)
			sker[iz].dealloc();
		delete [] sker;
	}
