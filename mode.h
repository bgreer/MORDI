#ifndef MODE_H
#define MODE_H
#include <complex>
#include "grid.h"

// Contains the measurements and kernels for all disk locations
// of a single wave mode. The kernels are translationally invariant,
// so you just need one kernel. The measurements are held in a grid

class mode
{
public:
	grid<double> vx, vy; // only held on proc0
	grid<complex<double> > vx_ft, vy_ft;
	bool has_sker;
	grid<double> *sker; // sensitivity kernel, exists on one processor
	grid<complex<double> > *sker_ft; // dimensions dimx_padded, dimy_padded, dimz_padded?
	int kval, nval;
	double err;
	// for parallelization, only one processor keeps the real-space data
	int proc_owner;
	// everyone holds the ft data together

	mode (int k, int n);

	~mode ()
	{
		vx.dealloc();
		vy.dealloc();
	}

	void allocateRealKernel (int x, int y, int z);
	void allocateFTKernel (int x, int y, int z, bool parallel);
	void deallocRealKernel (int z);

};
#endif
