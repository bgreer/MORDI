#include <iostream>
#include <fstream>
#include <complex>
#include <cmath>
#include <fftw3.h>
#include "mkl.h"
#include "omp.h"
#include "mpi.h"
#include "grid.h"
#include "measurement_set.h"
#include "settings.h"
#include "parallel.h"

using namespace std;

extern settings glob_set;

class inversion
{
public:
	int nx, ny, nz, nn;
	grid<complex<double> > *coefs; // a coefficients, one kx,ky grid for each mode
	grid<complex<double> > target_ft;
	grid<complex<double> > *target_3d; // target function, one grid per depth
	grid<double> vx, vy; // for storing answer
	grid<complex<double> > coefs_real; // real-space representation of coefs, temp usage
	complex<double> err;
	grid<double> *avgker; // one grid for each depth
	double *target_z;
	bool has_init, alloc_vel;

	// Function prototypes:

	inversion ();
	void init (int x, int y, int z, int num);
	void loadTarget ();
	void setTarget (double widthx, double widthz, measurement_set *mset);
	void loadCoefs (measurement_set *mset);
	void saveCoefs (measurement_set *mset);
	void solveCoefs (measurement_set *mset);
	void computeError (measurement_set *mset);
	void computeVelocity (measurement_set *mset);
	void clearVelocity ();
	void outputVelocity ();
	void outputAvgker (measurement_set *mset);
};
