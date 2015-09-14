#ifndef MEASUREMENT_SET_H
#define MEASUREMENT_SET_H
#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <fftw3.h>
#include "omp.h"
#include "mpi.h"
#include "parallel.h"
#include "grid.h"
#include "mode.h"

using namespace std;

extern settings glob_set;

// keeps the set of measurements AND their kernels in ordered vectors

class measurement_set
{
public:
	vector<mode*> set;
	grid<complex<double> > **overlap; // mode, mode
	grid<complex<double> > covar;

	double *depth, *dz;

	int numlons, numlats;
	int offsetx, offsety;
	int nx, ny;
	int nx_kers, ny_kers, nz_kers;
	double dx;
	bool alloc_ft;

	measurement_set ();
	void transformMeasurements ();
	void clearTransform ();
	void computeOverlap ();
	void transformKernels (bool delete_real = true);
	void transformSingleKernel (mode *m, grid<complex<double> > *data_ft);
	void createCovar ();
	void setPaddedSize ();
	void loadMeasurements ();
	int loadKernelSet (string fname, string zfname);

};
#endif
