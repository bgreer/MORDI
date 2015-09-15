using namespace std;
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <sys/time.h>
#include <ctime>
// MKL:
#include "mkl.h"
//#define MKL_Complex16 std::complex<double>

#include "settings.h"
#include "grid.h"
#include "measurement_set.h"
#include "inversion.h"
#include "iterate_params.h"

settings glob_set;

double getTime()
{
	double ret;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	ret = (tv.tv_sec) + tv.tv_usec*1e-6;
	return ret;
}

int main (int argc, char *argv[])
{
	int myid, nproc;
	int itapod, itpad, itreg1, itreg2, itdepth, itreg1z, itreg2z;
	int itsz, itszm, itsh, itshm;
	int ii, nummodes;
	double starttime, endtime;
	measurement_set mset;
	inversion inv;
	iterate_params looper;
	grid<double> grid1;
	grid<double> grid2;
	bool target_loaded = false;

	// set up MPI
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	// timer
	starttime = getTime();

	// parse global settings
	glob_set.parse(argc, argv); // set up and parse command-line args
	glob_set.setParallel(myid,nproc); // serial
	looper.init();
	// force number of threads in openmp
	omp_set_dynamic(0);
	omp_set_num_threads(glob_set.nthreads);
	glob_set.printDetails();


	mset.loadMeasurements();
	// loads whole kernels and stores them in real space, no padding
	// parallelization to distribute whole kernels among processes
	mset.loadKernelSet(glob_set.kernel_set_fname, glob_set.kernel_z_fname);

	// set padding size, does not reallocate memory yet
	// just sets some integers for later
	mset.setPaddedSize();
	// this requires knowing the padding size
	mset.createCovar();
	mset.transformKernels();
	mset.computeOverlap();
	mset.clearTransform();
	mset.transformMeasurements();
	// remove modes that don't have a kernel or don't have measurements
	inv.init(mset.nx, mset.ny, mset.nz_kers, mset.set.size());

	do
	{
		// do things based on which parameter changed last
		// answer can be -1, which means this is the first iter
		if (looper.keyparam >= 10)
		{
			mset.setPaddedSize();
			mset.createCovar();
			mset.transformKernels();
			mset.computeOverlap();
		}
		if (looper.keyparam >= 9)
		{
			mset.clearTransform();
			mset.transformMeasurements();
			inv.init(mset.nx, mset.ny, mset.nz_kers, mset.set.size());
		}

		// do the actual inversion
		if (glob_set.myid==0) cout << "<<INVERSION>>\n";

		// pick a target function
		if (glob_set.target_fname!="" && !target_loaded)
		{
			// load target from file, same format as avgker
			inv.loadTarget();
			target_loaded = true;
		} else {
			// create Gaussian target function
			inv.setTarget(
				max(glob_set.currsigmah*glob_set.currdepth,glob_set.currsigmah_min),
				max(glob_set.currsigmaz*glob_set.currdepth,glob_set.currsigmaz_min),
				&mset);
		}
		// solve for the a coefficients
		inv.solveCoefs(&mset);
		if (glob_set.coefs_save_fname != "") inv.saveCoefs(&mset);
		// create a solution by applying coefficients to measurements
		inv.clearVelocity();
		inv.computeVelocity(&mset);
		inv.computeError(&mset);
		inv.outputVelocity();
		if (glob_set.avgker_fname!="") inv.outputAvgker(&mset);

	} while (looper.nextParams());



	fftw_cleanup();
	endtime = getTime();
	if (glob_set.myid==0)
		cout << "Total Runtime: " << (endtime-starttime) << " seconds."<<endl;
	MPI_Finalize();

	return 0;
}
