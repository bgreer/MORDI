using namespace std;
#include <iostream>
#include <typeinfo>
#include <iomanip>
#include <sys/time.h>
#include <ctime>
#include "omp.h"
#include "mpi.h"
// MKL:
#include "mkl.h"
//#define MKL_Complex16 std::complex<double>
#include "settings.h"
#include "kernel_set.h"
#include "measurements.h"
#include "measurement_set.h"
#include "inversion.h"
#include "parallel.h"

settings glob_set;

double getTime()
{
	double ret;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	ret = (tv.tv_sec-1422700000) + tv.tv_usec*1e-6;
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
	// force number of threads in openmp
	omp_set_dynamic(0);
	omp_set_num_threads(glob_set.nthreads);
	glob_set.printDetails();


	mset.loadMeasurements();
	// loads whole kernels and stores them in real space, no padding
	// parallelization to distribute whole kernels among processes
	mset.loadKernelSet(glob_set.kernel_set_fname, glob_set.kernel_z_fname);

// BEGIN LOOPS
	for (itpad=0; itpad<glob_set.padding.size(); itpad++) { // padding loop
		glob_set.currpadding = glob_set.padding[itpad];
		if (glob_set.myid==0) cout << "ITERATION: Padding ("<<fixed<<glob_set.currpadding<<")\n";

		mset.setPaddedSize();
		mset.createCovar();
		mset.transformKernels();
		mset.computeOverlap();

	for (itapod=0; itapod<glob_set.apodization.size(); itapod++) { // apodization loop
		glob_set.currapodization = glob_set.apodization[itapod];
		if (glob_set.myid==0) cout << "ITERATION: Data Apodization ("<<fixed<<glob_set.currapodization<<")\n";

		mset.clearTransform();
		mset.transformMeasurements();
		// remove modes that don't have a kernel or don't have measurements
		inv.init(mset.nx, mset.ny, mset.nz_kers, mset.set.size());

	for (itreg1z=0; itreg1z<glob_set.reg1z.size(); itreg1z++) { // regularization 1z loop
		glob_set.currreg1z = glob_set.reg1z[itreg1z];
		if (glob_set.myid==0) cout << "ITERATION: Lambda-z ("<<fixed<<glob_set.currreg1z<<")\n";
	for (itreg2z=0; itreg2z<glob_set.reg2z.size(); itreg2z++) { // regularization 2z loop
		glob_set.currreg2z = glob_set.reg2z[itreg2z];
		if (glob_set.myid==0) cout << "ITERATION: Mu-z ("<<fixed<<glob_set.currreg2z<<")\n";
	for (itreg1=0; itreg1<glob_set.reg1.size(); itreg1++) { // regularization 1 loop
		glob_set.currreg1 = glob_set.reg1[itreg1];
		if (glob_set.myid==0) cout << "ITERATION: Lambda ("<<fixed<<glob_set.currreg1<<")\n";
	for (itreg2=0; itreg2<glob_set.reg2.size(); itreg2++) { // regularization 2 loop
		glob_set.currreg2 = glob_set.reg2[itreg2];
		if (glob_set.myid==0) cout << "ITERATION: Mu ("<<fixed<<glob_set.currreg2<<")\n";
	for (itsz=0; itsz<glob_set.sigmaz.size(); itsz++) { // sigmaz loop
		glob_set.currsigmaz = glob_set.sigmaz[itsz];
		if (glob_set.myid==0) cout << "ITERATION: SigmaZ {" << glob_set.currsigmaz<<")\n";
	for (itszm=0; itszm<glob_set.sigmaz_min.size(); itszm++) { // sigmaz_min loop
		glob_set.currsigmaz_min = glob_set.sigmaz_min[itszm];
		if (glob_set.myid==0) cout << "ITERATION: SigmaZ-min {" << glob_set.currsigmaz_min<<")\n";
	for (itsh=0; itsh<glob_set.sigmah.size(); itsh++) { // sigmah loop
		glob_set.currsigmah = glob_set.sigmah[itsh];
		if (glob_set.myid==0) cout << "ITERATION: SigmaH {" << glob_set.currsigmah<<")\n";
	for (itshm=0; itshm<glob_set.sigmah_min.size(); itshm++) { // sigmah_min loop
		glob_set.currsigmah_min = glob_set.sigmah_min[itshm];
		if (glob_set.myid==0) cout << "ITERATION: SigmaH-min {" << glob_set.currsigmah_min<<")\n";
	for (itdepth=0; itdepth<glob_set.depths.size(); itdepth++) { // depth loop
		glob_set.currdepth = glob_set.depths[itdepth];
		if (glob_set.myid==0) cout << "ITERATION: Depth ("<<fixed<<glob_set.currdepth<<")\n";


		// INVERSION
		if (glob_set.myid==0) cout << "<<INVERSION>>\n";

		if (glob_set.target_fname!="" && !target_loaded)
		{
			inv.loadTarget();
			target_loaded = true;
		} else
			inv.setTarget(
					max(glob_set.currsigmah*glob_set.currdepth,glob_set.currsigmah_min),
					max(glob_set.currsigmaz*glob_set.currdepth,glob_set.currsigmaz_min),
					&mset);
		inv.solveCoefs(&mset);
		if (glob_set.coefs_save_fname != "") inv.saveCoefs(&mset);
		inv.clearVelocity();
		inv.computeVelocity(&mset);
		inv.computeError(&mset);
		inv.outputVelocity();
		if (glob_set.avgker_fname!="") inv.outputAvgker(&mset);

	} // depth loop
	} // sigmah-min loop
	} // sigmah loop
	} // sigmaz-min loop
	} // sigmaz loop
	} // mu loop
	} // lambda loop
	} // mu-z loop
	} // lambda-z loop
	} // padding loop
	} // apodization loop
// END LOOPS

	fftw_cleanup();
	endtime = getTime();
	if (glob_set.myid==0)
		cout << "Total Runtime: " << (endtime-starttime) << " seconds."<<endl;
	MPI_Finalize();

	return 0;
}
