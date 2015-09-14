
extern settings glob_set;

// contains both measurement grid and kernel

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

	mode (int k, int n)
	{
		kval = k;
		nval = n;
		has_sker = false;
		proc_owner = -1;
		err = 0.0;
	}

	~mode ()
	{
		vx.dealloc();
		vy.dealloc();
	}

	void allocateRealKernel (int x, int y, int z)
	{
		int iz;
		sker = new grid<double> [z];
		for (iz=0; iz<z; iz++)
			sker[iz].init(x, y, false);
	}

	void allocateFTKernel (int x, int y, int z, bool parallel)
	{
		int iz, iy;
		sker_ft = new grid<complex<double> > [z];
		for (iz=0; iz<z; iz++)
			sker_ft[iz].init(x,y,parallel);
	}

	void deallocRealKernel (int z)
	{
		int iz;
		for (iz=0; iz<z; iz++)
			sker[iz].dealloc();
		delete [] sker;
	}

};
