
//#include "grid.h"

template <class Type>
grid<Type>::~grid ()
{
	dealloc();
}

template <class Type>
grid<Type>::grid ()
{
	alloc = false;
}

template <class Type>
grid<Type>::grid (int x, int y, bool doparallel)
{
	init(x,y,doparallel);
}

// given global location, return true if exists in this particular grid
template <class Type>
bool grid<Type>::contains (int x, int y)
{
	if (x >= x0 && x <= x1 && y >= y0 && y <= y1) return alloc;
	return false;
}


// direct array access ([][]) permits out-of-bounds and such
// so use this instead
template <class Type>
Type grid<Type>::get (int x, int y) // global coords
{
	if (!contains(x,y)) return 0.0; // TODO: better error response
	return data[y-y0][x-x0]; // shift to local grid
}

template <class Type>
void grid<Type>::set (int x, int y, Type val)
{
	if (!contains(x,y)) return;
	data[y-y0][x-x0] = val;
}

// things that are bad ideas:
template <class Type>
Type* grid<Type>::operator [] (int y)
	{
		if (y < y0 || y > y1) return NULL;
		return data[y-y0];
	}

template <class Type>
void grid<Type>::init (int x, int y, bool doparallel)
	{
		int ii, me, num;

		nx = x;
		ny = y;

		parallel = doparallel;

		if (doparallel)
		{
			me = glob_set.myid;
			num = glob_set.nproc;
		} else {
			me = 0;
			num = 1;
		}

		// only allow a certain kind of parallelization
		if (num > ny)
		{
			cout << "ERROR: too many processes for proper parallelization of grid.\n";
			return;
		}

		// divide up grid
		x0 = 0;
		x1 = nx-1;
		nx_local = x1-x0+1;
		y0 = me*ny/num;
		y1 = (me+1)*ny/num - 1;
		ny_local = y1-y0+1;

		// allocate local space
		data = new Type* [ny_local];
		for (ii=0; ii<ny_local; ii++)
		{
			data[ii] = new Type [nx_local];
			memset(data[ii], 0x00, nx_local*sizeof(Type));
		}

		alloc = true;
	}

template <class Type>
double grid<Type>::norm(double val) {return val;}

	// mode = 0, print value
	// mode = 1, print norm
template <class Type>
void grid<Type>::print(string fname, int mode)
	{
		int ii, ix, iy;
		ofstream file;


		if (parallel)
		{
			for (ii=0; ii<glob_set.nproc; ii++)
			{
				MPI_Barrier(MPI_COMM_WORLD);
				if (ii==glob_set.myid)
				{
					file.open(fname.c_str(), fstream::out | fstream::app);
					for (ix=x0; ix<=x1; ix++)
					{
						for (iy=y0; iy<=y1; iy++)
						{
							if (mode==0)
								file << ix << " " << iy << " " << data[iy-y0][ix-x0] << endl;
							else
								file << ix << " " << iy << " " << std::norm(data[iy-y0][ix-x0]) << endl;
						}
					}
					file.close();
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}
		} else {
			for (ii=0; ii<glob_set.nproc; ii++)
			{
				if (ii==glob_set.myid)
				{
					file.open(fname.c_str(), fstream::out | fstream::app);
					for (ix=0; ix<nx; ix++)
					{
						for (iy=0; iy<ny; iy++)
						{
							if (mode==0)
								file << ix << " " << iy << " " << data[iy][ix] << endl;
							else
								file << ix << " " << iy << " " << std::norm(data[iy][ix]) << endl;
						}
					}
					file.close();
				}
			}
		}
	}

template <class Type>
void grid<Type>::dealloc()
	{
		int ii;
		if (alloc)
		{
			for (ii=0; ii<ny_local; ii++)
				delete [] data[ii];
			delete [] data;
			alloc = false;
		}
	}
	

