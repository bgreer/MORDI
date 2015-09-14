#include "mode.h"

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

	measurement_set ()
	{
		alloc_ft = false;
	}
	void transformMeasurements ()
	{
		int ii, iz, ij, ik, ix, iy, count, ind, tind;
		int ind_src, ind_dest;
		grid<complex<double> > vx_ft, vy_ft;
		complex<double> *input, *output;
		double apod, mult, dist;
		fftw_plan plan;

		if (glob_set.myid==0) cout << "Transforming measurements.."<<endl;

		input = new complex<double> [nx*ny];
		output = new complex<double> [nx*ny];
		plan = fftw_plan_dft_2d(nx, ny, 
				reinterpret_cast<fftw_complex*>(input), 
				reinterpret_cast<fftw_complex*>(output), 
				FFTW_FORWARD, FFTW_ESTIMATE);
		// make space for a single whole map
		vx_ft.init(nx,ny,false);
		vy_ft.init(nx,ny,false);
		apod = glob_set.currapodization;

		// march through all modes
		tind = 0;
		for (ii=0; ii<set.size(); ii++)
		{
			// find next kernel to process
			while (set[tind]->proc_owner != glob_set.myid && tind < set.size())
				tind ++;
			if (tind < set.size())
			{
				// TRANSFORM X
				memset(input, 0x00, nx*ny*sizeof(complex<double>));
				for (ix=0; ix<numlons; ix++)
				{
					for (iy=0; iy<numlats; iy++)
					{
						mult = 1.0;
						if (glob_set.apod_circ)
						{
							// circular apodization
							dist = numlons*0.5*dx - 
								sqrt(pow(ix-numlons*0.5,2.0) + pow(iy-numlats*0.5,2.0))*dx;
							if (dist < apod) mult = pow(1.0-pow((dist-apod)/apod,2.0),2.0);
							if (dist < 0.0) mult = 0.0;
						} else {
							// rectangular apodization
							dist = min(ix, numlons-ix)*dx;
							if (dist < apod) mult *= pow(1.-pow((apod-dist)/apod,2.),2.);
							dist = min(iy, numlats-iy)*dx;
							if (dist < apod) mult *= pow(1.-pow((apod-dist)/apod,2.),2.);
						}
						ind_dest = (ix+offsetx)*ny + iy + offsety;
						input[ind_dest] = set[tind]->vx.get(ix,iy)*mult;
					}
				}
				fftw_execute(plan);
				for (ix=0; ix<nx; ix++)
				{
					for (iy=0; iy<ny; iy++)
					{
						ind_src = ix*ny + iy;
						vx_ft.set(ix,iy,output[ind_src]);
					}
				}
				// TRANSFORM Y
				memset(input, 0x00, nx*ny*sizeof(complex<double>));
				for (ix=0; ix<numlons; ix++)
				{
					for (iy=0; iy<numlats; iy++)
					{
						mult = 1.0;
						if (glob_set.apod_circ)
						{
							// circular apodization
							dist = numlons*0.5*dx - 
								sqrt(pow(ix-numlons*0.5,2.0) + pow(iy-numlats*0.5,2.0))*dx;
							if (dist < apod) mult = pow(1.0-pow((dist-apod)/apod,2.0),2.0);
							if (dist < 0.0) mult = 0.0;
						} else {
							// rectangular apodization
							dist = min(ix, numlons-ix)*dx;
							if (dist < apod) mult *= pow(1.-pow((apod-dist)/apod,2.),2.);
							dist = min(iy, numlats-iy)*dx;
							if (dist < apod) mult *= pow(1.-pow((apod-dist)/apod,2.),2.);
						}
						ind_dest = (ix+offsetx)*ny + iy + offsety;
						input[ind_dest] = set[tind]->vy.get(ix,iy)*mult;
					}
				}
				fftw_execute(plan);
				for (ix=0; ix<nx; ix++)
				{
					for (iy=0; iy<ny; iy++)
					{
						ind_src = ix*ny + iy;
						vy_ft.set(ix,iy,output[ind_src]);
					}
				}
				tind ++;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			// share
			for (ij=0; ij<glob_set.nproc; ij++)
			{
				// which mode did proc ij just process?
				// ind = (the ii-th mode belonging to proc ij)
				ind = -1;
				count = 0;
				for (ik=0; ik<set.size(); ik++)
				{
					if (set[ik]->proc_owner == ij) count ++;
					if (count == ii+1)
					{
						ind = ik;
						ik = set.size();
					}
				}

				// proc ij did indeed process mode ind
				if (ind >= 0)
				{
					set[ind]->vx_ft.init(nx,ny,true);
					set[ind]->vy_ft.init(nx,ny,true);
					distributeGrid(&(vx_ft),&(set[ind]->vx_ft),ij);
					distributeGrid(&(vy_ft),&(set[ind]->vy_ft),ij);
				}
				

			}
		}

		fftw_destroy_plan(plan);
		fftw_cleanup();
		if (glob_set.myid==0) cout << "Done transforming measurements.\n";
		alloc_ft = true;
	}

	void clearTransform ()
	{
		int ii;
		if (alloc_ft)
		{
			for (ii=0; ii<set.size(); ii++)
			{
				set[ii]->vx.dealloc();
				set[ii]->vy.dealloc();
			}
			alloc_ft = false;
		}
	}

	void computeOverlap ()
	{
		int ii,ij,ix,iy,iz;
		int x0,x1,y0,y1;
		complex<double> val, a, b, zero(0,0);
		grid<complex<double> > recv;
		double delta;

		// allocate space for overlap matrix
		overlap = new grid<complex<double> >* [set.size()];
		for (ii=0; ii<set.size(); ii++)
			overlap[ii] = new grid<complex<double> > [set.size()];

		MPI_Barrier(MPI_COMM_WORLD);

		// assume parallel division has been done the same everywhere
		x0 = 0;
		x1 = nx-1;
		y0 = glob_set.myid*ny/glob_set.nproc;
		y1 = (glob_set.myid+1)*ny/glob_set.nproc - 1;

		if (glob_set.myid==0) cout << "Computing Overlap.."<<endl;

		// everyone has some part of the ft of each kernel
#pragma omp parallel for private(ii,ij,ix,iy,iz,val,a,b) schedule(dynamic,1)
		for (ii=0; ii<set.size(); ii++)
		{
			for (ij=ii; ij<set.size(); ij++)
			{
				overlap[ii][ij].init(nx,ny,true);
				overlap[ij][ii].init(nx,ny,true);
				// compute using own part
				for (iy=y0; iy<=y1; iy++)
				{
					for (ix=x0; ix<=x1; ix++)
					{
						overlap[ii][ij].set(ix,iy,zero);
						for (iz=0; iz<nz_kers; iz++)
						{
							val = overlap[ii][ij].get(ix,iy);
							a = set[ii]->sker_ft[iz].get(ix,iy);
							b = set[ij]->sker_ft[iz].get(ix,iy);
							delta = dx * dx * dz[iz];
							overlap[ii][ij].set(ix,iy, val + (a*b*delta));
						}
						overlap[ij][ii].set(ix,iy,overlap[ii][ij].get(ix,iy));
					}
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (glob_set.myid==0) cout << "Done computing overlap."<<endl;
	}


	// given the grid size, interpolate each kernel and FT it
	// store result inside each kernel class
	// by default, clears out the real-space kernel and frees that memory
	void transformKernels (bool delete_real = true)
	{
		int ii, iz, ij, ik, count, ind, tind;
		grid<complex<double> > *data_ft; // dimensions dimx_padded, dimy_padded, dimz_padded?

		// everyone has a set of whole kernels
		// march through personal kernels
		// - create a whole transformed kernel in personal memory
		// - everyone takes turns
		// - - make room for distributed ftkernel
		// - - distribute data around
	
		if (glob_set.myid==0) cout << "Transforming kernels.."<<endl;

		// make space for a single whole kernel
		data_ft = new grid<complex<double> > [nz_kers];
		for (iz=0; iz<nz_kers; iz++)
			data_ft[iz].init(nx,ny,false);

		// march through all modes
		tind = 0;
		for (ii=0; ii<set.size(); ii++)
		{
			// find next kernel to process
			while (set[tind]->proc_owner != glob_set.myid && tind < set.size())
				tind ++;
			if (tind < set.size())
			{
				transformSingleKernel(set[tind], data_ft);
				if (delete_real)
					set[tind]->deallocRealKernel(nz_kers);
				tind ++;
			}
			MPI_Barrier(MPI_COMM_WORLD);
			// share
			for (ij=0; ij<glob_set.nproc; ij++)
			{
				// which mode did proc ij just process?
				// ind = (the ii-th mode belonging to proc ij)
				ind = -1;
				count = 0;
				for (ik=0; ik<set.size(); ik++)
				{
					if (set[ik]->proc_owner == ij) count ++;
					if (count == ii+1)
					{
						ind = ik;
						ik = set.size();
					}
				}

				// proc ij did indeed process mode ind
				if (ind >= 0)
				{
					set[ind]->allocateFTKernel(nx,ny,nz_kers,true);
					for (iz=0; iz<nz_kers; iz++)
						distributeGrid(&(data_ft[iz]),&(set[ind]->sker_ft[iz]),ij);
				}
				

			}
		}

		if (glob_set.myid==0) cout << "Done transforming kernels.\n";
	}

	void transformSingleKernel (mode *m, grid<complex<double> > *data_ft)
	{
#pragma omp parallel
	{
		int ix, iy, iz, index, x, y;
		double val, r, apod;
		complex<double> *input, *output;
		fftw_plan plan;
		// manual parallelism with openmp
		int startz, endz, myid;
		myid = omp_get_thread_num();
		startz = myid*nz_kers / glob_set.nthreads;
		endz = (myid+1)*nz_kers / glob_set.nthreads - 1;

		input = new complex<double> [nx*ny];
		output = new complex<double> [nx*ny];
		// openmp is not thread-safe when planning
#pragma omp critical
		plan = fftw_plan_dft_2d(nx, ny, 
				reinterpret_cast<fftw_complex*>(input), 
				reinterpret_cast<fftw_complex*>(output), 
				FFTW_FORWARD, FFTW_ESTIMATE);

		// march through each depth
		for (iz=startz; iz<=endz; iz++)
		{
			// interpolate existing kernel on to input array
			for (iy=0; iy<ny; iy++)
			{
				for (ix=0; ix<nx; ix++)
				{
					// positions relative to origin of original kernel grid
					x = ix - nx/2 + nx_kers/2;
					y = iy - ny/2 + ny_kers/2;

					// apodization
					r = sqrt(pow(x-nx_kers*0.5,2) + pow(y-ny_kers*0.5,2))
						/ (8./dx);

					apod = 1;
//					if (r > 0.85) apod = pow(1.0-pow((r-0.85)/(1.0-0.85),2.),2.);
					if (r > 0.9375) apod = pow(1.0-pow((r-0.9375)/(1.0-0.9375),2.),2.);
			//		apod = pow(cos(r*3.14159/2.0),2.0);
					if (r > 1.0) apod = 0.0;

					if (x >= 0 && x < nx_kers && y >= 0 && y < ny_kers)
						val = m->sker[iz].get(x,y)*apod;
					else
						val = 0.0;
					index = iy*nx + ix;
					input[index] = complex<double> (val,0.0);
				}
			}
			fftw_execute(plan);
			// unload into kernel
			for (iy=0; iy<ny; iy++)
			{
				for (ix=0; ix<nx; ix++)
				{
					index = iy*nx + ix;
					data_ft[iz].set(ix,iy,output[index]);
				}
			}
		}
		fftw_destroy_plan(plan);
		delete [] input;
		delete [] output;

	} // end openmp
		fftw_cleanup();
		// TODO: enforce normalization

	}

	// compute the mutual overlap of two apodization circles
	// in fourier space
	void createCovar ()
	{
		int ii, ix, iy, index;
		double x, y, cx, cy, r, val;
		const double apod_min = 0.9375;

		complex<double> *input, *output;
		fftw_plan plan, plan2;

		// make room for ffts
		input = new complex<double> [nx*ny];
		output = new complex<double> [nx*ny];
		// plan the fft
		plan = fftw_plan_dft_2d(nx, ny, 
				reinterpret_cast<fftw_complex*>(input), 
				reinterpret_cast<fftw_complex*>(output), 
				FFTW_FORWARD, FFTW_ESTIMATE);
		plan2 = fftw_plan_dft_2d(nx, ny, 
				reinterpret_cast<fftw_complex*>(output), 
				reinterpret_cast<fftw_complex*>(input), 
				FFTW_BACKWARD, FFTW_ESTIMATE);

		cx = nx*0.5 - 0.5;
		cy = ny*0.5 - 0.5;
		for (ix=0; ix<nx; ix++)
		{
			x = ((double)ix) - cx;
			for (iy=0; iy<ny; iy++)
			{
				y = ((double)iy) - cy;
				index = ix*ny + iy;
				r = sqrt(x*x + y*y) * dx;
				r = r / 8.0; // 16 degree tile
				if (r > 1.0) r = 1.0;
				val = pow(1.0-pow((r-apod_min)/(1.0-apod_min),2.0),2.0);

				if (r > apod_min)
					input[index] = complex<double>(val, 0.0);
				else
					input[index] = complex<double>(1.0,0.0);

			}
		}

		// COMPUTE FT (input -> output)
		fftw_execute(plan);
		
		// compute overlap
		for (ix=0; ix<nx; ix++)
		{
			for (iy=0; iy<ny; iy++)
			{
				index = ix*ny + iy;
				output[index] = output[index]*conj(output[index]);
			}
		}
		// go back to real space to square it again
		fftw_execute(plan2);
		for (ix=0; ix<nx; ix++)
		{
			for (iy=0; iy<ny; iy++)
			{
				index = ix*ny + iy;
				input[index] = input[index]*input[index];
			}
		}
		// once more back to fourier space
		fftw_execute(plan);
		
		delete [] input;

		covar.init(nx,ny,false);

		// load output into covar
		for (ix=0; ix<nx; ix++)
		{
			for (iy=0; iy<ny; iy++)
			{
				index = ix*ny + iy;
				covar.set(ix,iy,output[index]/real(output[0]));
			}
		}

		delete [] output;
		fftw_destroy_plan(plan);
		fftw_cleanup();
	}

	void setPaddedSize ()
	{
		// decide on a size
		offsetx = glob_set.currpadding/dx;
		offsety = offsetx;
		nx = numlons+offsetx*2;
		ny = numlats+offsety*2;

		if (glob_set.myid==0) cout << "Final grid has dimensions "<<nx<<","<<ny<<endl;
	}

	void loadMeasurements ()
	{
		ifstream file;
		long counter, numlines, li, lj, firstmode;
		int ii, ij, ind, indx, indy;
		char buffer[32];
		float fread[6];
		int iread[2];
		bool unique;
		float minlon, maxlon, minlat, maxlat, mindelta;
		float delta_lon, delta_lat;
		// temporary list data
		vector<int> list_k, list_n;
		vector<float> list_vx, list_vy, list_evx, list_evy;
		vector<float> list_lat, list_lon;
		int setsize;
		vector<int> set_ks, set_ns;
		grid<double> vx_local, vy_local;
		MPI_Status stat;
		double err;

		if (glob_set.myid==0)
		{
			file.open(glob_set.input_fname.c_str(), ios::in | ios::binary | ios::ate);
			if (!file.is_open())
			{
				cout << "ERROR: could not open input file.\n";
				return;
			}

			// count number of bytes
			counter = file.tellg();
			if (counter % 32 != 0)
			{
				cout << "ERROR: invalid input file size.\n";
				return;
			}
			numlines = counter / 32;

			cout << "Found " << numlines << " measurements in file.\n";
			
			// error checking
			if (numlines < 10)
			{
				cout << "ERROR: too few measurements in file.\n";
				return;
			}

			// make room in list arrays
			list_n.resize(numlines);
			list_k.resize(numlines);
			list_vx.resize(numlines);
			list_vy.resize(numlines);
			list_evx.resize(numlines);
			list_evy.resize(numlines);
			list_lat.resize(numlines);
			list_lon.resize(numlines);

			// rewind
			file.seekg(0);
			minlon = 1e4;
			maxlon = -1e4;
			minlat = 1e4;
			maxlat = -1e4;
			mindelta = 1e4;
			for (li=0; li<numlines; li++)
			{
				file.read(buffer, 32);
				// copy into temporary buffer
				memcpy(fread+0, buffer+0, sizeof(float));
				memcpy(fread+1, buffer+4, sizeof(float));
				memcpy(iread+0, buffer+8, sizeof(float));
				memcpy(iread+1, buffer+12, sizeof(float));
				memcpy(fread+2, buffer+16, sizeof(float));
				memcpy(fread+3, buffer+20, sizeof(float));
				memcpy(fread+4, buffer+24, sizeof(float));
				memcpy(fread+5, buffer+28, sizeof(float));
				// move into preallocated lists
				list_k[li] = iread[0];
				list_n[li] = iread[1];
				list_vx[li] = fread[2];
				list_vy[li] = fread[4];
				list_evx[li] = fread[3];
				list_evy[li] = fread[5];
				list_lon[li] = fread[0];
				list_lat[li] = fread[1];
				// check for min/max locations
				if (list_lon[li] > maxlon) maxlon = list_lon[li];
				if (list_lon[li] < minlon) minlon = list_lon[li];
				if (list_lat[li] > maxlat) maxlat = list_lat[li];
				if (list_lat[li] < minlat) minlat = list_lat[li];
				if (li > 1)
				{
					delta_lat = fabs(list_lat[li]-list_lat[li-1]);
					delta_lon = fabs(list_lon[li]-list_lon[li-1]);
					if (delta_lon < mindelta && delta_lon > 0.0) mindelta = delta_lon;
					if (delta_lat < mindelta && delta_lat > 0.0) mindelta = delta_lat;
				}
			}
			file.close();
			dx = mindelta;

			// find unique k,n combinations
			// find first allowed mode, we know it's unique
			firstmode = 0;
			if (glob_set.nvals.size() != 0)
				while (find(glob_set.nvals.begin(),glob_set.nvals.end(),list_n[firstmode])==glob_set.nvals.end())
					firstmode ++;
			set.push_back(new mode(list_k[firstmode], list_n[firstmode]));
			for (li=firstmode+1; li<numlines; li++)
			{
				// for speed, check if current meas is any different from prev
				if (list_k[li] != list_k[li-1]
						|| list_n[li] != list_n[li-1])
				{
					// check for mutual uniqueness in mode list
					unique = true;
					for (lj=0; lj<set.size(); lj++)
						if (set[lj]->kval==list_k[li] && set[lj]->nval==list_n[li])
							unique = false;

					// insert if unique
					// and if there is a mode restriction
					if (unique && (glob_set.nvals.size() == 0 || 
							find(glob_set.nvals.begin(),glob_set.nvals.end(),
								list_n[li])!=glob_set.nvals.end()))
						set.push_back(new mode(list_k[li], list_n[li]));
				}
			}
			cout << "Using " << set.size() << " unique modes.\n";
			// NOTE: these are not necessarily modes in common with the kernel set
			// that will be determined later

			// for the lat/lon grid, asume uniform grid within min/max from before
			numlons = (int)((maxlon-minlon)/mindelta + 1);
			numlats = (int)((maxlat-minlat)/mindelta + 1);
			// for informing everyone else
			setsize = set.size();
			for (ii=0; ii<set.size(); ii++)
			{
				set_ks.push_back(set[ii]->kval);
				set_ns.push_back(set[ii]->nval);
			}
		}
		MPI_Bcast(&dx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(&numlons,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&numlats,1,MPI_INT,0,MPI_COMM_WORLD);
		// share the set info
		MPI_Bcast(&setsize,1,MPI_INT,0,MPI_COMM_WORLD);
		set_ks.resize(setsize);
		set_ns.resize(setsize);
		MPI_Bcast(&(set_ks.front()),setsize,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&(set_ns.front()),setsize,MPI_INT,0,MPI_COMM_WORLD);
		// make local copies
		if (glob_set.myid!=0)
			for (ii=0; ii<setsize; ii++)
				set.push_back(new mode(set_ks[ii], set_ns[ii]));
		// everyone has a set with k,n filled out, but no real data
		// tell everyone who will eventually store the real-space data
		for (ii=0; ii<set.size(); ii++)
		{
			set[ii]->proc_owner = ii % glob_set.nproc;
			// go ahead and make space for velocity, root gets all for a bit
			if (set[ii]->proc_owner == glob_set.myid || glob_set.myid == 0)
			{
				set[ii]->vx.init(numlons,numlats,false);
				set[ii]->vy.init(numlons,numlats,false);
			}
		}


		if (glob_set.myid==0)
		{
			// load list measurements into grid
			for (li=0; li<numlines; li++)
			{
				// index in grid
				indx = (list_lon[li]-minlon)/mindelta;
				indy = (list_lat[li]-minlat)/mindelta;
				// find appropriate mode in set
				ind = -1;
				for (ii=0; ii<set.size(); ii++)
					if (list_n[li] == set[ii]->nval && list_k[li] == set[ii]->kval)
						ind = ii;
				if (ind >= 0)
				{
					set[ind]->vx.set(indx, indy, list_vx[li]);
					set[ind]->vy.set(indx, indy, list_vy[li]);
					set[ind]->err += list_evx[li];
				}
			}
			// all grids are loaded, send some away
			MPI_Barrier(MPI_COMM_WORLD);
			for (ii=0; ii<set.size(); ii++)
			{
				set[ii]->err /= (double)(numlons*numlats);
				MPI_Bcast(&(set[ii]->err),1,MPI_DOUBLE,0,MPI_COMM_WORLD);
				if (set[ii]->proc_owner != 0)
				{
					// send it away
					sendGrid(&(set[ii]->vx), set[ii]->proc_owner);
					sendGrid(&(set[ii]->vy), set[ii]->proc_owner);
					// delete local copy
					set[ii]->vx.dealloc();
					set[ii]->vy.dealloc();
				}
			}
		} else {
			// set up to receive velocities
			MPI_Barrier(MPI_COMM_WORLD);
			for (ii=0; ii<set.size(); ii++)
			{
				// which mode are we considering
				MPI_Bcast(&err,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
				set[ii]->err = err;
				if (set[ii]->proc_owner == glob_set.myid)
				{
					// recv from proc0
					recvGrid(&(set[ii]->vx), 0);
					recvGrid(&(set[ii]->vy), 0);
				}
			}
		}

	}





	// given a filename, read in an entire kernel set
	// the kernel set is a binary file with the following format:
	// nkers, nx, ny, nz (ints, 4 bytes each)
	// ker0(x,y,z) (double, 8 bytes * nx * ny * nz)
	// ker1(x,y,z)
	// ...
	// return number of kernels loaded
	// parallel: distribute whole kernels in round-robin sense
	int loadKernelSet (string fname, string zfname)
	{
		ifstream file;
		streamoff pos;
		int newz, count, nkers;
		int ii, ij, ind, thisk, thisn, ix, iy, iz, message;
		double norm, val;
		double *buffer;
		bool found;
		MPI_Status stat;

		if (glob_set.myid==0)
		{
			cout << "Loading kernel set.." << endl;
			file.open(fname.c_str(), ios::in | ios::binary);

			if (!file.is_open())
			{
				cout << "ERROR: could not open kernel file.\n";
				return -1;
			}

			nkers = 0;
			nx_kers = 0;
			ny_kers = 0;
			nz_kers = 0;

			file.read(reinterpret_cast<char*>(&nkers), sizeof(int));
			file.read(reinterpret_cast<char*>(&nx_kers), sizeof(int));
			file.read(reinterpret_cast<char*>(&ny_kers), sizeof(int));
			file.read(reinterpret_cast<char*>(&nz_kers), sizeof(int));
		}

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&nx_kers,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&ny_kers,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&nz_kers,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&nkers,1,MPI_INT,0,MPI_COMM_WORLD);
		depth = new double [nz_kers];
		dz = new double [nz_kers];

		buffer = new double [nx_kers*ny_kers];


		if (glob_set.myid==0)
		{
			// read in kernels one by one
			for (ii=0; ii<nkers; ii++)
			{
				// read in what k,n we are going to read
				file.read(reinterpret_cast<char*>(&thisk), sizeof(int));
				file.read(reinterpret_cast<char*>(&thisn), sizeof(int));
				// look for it in the existing set
				ind = -1;
				for (ij=0; ij<set.size(); ij++)
					if (set[ij]->kval==thisk && set[ij]->nval==thisn)
						ind = ij;
//				cout <<ii<< " k="<<thisk<<" n="<<thisn<<" ind="<<ind<<endl;

				MPI_Bcast(&ind,1,MPI_INT,0,MPI_COMM_WORLD);
				if (ind >= 0)
				{
					set[ind]->has_sker = true;
					if (set[ind]->proc_owner==glob_set.myid)
					{
						// make space for kernel
						set[ind]->allocateRealKernel(nx_kers, ny_kers, nz_kers);
						// load into kernel
						// file is (x,y,z) where x is contiguous
						for (iz=0; iz<nz_kers; iz++)
							for (iy=0; iy<ny_kers; iy++)
								file.read(reinterpret_cast<char*>(set[ind]->sker[iz][iy]), 
										sizeof(double)*nx_kers);
					} else {
						// read and send data
						for (iz=0; iz<nz_kers; iz++)
						{
							for (iy=0; iy<ny_kers; iy++)
								file.read(reinterpret_cast<char*>(&(buffer[iy*nx_kers])), sizeof(double)*nx_kers);
							MPI_Send(buffer,nx_kers*ny_kers,MPI_DOUBLE,set[ind]->proc_owner,1,MPI_COMM_WORLD);
						}
					}
				} else { // CLI arguments tell us not to use this mode
					// skip to next kernel
					pos = file.tellg();
					pos += nx_kers*ny_kers*nz_kers*sizeof(double);
					file.seekg(pos);
				}
			}

			file.close();

		} else {
			for (ii=0; ii<nkers; ii++)
			{
				// hear about a new mode
				MPI_Bcast(&ind,1,MPI_INT,0,MPI_COMM_WORLD);
				// learn who gets this kernel
				if (ind >= 0)
				{
					set[ind]->has_sker = true;
					if (set[ind]->proc_owner == glob_set.myid)
					{
						// set up to receive a kernel
						set[ind]->allocateRealKernel(nx_kers, ny_kers, nz_kers);
						for (iz=0; iz<nz_kers; iz++)
						{
							MPI_Recv(buffer,nx_kers*ny_kers,MPI_DOUBLE,0,1,MPI_COMM_WORLD,&stat);
							for (iy=0; iy<ny_kers; iy++)
								memcpy(set[ind]->sker[iz][iy], &(buffer[iy*nx_kers]), nx_kers*sizeof(double));
						}
					}
				}
			}
		}

		nkers = set.size();

		if (glob_set.myid==0)
		{

			// read in vertical structure
			file.open(zfname.c_str(), ios::in | ios::binary);
			if (!file.is_open())
			{
				cout << "ERROR: could not open z-profile file for kernels.\n";
				return -1;
			}
			file.read(reinterpret_cast<char*>(&newz), sizeof(int));
			if (newz != nz_kers)
			{
				cout << "ERROR: z-profile doesn't have the same number of depths as the kernels.\n";
				return -1;
			}
			file.read(reinterpret_cast<char*>(depth), sizeof(double)*nz_kers);
			file.read(reinterpret_cast<char*>(dz), sizeof(double)*nz_kers);
			file.close();
		}

		// share vertical stuff
		MPI_Bcast(depth,nz_kers,MPI_DOUBLE,0,MPI_COMM_WORLD);
		MPI_Bcast(dz,nz_kers,MPI_DOUBLE,0,MPI_COMM_WORLD);

		// get rid of modes that don't have kernels at this point
		for (ii=0; ii<set.size(); ii++)
		{
			if (!set[ii]->has_sker)
			{
				delete set[ii];
				set.erase(set.begin()+ii);
				ii--;
			}
		}

		// don't trust input kernels, enforce normalization
		MPI_Barrier(MPI_COMM_WORLD);
		if (glob_set.myid==0) cout << "Enforcing Normalization.." << endl;
#pragma omp parallel for private(ii,norm,iz,iy,ix,val) schedule(dynamic,1)
		for (ii=0; ii<set.size(); ii++)
		{
			if (set[ii]->proc_owner == glob_set.myid && set[ii]->has_sker)
			{
				norm = 0.0;
				for (iz=0; iz<nz_kers; iz++)
					for (iy=0; iy<ny_kers; iy++)
						for (ix=0; ix<nx_kers; ix++)
							norm += set[ii]->sker[iz].get(ix,iy) * dz[iz] * dx * dx;
				for (iz=0; iz<nz_kers; iz++)
					for (iy=0; iy<ny_kers; iy++)
						for (ix=0; ix<nx_kers; ix++)
						{
							val = set[ii]->sker[iz].get(ix,iy);
							set[ii]->sker[iz].set(ix,iy,val/norm);
						}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if (glob_set.myid==0) cout << "Done loading kernels.\n";

		
		return nkers;
	}


};
