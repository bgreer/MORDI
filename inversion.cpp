#include "inversion.h"


// CONSTRUCTOR //
inversion::inversion ()
{
	has_init = false;
	alloc_vel = false;
}

void inversion::init (int x, int y, int z, int num)
{
	int ii;
	if (has_init && (x != nx || y != ny || num != nn))
	{
		// deallocate first
		cout << "ERROR: reallocating inversion not yet implemented!\n";
		exit(EXIT_FAILURE);
	}
	nx = x;
	ny = y;
	nz = z;
	nn = num;
	// space for the coefficients
	coefs = new grid<complex<double> > [nn];
	for (ii=0; ii<nn; ii++)
		coefs[ii].init(nx,ny,true);
	target_ft.init(nx,ny,true);
	target_z = new double [nz];
	has_init = true;
}

void inversion::loadTarget ()
	{
		int ix, iy, iz;
		complex<double> *input, *output;
		grid<complex<double> > local;
		double *buffer;
		fftw_plan plan;
		stringstream fname;
		ifstream file;

		if (glob_set.myid==0)
		{
			input = new complex<double> [nx*ny];
			output = new complex<double> [nx*ny];
			plan = fftw_plan_dft_2d(nx, ny,
					reinterpret_cast<fftw_complex*>(input),
					reinterpret_cast<fftw_complex*>(output),
					FFTW_FORWARD, FFTW_ESTIMATE);

			cout << "Loading Target from file.." << endl;
			// make filename based on just depth
			fname << glob_set.target_fname << "_" << fixed << setprecision(2) << glob_set.currdepth;
			file.open(fname.str().c_str(), ios::in | ios::binary);
			file.read(reinterpret_cast<char*>(&ix), sizeof(nx));
			file.read(reinterpret_cast<char*>(&iy), sizeof(ny));
			file.read(reinterpret_cast<char*>(&iz), sizeof(nz));
			if (ix!=nx || iy!=ny || iz!=nz)
			{
				cout << "ERROR: Target Dimensions don't match inversion dimensions!" << endl;
				cout << "nx = " << nx << ", target nx = " << ix << endl;
				cout << "ny = " << ny << ", target ny = " << iy << endl;
				cout << "nz = " << nz << ", target nz = " << iz << endl;
				exit(-1);
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		buffer = new double [nx];
		target_3d = new grid<complex<double> > [nz];
		for (iz=0; iz<nz; iz++)
			target_3d[iz].init(nx,ny,true);

		if (glob_set.myid==0)
		{
			local.init(nx,ny,false);
			buffer = new double [nx];
			for (iz=0; iz<nz; iz++)
			{
				// read slice from file
				for (iy=0; iy<ny; iy++)
				{
					file.read(reinterpret_cast<char*>(buffer), nx*sizeof(buffer[0]));
					for (ix=0; ix<nx; ix++)
						input[iy*nx + ix] = buffer[ix];
				}
				// ft
				fftw_execute(plan);
				// unload to target
				for (iy=0; iy<ny; iy++)
					for (ix=0; ix<nx; ix++)
						local.set(ix,iy,output[iy*nx+ix]);
				// share
				distributeGrid(&local, &(target_3d[iz]), 0);
			}
			delete [] buffer;
			file.close();
			fftw_destroy_plan(plan);
		} else {
			for (iz=0; iz<nz; iz++)
				distributeGrid(&local, &(target_3d[iz]), 0);
		}

	}


void inversion::setTarget (double widthx, double widthz, measurement_set *mset)
	{
		int ix, iy, iz, ind;
		complex<double> *input, *output;
		grid<complex<double> > local;
		fftw_plan plan;
		double x, y, tot;

		// depth
		for (iz=0; iz<nz; iz++)
		{
			target_z[iz] = exp(-pow((mset->depth[iz] - glob_set.currdepth)/widthz,2.0));
		}

		// horizontal
		input = new complex<double> [nx*ny];
		output = new complex<double> [nx*ny];
		plan = fftw_plan_dft_2d(nx, ny, 
				reinterpret_cast<fftw_complex*>(input), 
				reinterpret_cast<fftw_complex*>(output), 
				FFTW_FORWARD, FFTW_ESTIMATE);

		// load
		tot = 0.0;
		for (ix=0; ix<nx; ix++)
		{
			x = (ix - nx*0.5 + 0.5)*mset->dx;
			for (iy=0; iy<ny; iy++)
			{
				y = (iy - ny*0.5 + 0.5)*mset->dx;
				ind = ix*ny + iy;
				input[ind] = exp(-pow(x/widthx,2.0) - pow(y/widthx,2.0));
				for (iz=0; iz<nz; iz++)
					tot += real(input[ind]) * target_z[iz] * mset->dz[iz] * mset->dx * mset->dx;
			}
		}
		// execute
		fftw_execute(plan);
		// unload
		for (ix=0; ix<nx; ix++)
		{
			for (iy=0; iy<ny; iy++)
			{
				ind = ix*ny + iy;
				target_ft.set(ix,iy,output[ind]/tot);
			}
		}

		fftw_destroy_plan(plan);
	}

void inversion::loadCoefs (measurement_set *mset)
	{
		int ii, ix, iy;
		int nx_file, ny_file, nn_file;
		int *nvals, *kvals;
		ifstream file;
		grid<complex<double> > local;

		// load header
		if (glob_set.myid==0)
		{
			cout << "Loading coefficients.." << endl;
			file.open(glob_set.coefs_load_fname.c_str(), ios::binary);
			file.read(reinterpret_cast<char*>(&nx_file), sizeof(nx_file));
			file.read(reinterpret_cast<char*>(&ny_file), sizeof(ny_file));
			file.read(reinterpret_cast<char*>(&nn_file), sizeof(nn_file));
			// error checking
			if (nn_file != nn)
				cout << "WARNING: coef file has " << nn_file << " modes, current run has " << nn << endl;
			if (nx_file != nx || ny_file != ny)
			{
				cout << "ERROR: coef file dims do not match current dims!" << endl;
				cout << "	nx_coefs = " << nx_file << "   nx = " << nx << endl;
				cout << "   ny_coefs = " << ny_file << "   ny = " << ny << endl;
				exit(-1);
			}
			// load n/k vals
			nvals = new int [nn_file];
			kvals = new int [nn_file];
			for (ii=0; ii<nn_file; ii++)
			{
				file.read(reinterpret_cast<char*>(&(nvals[ii])), sizeof(int));
				file.read(reinterpret_cast<char*>(&(kvals[ii])), sizeof(int));
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);


		if (glob_set.myid==0)
		{
			file.close();
			delete [] nvals;
			delete [] kvals;
		}
	}

void inversion::saveCoefs (measurement_set *mset)
	{
		int ii, ix, iy;
		double lambda, mu, sh, sz;
		ofstream file;
		stringstream fname;
		grid<complex<double> > local;

		if (glob_set.myid==0)
		{
			cout << "Saving Coefficients.." << endl;
			local.init(nx,ny,false); // local grid
			lambda = glob_set.currreg1 + glob_set.currreg1z * glob_set.currdepth;
			mu = glob_set.currreg2 + glob_set.currreg2z * glob_set.currdepth;
			sz = max(glob_set.currsigmaz_min, glob_set.currsigmaz * glob_set.currdepth);
			sh = max(glob_set.currsigmah_min, glob_set.currsigmah * glob_set.currdepth);
			fname << glob_set.coefs_save_fname
				 << "_" << fixed << setprecision(2) << glob_set.currdepth
				 << "_" << setprecision(3) << lambda
				 << "_" << setprecision(3) << mu
				 << "_" << setprecision(1) << glob_set.currpadding
				 << "_" << setprecision(1) << glob_set.currapodization
				 << "_" << setprecision(3) << sz
				 << "_" << setprecision(3) << sh;
			file.open(fname.str().c_str(), ios::binary);
			// write a basic header
			// nx ny nn
			file.write((char*)(&nx), sizeof(int));
			file.write((char*)(&ny), sizeof(int));
			file.write((char*)(&nn), sizeof(int));
			// write a list of n,k pairs
			for (ii=0; ii<nn; ii++)
			{
				file.write((char*)(&(mset->set[ii]->nval)), sizeof(int));
				file.write((char*)(&(mset->set[ii]->kval)), sizeof(int));
			}
		}

		for (ii=0; ii<nn; ii++)
		{
			cout << ii << cout;
			// everyone wait
			MPI_Barrier(MPI_COMM_WORLD);
			// send grid for mode ii to root
			collectGrid(&(coefs[ii]), &local, 0);
			if (glob_set.myid==0)
				for (iy=0; iy<ny; iy++)
					file.write((char*)(local.data[iy]), nx*sizeof(complex<double>));
		}

		if (glob_set.myid==0)
			file.close();

	}

void inversion::solveCoefs (measurement_set *mset)
	{
		int x0, x1, y0, y1;
		
		// each proc has their own block of wavenumbers
		// use a kernel_ft to figure out how much that is
		x0 = mset->set[0]->sker_ft[0].x0;
		x1 = mset->set[0]->sker_ft[0].x1;
		y0 = mset->set[0]->sker_ft[0].y0;
		y1 = mset->set[0]->sker_ft[0].y1;

		// manual parallelization
#pragma omp parallel
	{
		int ix, iy, iz, ii, ij;
		int modei, modej, start, end, myid, num;
		complex<double> *rhs, *ktk, *work;
		complex<double> zero(0.0,0.0), i(0.0,1.0);
		complex<double> val, target_val;
		double lambda, mu;
		int info, *piv, worksize;

		ktk = new complex<double> [nn*nn];
		rhs = new complex<double> [nn];
		piv = new int [nn];
		worksize = nn*32;
		work = new complex<double> [worksize];

		lambda = pow(10.0, glob_set.currreg1 + glob_set.currreg1z * glob_set.currdepth);
		mu = pow(10.0, glob_set.currreg2 + glob_set.currreg2z * glob_set.currdepth);

		// split up tasks
		num = y1-y0+1;
		myid = omp_get_thread_num();
		start = myid*num / glob_set.nthreads;
		end = (myid+1)*num / glob_set.nthreads - 1;
		start += y0;
		end += y0;

		// loop through each wavenumber
		for (iy=start; iy<=end; iy++)
		{
			for (ix=x0; ix<=x1; ix++)
			{
				// compute ktk + r
				memset(ktk,0x00,nn*nn*sizeof(complex<double>));
				for (ii=0; ii<nn; ii++)
				{
					for (ij=0; ij<nn; ij++)
					{
						ktk[ii*nn+ij] += mset->overlap[ii][ij].get(ix,iy);
						// mode-mode covariance
					/*	
						if (kset->set[keri]->n==kset->set[kerj]->n &&
								abs(kset->set[keri]->k - kset->set[kerj]->k) == 1)
							ktk[ii*nn+ij] += meas->err[modei] * meas->err[modej] * 
								(mu*(1.0+i) + lambda*kset->covar.get(ix,iy)) * 
								(1.6*pow(kset->set[keri]->k/46.,0.25)/(kset->set[keri]->n+2.0));
					*/	
					}
					// regularization terms
					ktk[ii*nn+ii] += mset->set[ii]->err * mset->set[ii]->err * 
						(mu*(1.0) + lambda*mset->covar.get(ix,iy));
				}
						
				// LU decomposition
				zgetrf(&nn, &nn, reinterpret_cast<MKL_Complex16*>(ktk), &nn, piv, &info);
				// Inverse
				zgetri(&nn, reinterpret_cast<MKL_Complex16*>(ktk), &nn, piv, 
						reinterpret_cast<MKL_Complex16*>(work), &worksize, &info);

				if (info!=0)
				{
					cout << "Proc"<<glob_set.myid << " inv status "<<info<<endl;
					for (ii=0; ii<nn; ii++)
						for (ij=0; ij<nn; ij++)
							cout << ii << " " << ij << " " << ktk[ii*nn+ij]<<endl;
					exit(EXIT_FAILURE);
				}
				// compute RHS
				for (ii=0; ii<nn; ii++)
				{
					rhs[ii] = zero;
					for (iz=0; iz<nz; iz++)
					{
						target_val = target_z[iz] * target_ft.get(ix,iy) * mset->dx * mset->dx * mset->dz[iz];
						if (glob_set.target_fname!="")
							target_val = target_3d[iz].get(ix,iy) * mset->dx * mset->dx * mset->dz[iz];
						rhs[ii] += mset->set[ii]->sker_ft[iz].get(ix,iy) * target_val;
					}
				}
				// compute coefs
				for (ii=0; ii<nn; ii++)
				{
					val = zero;
					for (ij=0; ij<nn; ij++)
						val += ktk[ii*nn+ij] * rhs[ij];
					coefs[ii].set(ix,iy,val);
				}
			}
		}
	}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	// to get error, loop through each mode
	//   multiply coefs_ft * covariance
	//   ft into real space
	//   multiply by real space coefs
	//   sum over all space
void inversion::computeError (measurement_set *mset)
	{
		int ii, ix, iy, x0, x1, y0, y1;
		grid<complex<double> > local, recv, recv2;
		complex<double> val, zero(0.0,0.0);
		complex<double> *input, *output, *input2, *output2;
		fftw_plan plan, plan2;
		double normalize;

		err = zero;

		local.init(nx,ny,true);
		if (glob_set.myid==0)
		{
			coefs_real.init(nx,ny,false);
			recv.init(nx,ny,false);
			recv2.init(nx,ny,false);
			input = new complex<double> [nx*ny];
			output = new complex<double> [nx*ny];
			input2 = new complex<double> [nx*ny];
			output2 = new complex<double> [nx*ny];
			plan = fftw_plan_dft_2d(nx, ny, 
					reinterpret_cast<fftw_complex*>(input), 
					reinterpret_cast<fftw_complex*>(output), 
					FFTW_BACKWARD, FFTW_ESTIMATE);
			plan2 = fftw_plan_dft_2d(nx, ny, 
					reinterpret_cast<fftw_complex*>(input2), 
					reinterpret_cast<fftw_complex*>(output2), 
					FFTW_BACKWARD, FFTW_ESTIMATE);
			for (iy=0; iy<ny; iy++)
				for (ix=0; ix<nx; ix++)
					coefs_real.set(ix,iy,zero);
			normalize = 0.0;
			for (ii=0; ii<nn; ii++)
				normalize += real(coefs[ii].get(0,0))*mset->dx;
		}

		// get local bounds
		x0 = mset->set[0]->sker_ft[0].x0;
		x1 = mset->set[0]->sker_ft[0].x1;
		y0 = mset->set[0]->sker_ft[0].y0;
		y1 = mset->set[0]->sker_ft[0].y1;

		for (ii=0; ii<nn; ii++)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			// bring raw coefs together
			collectGrid(&(coefs[ii]),&recv,0);

			for (iy=y0; iy<=y1; iy++)
			{
				for (ix=x0; ix<=x1; ix++)
				{
					val = coefs[ii].get(ix,iy) * mset->covar.get(ix,iy);
					local.set(ix,iy,val);
				}
			}

			collectGrid(&local, &recv2, 0);

			if (glob_set.myid==0)
			{
				// root proc fts both
				for (iy=0; iy<ny; iy++)
				{
					for (ix=0; ix<nx; ix++)
					{
						input[iy*nx+ix] = recv.get(ix,iy)/normalize;
						input2[iy*nx+ix] = recv2.get(ix,iy)/normalize;
					}
				}
				fftw_execute(plan);
				fftw_execute(plan2);
				for (iy=0; iy<ny; iy++)
				{
					for (ix=0; ix<nx; ix++)
					{
						// store coefs for output
						val = coefs_real.get(ix,iy) + output[iy*nx+ix];
						coefs_real.set(ix,iy,val);
	
						err += output[iy*nx+ix] * output[iy*nx+ix]
							* mset->set[ii]->err * mset->set[ii]->err / ((double)(nx*ny));
					}
				}
			}
		}
		
		MPI_Barrier(MPI_COMM_WORLD);
	}

	// the velocity field is the sum over modes of the coefs convolved with the measurements
	// so compute the multiplication of the fts in parallel,
	// adding each to a cumulative parallel grid
	// then condense to proc0 for output
void inversion::computeVelocity (measurement_set *mset)
	{
		int ii, ind, ix, iy, x0, x1, y0, y1;
		grid<complex<double> > localx, localy;
		grid<complex<double> > recvx, recvy;
		complex<double> val, zero(0.0,0.0), *input, *output;
		double normalize;
		fftw_plan plan;

		localx.init(nx,ny,true);
		localy.init(nx,ny,true);

		if (glob_set.myid==0)
		{
			recvx.init(nx,ny,false);
			recvy.init(nx,ny,false);
			if (!alloc_vel)
			{
				vx.init(nx,ny,false);
				vy.init(nx,ny,false);
				alloc_vel = true;
			}
			normalize = 0.0;
			for (ii=0; ii<nn; ii++)
				normalize += real(coefs[ii].get(0,0))*mset->dx;
		}
		// get local bounds
		x0 = mset->set[0]->sker_ft[0].x0;
		x1 = mset->set[0]->sker_ft[0].x1;
		y0 = mset->set[0]->sker_ft[0].y0;
		y1 = mset->set[0]->sker_ft[0].y1;
#pragma omp parallel for private(ii,ind,ix,iy,val) schedule(dynamic,1)
		for (ii=0; ii<nn; ii++)
		{
			// multiply coefs with measurements
			for (ix=x0; ix<=x1; ix++)
			{
				for (iy=y0; iy<=y1; iy++)
				{
					if (ii==0)
					{
						localx.set(ix,iy,zero);
						localy.set(ix,iy,zero);
					}
					// TODO: separate coefs for vy?
					val = localx.get(ix,iy);
					localx.set(ix,iy,val + (mset->set[ii]->vx_ft.get(ix,iy) * coefs[ii].get(ix,iy)));
					val = localy.get(ix,iy);
					localy.set(ix,iy,val + (mset->set[ii]->vy_ft.get(ix,iy) * coefs[ii].get(ix,iy)));
				}
			}
		}

		// collect grids
		collectGrid(&localx,&recvx,0);
		collectGrid(&localy,&recvy,0);

		// take real part
		if (glob_set.myid==0)
		{
			input = new complex<double> [nx*ny];
			output = new complex<double> [nx*ny];
			plan = fftw_plan_dft_2d(nx, ny, 
					reinterpret_cast<fftw_complex*>(input), 
					reinterpret_cast<fftw_complex*>(output), 
					FFTW_BACKWARD, FFTW_ESTIMATE);
			
			// load data
			for (ix=0; ix<nx; ix++)
			{
				for (iy=0; iy<ny; iy++)
				{
					ind = ix*ny+iy;
					input[ind] = recvx.get(ix,iy)/normalize;
				}
			}
			// execute
			fftw_execute(plan);
			// unload
			for (ix=0; ix<nx; ix++)
			{
				for (iy=0; iy<ny; iy++)
				{
					ind = ix*ny + iy;
					vx.set(ix,iy,real(output[ind])/(nx*ny));
				}
			}
			// load data
			for (ix=0; ix<nx; ix++)
			{
				for (iy=0; iy<ny; iy++)
				{
					ind = ix*ny+iy;
					input[ind] = recvy.get(ix,iy)/normalize;
				}
			}
			// execute
			fftw_execute(plan);
			// unload
			for (ix=0; ix<nx; ix++)
			{
				for (iy=0; iy<ny; iy++)
				{
					ind = ix*ny + iy;
					vy.set(ix,iy,real(output[ind])/(nx*ny));
				}
			}
			fftw_destroy_plan(plan);
		}
	}

void inversion::clearVelocity ()
	{
		if (glob_set.myid==0 && alloc_vel)
		{
			vx.dealloc();
			vy.dealloc();
			alloc_vel = false;
		}
	}

void inversion::outputVelocity ()
	{
		int ix, iy, ind;
		stringstream fname;
		double lambda, mu, sz, sh;
		ofstream file;
		float temp;
		char *buffer;
	
		if (glob_set.myid==0)
		{
			lambda = glob_set.currreg1 + glob_set.currreg1z * glob_set.currdepth;
			mu = glob_set.currreg2 + glob_set.currreg2z * glob_set.currdepth;
			sz = max(glob_set.currsigmaz_min, glob_set.currsigmaz * glob_set.currdepth);
			sh = max(glob_set.currsigmah_min, glob_set.currsigmah * glob_set.currdepth);

			// construct fname
			fname << glob_set.output_fname
				 << "_" << fixed << setprecision(2) << glob_set.currdepth
				 << "_" << setprecision(3) << lambda
				 << "_" << setprecision(3) << mu
				 << "_" << setprecision(1) << glob_set.currpadding
				 << "_" << setprecision(1) << glob_set.currapodization
				 << "_" << setprecision(3) << sz
				 << "_" << setprecision(3) << sh;

			cout << "Writing result to file " << fname.str() << endl;

			// write everything as floats?
			buffer = new char [ny * 8 * sizeof(float)];
			
			file.open(fname.str().c_str(), ios::binary);
			// write header as ints
			memcpy(buffer+0, &nx, sizeof(int));
			memcpy(buffer+sizeof(int), &ny, sizeof(int));
			file.write(buffer, 2*sizeof(int));
			
			for (ix=0; ix<nx; ix++)
			{
				memset(buffer, 0x00, ny*8*sizeof(float));
				for (iy=0; iy<ny; iy++)
				{
					temp = ix;
					memcpy(buffer+(iy*8 + 0)*sizeof(float), &temp, sizeof(float));
					temp = iy;
					memcpy(buffer+(iy*8 + 1)*sizeof(float), &temp, sizeof(float));
					temp = glob_set.currdepth;
					memcpy(buffer+(iy*8 + 2)*sizeof(float), &temp, sizeof(float));
					temp = vx.get(ix,iy);
					memcpy(buffer+(iy*8 + 3)*sizeof(float), &temp, sizeof(float));
					temp = sqrt(real(err));
					memcpy(buffer+(iy*8 + 4)*sizeof(float), &temp, sizeof(float));
					temp = vy.get(ix,iy);
					memcpy(buffer+(iy*8 + 5)*sizeof(float), &temp, sizeof(float));
					temp = 0.0;
					memcpy(buffer+(iy*8 + 6)*sizeof(float), &temp, sizeof(float));
					temp = real(coefs_real.get(ix,iy));
					memcpy(buffer+(iy*8 + 7)*sizeof(float), &temp, sizeof(float));
					/*
					file << ix << " " << iy << " " << setprecision(3) << fixed
						<< glob_set.currdepth << " "  << scientific
						<< vx.get(ix,iy) << " " << sqrt(real(err)) << " "
						<< vy.get(ix,iy) << " 0.0 "
						<< real(coefs_real.get(ix,iy)) << endl;
					*/
				}
				file.write(buffer, ny*8*sizeof(float));
			}
			delete [] buffer;
			file.close();
		}
	}

void inversion::outputAvgker (measurement_set *mset)
	{
		int ii, ix, iy, iz;
		grid<complex<double> > avgker_local, avgker;
		complex<double> zero(0.0,0.0);
		complex<double> val, *input, *output;
		double mu, lambda, sz, sh, *buffer;
		ofstream file;
		stringstream fname;
		fftw_plan plan;

		// everyone has some set of wavenumbers for the coefs and the kers
		// march through each depth
		// for each mode, multiply coefs and kers, collect on proc0, ft and output

		avgker_local.init(nx,ny,true);
		if (glob_set.myid==0)
		{
			avgker.init(nx,ny,false);
			buffer = new double [nx];
			// construct fname
			lambda = glob_set.currreg1 + glob_set.currreg1z * glob_set.currdepth;
			mu = glob_set.currreg2 + glob_set.currreg2z * glob_set.currdepth;
			sz = max(glob_set.currsigmaz_min, glob_set.currsigmaz * glob_set.currdepth);
			sh = max(glob_set.currsigmah_min, glob_set.currsigmah * glob_set.currdepth);
			fname << glob_set.avgker_fname
				 << "_" << fixed << setprecision(2) << glob_set.currdepth
				 << "_" << setprecision(3) << lambda
				 << "_" << setprecision(3) << mu
				 << "_" << setprecision(1) << glob_set.currpadding
				 << "_" << setprecision(1) << glob_set.currapodization
				 << "_" << setprecision(3) << sz
				 << "_" << setprecision(3) << sh;

			cout << "Writing averaging kernel to file " << fname.str() << endl;
			input = new complex<double> [nx*ny];
			output = new complex<double> [nx*ny];
			plan = fftw_plan_dft_2d(nx, ny, 
					reinterpret_cast<fftw_complex*>(input), 
					reinterpret_cast<fftw_complex*>(output), 
					FFTW_BACKWARD, FFTW_ESTIMATE);

			file.open(fname.str().c_str(), ios::out | ios::binary);
			// output dimensions as 4-byte integers
			file.write(reinterpret_cast<char*>(&nx), sizeof(nx));
			file.write(reinterpret_cast<char*>(&ny), sizeof(ny));
			file.write(reinterpret_cast<char*>(&nz), sizeof(nz));
		}
		MPI_Barrier(MPI_COMM_WORLD);

		for (iz=0; iz<nz; iz++)
		{
			// sum up each mode
			for (ii=0; ii<nn; ii++)
			{
				for (iy=0; iy<ny; iy++)
				{
					for (ix=0; ix<nx; ix++)
					{
						if (ii==0) avgker_local.set(ix,iy,zero);
						val = avgker_local.get(ix,iy);
						val += (coefs[ii].get(ix,iy) * mset->set[ii]->sker_ft[iz].get(ix,iy));
						avgker_local.set(ix,iy,val);

					}
				}
			}
			// share this result
			collectGrid(&avgker_local, &avgker, 0);
			if (glob_set.myid==0)
			{
				// load data
				for (iy=0; iy<ny; iy++)
					for (ix=0; ix<nx; ix++)
						input[iy*nx+ix] = avgker.get(ix,iy);
				// execute
				fftw_execute(plan);
				// unload
				for (iy=0; iy<ny; iy++)
				{
					for (ix=0; ix<nx; ix++)
						buffer[ix] = real(output[iy*nx+ix]);
					file.write(reinterpret_cast<char*>(buffer), nx*sizeof(buffer[0]));
				}
			}
		}

		if (glob_set.myid==0)
		{
			file.close();
			fftw_destroy_plan(plan);
		}

	}
