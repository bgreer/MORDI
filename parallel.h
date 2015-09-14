#ifndef PARALLEL_H
#define PARALLEL_H

extern settings glob_set;

template<typename T> MPI_Datatype convertType(T thing)
{
	float fl;
	double dbl;
	int in;
	long lon;
	complex<double> cmplxdbl;

	if (typeid(thing) == typeid(fl)) return MPI_FLOAT;
	if (typeid(thing) == typeid(dbl)) return MPI_DOUBLE;
	if (typeid(thing) == typeid(in)) return MPI_INT;
	if (typeid(thing) == typeid(lon)) return MPI_LONG;
	if (typeid(thing) == typeid(cmplxdbl)) return MPI::DOUBLE_COMPLEX;
	return MPI_BYTE;
}

template<typename Type>
void distributeGrid(grid<Type> *src, grid<Type> *dest, int root)
{
	Type a;
	int iy, recipient;
	MPI_Status stat;

	MPI_Barrier(MPI_COMM_WORLD);

	// root has the entire src grid
	// everyone has a dest grid
	// root needs to send 
	
	if (root == glob_set.myid)
	{
		// check dimensions
		if (src->nx != dest->nx || src->ny != dest->ny 
				|| src->parallel == true || dest->parallel == false
				|| src->alloc == false || dest->alloc == false)
		{
			cout << "ERROR: grid distribution failed.\n";
			return;
		}

		// send data to everyone, copy if local
		for (iy=0; iy<src->ny; iy++)
		{
			// if local
			if (iy >= dest->y0 && iy <= dest->y1)
				memcpy(dest->data[iy-dest->y0], src->data[iy], src->nx * sizeof(Type));
			else
			{
				recipient = ((iy+1)*glob_set.nproc-1)/src->ny;
				MPI_Send(src->data[iy], src->nx, convertType(a), recipient, 0, MPI_COMM_WORLD);
			}
		}
	} else {
		for (iy=0; iy<dest->ny_local; iy++)
		{
			MPI_Recv(dest->data[iy], dest->nx, convertType(a), root, 0, MPI_COMM_WORLD, &stat);
		}
	}
}

// everyone has the src grid, root has the dest grid
template<typename Type>
void collectGrid(grid<Type> *src, grid<Type> *dest, int root)
{
	Type a;
	int iy, sender;
	MPI_Status stat;

	MPI_Barrier(MPI_COMM_WORLD);

	if (root == glob_set.myid)
	{
		// check dimensions
		if (src->nx != dest->nx || src->ny != dest->ny 
				|| src->parallel == false || dest->parallel == true
				|| src->alloc == false || dest->alloc == false)
		{
			cout << "ERROR: grid collection failed.\n";
			return;
		}

		for (iy=0; iy<src->ny; iy++)
		{
			// if local
			if (iy >= src->y0 && iy <= src->y1)
				memcpy(dest->data[iy-src->y0], src->data[iy], src->nx * sizeof(Type));
			else
			{
				sender = ((iy+1)*glob_set.nproc-1)/src->ny;
				MPI_Recv(dest->data[iy], src->nx, convertType(a), sender, 0, MPI_COMM_WORLD, &stat);
			}
		}
	} else {
		for (iy=0; iy<src->ny_local; iy++)
			MPI_Send(src->data[iy], src->nx, convertType(a), root, 0, MPI_COMM_WORLD);
	}
}

template<typename Type>
void sendGrid(grid<Type> *src, int proc)
{
	Type a;
	int iy;

	for (iy=0; iy<src->ny; iy++)
		MPI_Send(src->data[iy], src->nx, convertType(a), proc, 0, MPI_COMM_WORLD);
	
}

template<typename Type>
void recvGrid(grid<Type> *dest, int proc)
{
	Type a;
	int iy;
	MPI_Status stat;

	for (iy=0; iy<dest->ny; iy++)
		MPI_Recv(dest->data[iy], dest->nx, convertType(a), proc, 0, MPI_COMM_WORLD, &stat);
}

#endif
