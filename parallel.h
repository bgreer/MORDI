#ifndef PARALLEL_H
#define PARALLEL_H
#include <typeinfo>
#include "mpi.h"
#include "grid.h"
//extern settings glob_set;


template<typename T> MPI_Datatype convertType(T thing);

template<typename Type>
void distributeGrid(grid<Type> *src, grid<Type> *dest, int root);

template<typename Type>
void collectGrid(grid<Type> *src, grid<Type> *dest, int root);

template<typename Type>
void sendGrid(grid<Type> *src, int proc);

template<typename Type>
void recvGrid(grid<Type> *dest, int proc);

#include "parallel.tcc"

#endif
