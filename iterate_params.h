#ifndef ITERATE_H
#define ITERATE_H
#include "settings.h"
#define NUM_PARAMS 11

// The inversion code loops over many different variables to do a 
// parameter space search. This bit of code handles that looping
// so that it doesn't clutter up main.cpp

extern settings glob_set;

class iterate_params
{
public:
	int num, ind, keyparam;
	// sub-loop stuff:
	int sl_ind[NUM_PARAMS];
	int sl_num[NUM_PARAMS];

	iterate_params();

	void init();
	void reset();
	void setParams();
	bool nextParams();
	
};

#endif
