#ifndef ITERATE_H
#define ITERATE_H
#include "settings.h"
#define NUM_PARAMS 11

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
