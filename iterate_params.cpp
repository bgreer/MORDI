
#include "iterate_params.h"

iterate_params::iterate_params ()
{
	num = 0;
	keyparam = -1;
}

void iterate_params::init ()
{
	int ii;

	sl_num[0] = glob_set.depths.size();
	sl_num[1] = glob_set.sigmah_min.size();
	sl_num[2] = glob_set.sigmah.size();
	sl_num[3] = glob_set.sigmaz_min.size();
	sl_num[4] = glob_set.sigmaz.size();
	sl_num[5] = glob_set.reg2.size();
	sl_num[6] = glob_set.reg1.size();
	sl_num[7] = glob_set.reg2z.size();
	sl_num[8] = glob_set.reg1z.size();
	sl_num[9] = glob_set.apodization.size();
	sl_num[10] = glob_set.padding.size();

	num = 1;
	ind = 0;
	keyparam = -1;
	for (ii=0; ii<NUM_PARAMS; ii++)
	{
		sl_ind[ii] = 0;
		num *= sl_num[ii];
	}
	setParams();
}

void iterate_params::reset ()
{
	int ii;
	keyparam = 0;
	ind = 0;
	for (ii=0; ii<NUM_PARAMS; ii++)
		sl_ind[ii] = 0;
	setParams();
}

void iterate_params::setParams()
{
	glob_set.currdepth = glob_set.depths[sl_ind[0]];
	glob_set.currsigmah_min = glob_set.sigmah_min[sl_ind[1]];
	glob_set.currsigmah = glob_set.sigmah[sl_ind[2]];
	glob_set.currsigmaz_min = glob_set.sigmaz_min[sl_ind[3]];
	glob_set.currsigmaz = glob_set.sigmaz[sl_ind[4]];
	glob_set.currreg2 = glob_set.reg2[sl_ind[5]];
	glob_set.currreg1 = glob_set.reg1[sl_ind[6]];
	glob_set.currreg2z = glob_set.reg2z[sl_ind[7]];
	glob_set.currreg1z = glob_set.reg1z[sl_ind[8]];
	glob_set.currapodization = glob_set.apodization[sl_ind[9]];
	glob_set.currpadding = glob_set.padding[sl_ind[10]];
}

bool iterate_params::nextParams()
{
	int ii, ij;
	bool attr = false;

	if (ind >= num-1) return false;

	for (ii=0; ii<NUM_PARAMS && !attr; ii++)
	{
		if (sl_ind[ii] < sl_num[ii]-1)
		{
			keyparam = ii;
			sl_ind[ii] ++;
			attr = true;
			// reset others
			for (ij=0; ij<ii; ij++)
				sl_ind[ij] = 0;
		}
	}
	ind ++;

	setParams();

	switch (keyparam)
	{
		case 10: // padding changed
			if (glob_set.myid==0) 
				cout << "ITERATION: Padding ("<<
					fixed<<glob_set.currpadding<<")\n";
			break;
		case 9: // apodization changed
			if (glob_set.myid==0) 
				cout << "ITERATION: Data Apodization ("<<
					fixed<<glob_set.currapodization<<")\n";
			break;
		case 8:
			if (glob_set.myid==0) 
				cout << "ITERATION: Lambda-z ("<<fixed<<
					glob_set.currreg1z<<")\n"; break;
		case 7:
			if (glob_set.myid==0) 
				cout << "ITERATION: Mu-z ("<<fixed<<
					glob_set.currreg2z<<")\n"; break;
		case 6:
			if (glob_set.myid==0) 
				cout << "ITERATION: Lambda ("<<fixed<<
					glob_set.currreg1<<")\n"; break;
		case 5:
			if (glob_set.myid==0) 
				cout << "ITERATION: Mu ("<<fixed<<
					glob_set.currreg2<<")\n"; break;
		case 4:
			if (glob_set.myid==0) 
				cout << "ITERATION: SigmaZ {" << 
					glob_set.currsigmaz<<")\n"; break;
		case 3:
			if (glob_set.myid==0) 
				cout << "ITERATION: SigmaZ-min {" << 
					glob_set.currsigmaz_min<<")\n"; break;
		case 2:
			if (glob_set.myid==0) 
				cout << "ITERATION: SigmaH {" << 
					glob_set.currsigmah<<")\n"; break;
		case 1:
			if (glob_set.myid==0) cout << 
				"ITERATION: SigmaH-min {" << 
					glob_set.currsigmah_min<<")\n"; break;
		case 0:
			if (glob_set.myid==0) 
				cout << "ITERATION: Depth ("<<fixed<<
					glob_set.currdepth<<")\n"; break;
	}

	return true;
}


