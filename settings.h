#ifndef SETTINGS_H
#define SETTINGS_H
#include <algorithm>
#include <sstream>
#include "tclap/CmdLine.h"

using namespace std;

class settings
{
public:
	string kernel_set_fname, kernel_z_fname, input_fname, output_fname;
	string avgker_fname, target_fname, coefs_save_fname, coefs_load_fname;
	vector<double> depths;
	vector<int> nvals;
	vector<double> padding; // in degrees
	vector<double> apodization; // in degrees
	vector<double> reg1, reg1z, reg2, reg2z; // regularization terms
	vector<double> sigmaz, sigmaz_min, sigmah, sigmah_min;

	double currpadding, currapodization, currdepth;
	double currreg1, currreg1z, currreg2, currreg2z;
	double currsigmaz, currsigmaz_min, currsigmah, currsigmah_min;
	bool apod_circ;

	int verbose;
	// parallelization
	int nthreads, myid, nproc;

	settings ()
	{
		// do nothing yet?
		myid = 0;
		nproc = 1;
		verbose = 9001;
		apod_circ = false;
	}

	void setParallel (int me, int num)
	{
		myid = me;
		nproc = num;
	}

	void parse (int argc, char *argv[])
	{
		// temporary variables needed for init
		string depthlist, nlist, token, plist, alist;
		string sreg1, sreg1z, sreg2, sreg2z, ssz, sszm, ssh, sshm;
		double tempd;
		int tempi;

		try{
		TCLAP::CmdLine cmd("Multi-Channel Deconvolution Code", ' ', "0.9");


		// kernel set stuff
		TCLAP::ValueArg<std::string> kset("k","kernel_set",
				"Kernel set to load (binary file)",false,"input/kers_binary","filename");
		cmd.add(kset);
		TCLAP::ValueArg<std::string> ksetz("z","kernel_set_depth",
				"Kernel set depths and dz (binary file)",false,"input/kers_z","filename");
		cmd.add(ksetz);
		// input data
		TCLAP::ValueArg<std::string> inputdata("i","input",
				"Input fits data (binary file)",true,"","filename");
		cmd.add(inputdata);
		// depths
		TCLAP::ValueArg<string> depths("d","depths",
				"Set of depths to do, in megameters",true,"0.0","comma separated list");
		cmd.add(depths);
		// restrictions on processing
		TCLAP::ValueArg<string> nrange("n","n-values",
				"Set of radial orders to use",false,"","comma separated list");
		cmd.add(nrange);
		// output filename base
		TCLAP::ValueArg<std::string> outputdata("o","output",
				"Output filename",true,"","filename");
		cmd.add(outputdata);
		// regularization
		TCLAP::ValueArg<string> lambda("l","lambda",
				"log10 regularization applied to the full error term",false,
				"-1.0","comma separated list");
		cmd.add(lambda);
		TCLAP::ValueArg<string> lambda2("L","lambda-z",
				"how the regularization on the full error term changes in depth, in log10 per Mm depth",
				false,"-1.0","comma separated list");
		cmd.add(lambda2);
		TCLAP::ValueArg<string> mu("m","mu",
				"log10 regularization applied to the truncated error term",false,
				"-1.0","comma separated list");
		cmd.add(mu);
		TCLAP::ValueArg<string> mu2("M","mu-z",
				"how the regularization on the truncated error term changes in depth, in log10 per Mm depth",
				false,"-1.0","comma separated list");
		cmd.add(mu2);
		// averaging kernel
		TCLAP::ValueArg<std::string> avgker("a","avgker",
				"Averaging kernel filename",false,"","filename");
		cmd.add(avgker);
		// threading
		TCLAP::ValueArg<int> threads("t","threads",
				"Number of OpenMP threads to use",false,1,"number");
		cmd.add(threads);
		// tuning (apodize data, padding)
		TCLAP::ValueArg<string> apod("A","apod",
				"Amount of apodization to apply to the measurements (degrees)",false,"5.0","comma separated list");
		cmd.add(apod);
		TCLAP::SwitchArg circapod("C","circular","Circular apodization of data instead of rectangular", cmd, false);
		TCLAP::ValueArg<string> pad("p","padding",
				"Amount of padding to use around domain (degrees)",false,"16.0","comma separated list");
		cmd.add(pad);
		// averaging kernel
		TCLAP::ValueArg<std::string> target_file("T","target",
				"Target filename to override depth",false,"","filename");
		cmd.add(target_file);
		// target widths
		TCLAP::ValueArg<string> sz("s","sigmaz",
				"Width in depth of target per Mm depth (degrees)",false,"0.0","comma separated list");
		cmd.add(sz);
		TCLAP::ValueArg<string> szm("S","sigmaz_min",
				"Minimum width in depth of target (degrees)",false,"0.1","comma separated list");
		cmd.add(szm);
		TCLAP::ValueArg<string> sh("g","sigmah",
				"Horizontal width of target per Mm depth (degrees)",false,"0.0","comma separated list");
		cmd.add(sh);
		TCLAP::ValueArg<string> shm("G","sigmah_min",
				"Minimum horizontal width of target (degrees)",false,"0.1","comma separated list");
		cmd.add(shm);

		// a-coef saving / loading
		TCLAP::ValueArg<std::string> coef_save_file("b","coef_save",
				"Filename for saving inversion coefficients",false,"","filename");
		cmd.add(coef_save_file);
		TCLAP::ValueArg<std::string> coef_load_file("B","coef_load",
				"Filename for loading inversion coefficients (no inversion will be done)",false,"","filename");
		cmd.add(coef_load_file);


		// parse things!
		cmd.parse( argc, argv );

		// get results
		kernel_set_fname = kset.getValue();
		kernel_z_fname = ksetz.getValue();
		depthlist = depths.getValue();
		nlist = nrange.getValue();
		input_fname = inputdata.getValue();
		output_fname = outputdata.getValue();
		avgker_fname = avgker.getValue();
		target_fname = target_file.getValue();
		coefs_save_fname = coef_save_file.getValue();
		coefs_load_fname = coef_load_file.getValue();
		nthreads = threads.getValue();
		alist = apod.getValue();
		plist = pad.getValue();
		sreg1 = lambda.getValue();
		sreg1z = lambda2.getValue();
		sreg2 = mu.getValue();
		sreg2z = mu2.getValue();
		apod_circ = circapod.getValue();
		ssz = sz.getValue();
		sszm = szm.getValue();
		ssh = sh.getValue();
		sshm = shm.getValue();

		if (coefs_save_fname != "" && coefs_load_fname != "")
		{
			cout << "ERROR: doesn't make sense to both load and safe coefs. Try again." << endl;
			exit(-1);
		}

		} catch (TCLAP::ArgException &e)  // catch any exceptions
		{ std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl; }


		// do extra parsing on the lists
		istringstream ss(depthlist);
		while (getline(ss,token,','))
			depths.push_back(atof(token.c_str()));
		istringstream ss2(nlist);
		while (getline(ss2,token,','))
			nvals.push_back(atoi(token.c_str()));
		istringstream ss3(plist);
		while (getline(ss3,token,','))
			padding.push_back(atof(token.c_str()));
		istringstream ss4(alist);
		while (getline(ss4,token,','))
			apodization.push_back(atof(token.c_str()));

		istringstream ss5(sreg1);
		while (getline(ss5,token,','))
			reg1.push_back(atof(token.c_str()));
		istringstream ss6(sreg1z);
		while (getline(ss6,token,','))
			reg1z.push_back(atof(token.c_str()));
		istringstream ss7(sreg2);
		while (getline(ss7,token,','))
			reg2.push_back(atof(token.c_str()));
		istringstream ss8(sreg2z);
		while (getline(ss8,token,','))
			reg2z.push_back(atof(token.c_str()));
		istringstream ss9(ssz);
		while (getline(ss9,token,','))
			sigmaz.push_back(atof(token.c_str()));
		istringstream ss10(sszm);
		while (getline(ss10,token,','))
			sigmaz_min.push_back(atof(token.c_str()));
		istringstream ss11(ssh);
		while (getline(ss11,token,','))
			sigmah.push_back(atof(token.c_str()));
		istringstream ss12(sshm);
		while (getline(ss12,token,','))
			sigmah_min.push_back(atof(token.c_str()));
	}

	// print current settings to stdout
	void printDetails ()
	{
		int ii;
		if (myid == 0)
		{
		cout << endl;
		cout << "--- 3D Multi-Channel Deconvolution ---" << endl;
		cout << "    Verbose Level: " << verbose << endl;
		cout << "    " << endl;
		if (verbose > 0) cout << " -Run Details-" << endl;
		if (verbose > 1) cout << "    Kernel Set: " << kernel_set_fname << endl;
		if (verbose > 2) cout << "    Kernel Depths: " << kernel_z_fname << endl;
		if (verbose > 0) cout << "    Input Data: " << input_fname << endl;
		if (verbose > 0) cout << "    Output Name: " << output_fname << endl;
		if (verbose > 0 && avgker_fname=="") cout << "    No Averaging Kernel"<< endl;
		if (verbose > 0 && avgker_fname!="") cout << "    Averaging Kernel Name: "<<avgker_fname<< endl;
		if (verbose > 0) cout << "    Number of Depths: " << depths.size() << endl;
		if (verbose > 0 && nvals.size()>0) cout << "    Radial Orders: " << nvals.size() << endl;
		if (verbose > 0 && nvals.size()==0) cout << "    No Radial Order Restriction" << endl;
		if (verbose > 0) cout << "    Amount of Apodization: ";
		setprecision(2);
		for (ii=0; ii<apodization.size(); ii++) {cout << fixed << apodization[ii]; if (ii<apodization.size()-1) cout << ", ";}
		cout << endl;
		if (verbose > 0) cout << "    Amount of Padding: ";
		for (ii=0; ii<padding.size(); ii++) {cout << fixed << padding[ii]; if (ii<padding.size()-1) cout << ", ";}
		cout << endl;
		if (verbose > 0) cout << "    Number of Processes: " << nproc << endl;
		if (verbose > 0) cout << "    Number of Threads: " << nthreads << endl;
		if (verbose > 0 && coefs_save_fname != "") cout << "    Saving coefs to: " << coefs_save_fname << endl;
		if (verbose > 0 && coefs_load_fname != "") cout << "    Loading coefs from: " << coefs_load_fname << endl;
		if (verbose > 0) cout << endl;
		}
	}
};
#endif
