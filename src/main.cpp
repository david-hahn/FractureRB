#include <fstream>
#include <tclap/CmdLine.h>

#include "BulletWrapper.h"

using namespace FractureSim;
using namespace std;

int main(int argc, char* argv[]){
	setbuf(stdout, NULL);
	try {
		printf("\n%% parsing command line: \n%% %s ...\n%%  ",argv[0]);
		for(int i=1; i<argc; ++i) printf(" %s",argv[i]);
		
		TCLAP::CmdLine cmd("Rigid Body and Fracture simulation", ' ', "using HyENA 2.0, Bullet 2.83 and OpenVDB 2.2 backend");
		TCLAP::UnlabeledValueArg<string> bulletFile( "input-file", "Bullet format input scene", true, "default", "input-file", cmd );
		TCLAP::UnlabeledValueArg<string> paramFile(  "param-file", "parameter file (CSV format)", true, "default", "param-file", cmd );
        TCLAP::ValueArg<string> outDir(              "o","out","output directory and file prefix", false, "","out-dir", cmd);
		TCLAP::ValueArg<double> rbTimestep (         "d","time-step" ,"rigid body time step size", false, 1.0/250.0,"time-step", cmd);
		TCLAP::ValueArg<int>    rbTimesteps(         "n","time-steps","number of rigid body time steps", false, 1000,"time-steps", cmd);
		TCLAP::ValueArg<double> forceThreshold(      "f","min-force","minimum total deformational force per rigid body timestep that triggers a fracture simulation", false, DBL_MAX,"impulse", cmd);
		TCLAP::ValueArg<double> impulseThreshold(    "i","min-impulse","minimum total impulse per rigid body timestep that triggers a fracture simulation", false, DBL_MAX,"impulse", cmd);
		TCLAP::ValueArg<double> impulseSplit    (    "s","split-impulse","depth in units of the minimal voxel size of all input objects where impulse splitting occurs in the RB solver (default disabled)", false, -1.0,"depth", cmd);
		TCLAP::SwitchArg        useDefaultRBsolver(  "" ,"default-solver","use Bullet's default impulse solver instead of Danzig MLCP",cmd);
		TCLAP::ValueArg<int>    estFromElems(        "e","est-elems","if a fracture simulation contains more than the specified number of elements switch to faster estimated SIF evaluation (default -1=never) ", false, -1,"est-elems", cmd);
		
		cmd.parse( argc, argv );

		BulletWrapper rbsim(useDefaultRBsolver.getValue());
		rbsim.setOutputDir(outDir.getValue()); //ToDo: make this a constructor param. of BulletWrapper (?)
		printf("\n%% reading scene %s with parameters %s\n", bulletFile.getValue().c_str(), paramFile.getValue().c_str());
		rbsim.initSceneFromFile(
			bulletFile.getValue(),
			paramFile.getValue(),
			impulseSplit.getValue(),
			estFromElems.getValue()
		);
		printf("\n%% scene import done\n");

		string outfile = outDir.getValue()+"sim.mel";
//		string outfile = outDir.getValue()+"sim.py"; //ToDo? add cli-param to choose between MEL and BPY output
		ofstream outSim(outfile.c_str());

		printf("\n%% timestepping the rigid body sim ...\n");
		for (int i = 0; i < rbTimesteps.getValue(); i++) {
			rbsim.stepSimulation(
				&outSim, BulletWrapper::OUT_MEL,
				rbTimestep.getValue(),
				impulseThreshold.getValue(),
				forceThreshold.getValue());
		}

		outSim.close();

		//# we can do this in Maya's python interface (probably similar for Blender?):
		//# which is more stable and less memory consuming than the "Source script" command
		//import os
		//import maya.mel
		//from pymel.core import sceneName
		//os.chdir(os.path.abspath(os.path.join(sceneName(), os.pardir)))
		//f  = open("_out/t2_sim.mel")
		//for line in f:
		// maya.mel.eval(line)
		//f.close()
		
		printf("\n%% simulation done!\n");
	}catch (TCLAP::ArgException &e){
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return -1;
	}
	return 0;
}
