#include "Reader.h"
#include "FractureBEM.h"
#include "MaterialModel.h"
#include "SubsampledCrackTip.h"

#include <tclap/CmdLine.h>
#include <omp.h>

#include <sstream>

using namespace FractureSim;
using namespace std;

int main(int argc, char* argv[]){
	double t0=omp_get_wtime(), t1;
	double startTime=t0;
	MaterialModel* material=NULL;
	try {
		int ret;
		printf("\n%% parsing command line ... ");
		TCLAP::CmdLine cmd("Fracture simulation with symmetric Galerkin BEM", ' ', "using HyENA 2.0 and OpenVDB 2.2 backend");
		// simulation params (input file, boundary conditions, timesteps, output, stuff...)
		TCLAP::UnlabeledValueArg<string> inpFile("input-file", 
			"Elmer format input files containing the BEM mesh", true, "default", "input-file",
			cmd
		);
		TCLAP::UnlabeledMultiArg<string> bcStrings("boundary-condition",
			"specifies a boundary condition on mesh-region i either i(fixed) for hom. Dirichlet, i(crack) for cracks or i(f_x,f_y,f_z) for tractions",
			true, "boundary-condition", cmd
		);
		TCLAP::MultiArg<unsigned int> crackTipIDs( "c", "crack-tip",
			"Boundary IDs of crack tips (must be line elements in the boundary file).", false,
			"crack-tip", cmd
		);
		TCLAP::MultiArg<unsigned int> staticCracks( "", "static-crack",
			"Region-IDs of cracks that are not allowed to propagate.", false,
			"static-crack", cmd
		);
		// option to start cracks at user specified locations before the first timestep
		TCLAP::MultiArg<string>   startCracks("X", "start-crack",
			"Start a crack at (p) with face normal (n) and tip normal (N)", false,
			"(px,py,pz,nx,ny,nz,Nx,Ny,Nz)", cmd
		);
		TCLAP::SwitchArg is2D("", "two-d", "pseudo-2d computation, locks all displacements along z-axis" /*, cmd*/); // disabled
		TCLAP::SwitchArg noSI("", "no-si", "disable crack self intersection testing", cmd);
		TCLAP::ValueArg<unsigned int> crackSteps("s","steps","number of crack propagation steps", false, 0,"steps", cmd);
		TCLAP::ValueArg<unsigned int>  cpVersion("v","cp-version","select crack propagation algorithm (0..original, 1..improved, 2..experimental) default is 0", false, 0,"cp-version", cmd);
		//output params
        TCLAP::ValueArg<string> outFile("o","out","name of VTK/VDB output files", false, "","out-file", cmd);
		TCLAP::SwitchArg outVDBsubstep("", "out-sub", "write VDB output for every substep (as opposed to full steps only)", cmd);
		TCLAP::ValueArg<double> outMeshQual("","vis-qual","quality parameter for hi-res output, in [0,1] use VDB adaptive mesh, < 0 decimate to quadric error tolerance, > 1 no hi-res output (default)", false, 2.0,"vis-qual", cmd);
		TCLAP::SwitchArg outMeshNoU("", "vis-no-disp", "do NOT apply any displacements when producing hi-res output (use with --vis-qual)", cmd);
		TCLAP::SwitchArg outMeshNoCOD("", "vis-no-cod", "do NOT apply crack opening displacements when producing hi-res output (use with --vis-qual)", cmd);
		TCLAP::SwitchArg outMeshClose("", "vis-close", "make crack faces coincident for 0 COD in hi-res output (experimental) (use with --vis-qual)", cmd);
		TCLAP::SwitchArg outMeshOBJ("", "vis-obj", "write hi-res output as .OBJ file (otherwise as .STL) (use with --vis-qual)", cmd);
		// resolution params for BEM mesh and level-set
		TCLAP::ValueArg<double> resMesh("r","res-mesh","resolution of BEM mesh, default is 0.1", false, 0.1,"res-mesh", cmd);
		TCLAP::ValueArg<double> resGrid("R","res-grid","resolution of VDB grid, default is 0.01", false, 0.01,"res-grid", cmd);
		TCLAP::ValueArg<double> cfSampl("", "cf-sampl","sampling density of high-res crack-fronts relative to grid size, default is 0.8", false, 0.8,"density", cmd);
		TCLAP::ValueArg<int> remeshTris("T","remesh","remesh input from level-set decimating to the given number of triangles; uses boundary regions definition from file 'input-file.region'", false, -1,"num-tris", cmd);
		TCLAP::ValueArg<int>    nCracks("x","max-cracks","seeding cut-off for new cracks, set to 0 for no seeding and -1 (default) for unlimited seeding", false, -1,"max-cracks", cmd);
		TCLAP::ValueArg<double> offsetMesh("O","offset","offset coarse mesh by n-voxels along the vertex normals after remeshing (use with -T, default is 0)", false, 0.0,"n-voxels", cmd);
		// material params (elasticity constants, fracture toughnes = max. stress intensity, tensile strength = max. principal stress)
		TCLAP::ValueArg<double>     youngsMod("E","youngs-modulus","elastic modulus of material in Pa, default is 1.0", false, 1.0,"youngs-modulus", cmd);
		TCLAP::ValueArg<double> poissonsRatio("n","poissons-ratio","poisson's ratio of material, default is 0.2", false, 0.2,"poissons-ratio", cmd);
		TCLAP::ValueArg<double>       density("d","density","density of material, default is 1.0", false, 1.0,"density", cmd);
		TCLAP::ValueArg<double>     toughness("K","toughness","fracture toughness, ie. max. stress intensity at crack-tips", false, 1.0,"toughness", cmd);
		TCLAP::ValueArg<double>      strength("S","strength", "tensile strength, ie. max. principal stress at surfaces", false, 1.0,"strength", cmd);
		TCLAP::ValueArg<double>      compress("C","compressive", "factor by which the material is stronger and tougher under compression", false, 3.0,"compressive", cmd);
		TCLAP::ValueArg<string>        matMdl("M","material","select a material model, default \"h\" (homogeneous), use \"?\" to print a list of options", false, "h","material", cmd);
		// other stuff
		TCLAP::ValueArg<int> method(        "m","method","choose 0 .. full BEM solution (default) or 1 .. initial BEM solution with estimated SIFs", false, FULL_BEM,"method", cmd);
		TCLAP::ValueArg<int> estFromElems(  "e","est-elems","once the simulation contains more than the specified number of elements switch to faster estimated SIF evaluation (default -1=never) ", false, -1,"est-elems", cmd);

		cmd.parse( argc, argv );

		printf("\n%% ... input files are \"%s.*\"", inpFile.getValue().c_str());
		printf("\n%% ... material parameters (E,nu) are (%.2le, %.2lf)", youngsMod.getValue(), poissonsRatio.getValue());
		if(is2D.getValue()) printf("\n%% ... computing only 2d displacements in x-y-plane");
		if(outFile.isSet()) printf("\n%% ... output files are \"%s_*.vtk/.vdb\"", outFile.getValue().c_str());
		t1=omp_get_wtime();
		printf("\n%% ... %.4lfs CLI parsed",t1-t0);
		t0=t1;

		printf("\n%% initializing simulation %s ... ", (method.getValue()==INITIAL_BEM)?"(using estimated SIFs from initial BEM solution)":"");
		FractureBEM sim( resMesh.getValue(), (FRACTURE_METHOD) (method.getValue()) ); //resolution of crack-meshes (crack-tips will be subdivided, also "time-step" of propagation is chosen to match this)
		material = createMaterialModel(  matMdl.getValue(),
			youngsMod.getValue(), poissonsRatio.getValue(), density.getValue(),
			strength.getValue(),  toughness.getValue(),     compress.getValue()
		);
		if( material==NULL ) return -1;
		sim.setCPVersion(cpVersion.getValue());

		printf("\n%% reading model ...      ");
		{
			BEMReader reader(inpFile.getValue());
			// reading model
			set<unsigned int> region_ids/*leave empty to read all regions*/, crack_tip_ids/*read these line-elements as crack tips*/;
			crack_tip_ids.insert(crackTipIDs.getValue().begin(), crackTipIDs.getValue().end()); //fill from CLI
			ret=reader.readModel(
				sim.getElems(), sim.getRegions(), sim.getCrackTips(), sim.getCrackTipParents(), sim.getNodes(),
				FractureSim::TRI, region_ids, FractureSim::LINE, crack_tip_ids
			);
			t1=omp_get_wtime();
			printf("\t%.4lfs - have %d elements and %d nodes",t1-t0, sim.getElems().size() ,sim.getNodes().size());
			t0=t1;
		}

		if(remeshTris.isSet()){ // do remeshing based on the level-set representation of the initial model to build the BEM model
			RegionHandler regionHandler;
            // if using remeshing, load the definitions and apply them to the input mesh
            // as we remesh from a level-set, there can NOT be any initial cracks!
			printf("\n%% initializing mesh regions from file %s.regions ...",inpFile.getValue().c_str());
			sim.getRegions().clear();
			
			regionHandler.readRegionDefinitions(inpFile.getValue()+".regions");
			regionHandler.assignRegions(sim.getRegions(), sim.getNodes(),sim.getElems());

            printf("\n%% initializing VDB grids  ...");
			sim.initVDB(resGrid.getValue(), noSI.getValue() );
            t1=omp_get_wtime();
            printf("\t%.4lfs",t1-t0);
            t0=t1;

            if(remeshTris.getValue()>3){ // remeshing to the desired resolution
                printf("\n%% remeshing            ... ");
				sim.remesh(remeshTris.getValue(), 0.1, offsetMesh.getValue() );
                regionHandler.assignRegions(sim.getRegions(), sim.getNodes(),sim.getElems());
                t1=omp_get_wtime();
                printf("\t%.4lfs - have %d elements",t1-t0,sim.getElems().size());
                t0=t1;
            }

            printf("\n%% initializing BEM solver ... ");
            ret=sim.initBEM(*material, is2D.getValue(), bcStrings.getValue());
			if(ret!=0){ printf("\n%% initBEM returned %d!",ret); return ret;}
		}else{ // work directly with the input mesh, supporting initial cracks
            // since initial cracks are set by boundary conditions
            // initialize the BEM first - parsing BCs in the process
            printf("\n%% initializing BEM solver ... ");
            ret=sim.initBEM(*material, is2D.getValue(), bcStrings.getValue());
			if(ret!=0){ printf("\n%% initBEM returned %d!",ret); return ret;}
            t1=omp_get_wtime();
            printf("\t%.4lfs",t1-t0);
            t0=t1;
			// if we have initial cracks, the BEM must be initialized first, because this will parse the boundary conditions and mark these regions as cracks
			// once this has been done, the level-set representation can be built accordingly
            printf("\n%% initializing VDB grids  ...");
            sim.initVDB(resGrid.getValue(), noSI.getValue() );
        }
		t1=omp_get_wtime();
		printf("\t%.4lfs",t1-t0);
		t0=t1;	

		// allow the user to pre-arrest some cracks (for testing)
        {
			id_set staticCrackSet;
			staticCrackSet.insert(staticCracks.getValue().begin(),staticCracks.getValue().end());
			sim.setCracksInactive(staticCrackSet);
		}
        
		printf("\n%% initializing fracture model ... ");
		sim.initFractureModel(*material, outVDBsubstep.isSet(),cfSampl.getValue());

		// can add cracks manually here ...
		//sim.startCrack( Eigen::Vector3d(0.0,0.8,0.2), Eigen::Vector3d(0.0,1.0,0.0), Eigen::Vector3d( 1.0,0.0,0.0) );
		//sim.startCrack( Eigen::Vector3d(0.0,0.4,0.1), Eigen::Vector3d(0.0,1.0,0.0), Eigen::Vector3d(-1.0,0.0,0.0) );
		//sim.computeAddedCracks();
		if( startCracks.isSet() ){
			bool needUpdate=false;
			for(vector<string>::const_iterator it = startCracks.getValue().begin();
				it != startCracks.getValue().end(); ++it){
				Eigen::Vector3d p,n,m;
				char c,d;
				d=sscanf(it->c_str(), "(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf%c",
					&(p[0]),&(p[1]),&(p[2]), &(n[0]),&(n[1]),&(n[2]), &(m[0]),&(m[1]),&(m[2]), &c
				);
				if( c==')' && d==10 ){
					printf("\n%% adding crack \"%s\"", it->c_str());
					sim.startCrack(p, n.normalized(), m.normalized());
					needUpdate=true;
				}
			}
			if( needUpdate ) sim.computeAddedCracks();
		}
		sim.getHiResCrackTip().setSIFsOutputFile("sifs");

		// *********** the simulation loop ***************
		int newCracks=0, steps = crackSteps.getValue(); // number of crack propagation steps
		bool doFracture=true; int k;
		for(k=0; k<=steps && doFracture; ++k){
			if( (sim.getMethod()==FULL_BEM) && estFromElems.getValue() > 0 )
				if( sim.getElems().size() > estFromElems.getValue() ) sim.setMethodEstimated();

			printf("\n%% computing SIFs           ...");
			sim.computeBEM();
			t1=omp_get_wtime();
			printf("\t%.4lfs",t1-t0);
			t0=t1;

			//if(k==0)
			//	sim.dumpDebugData();

			if( outFile.isSet() ){
				printf("\n%% writing VTK & VDB output ...");
				stringstream fileName;
				fileName << outFile.getValue() << "_";
				fileName.fill('0'); fileName.width(4);
				fileName << k;
				sim.writeMesh(fileName.str(),
                    outMeshQual.getValue(), !outMeshNoU.isSet(), !outMeshNoCOD.isSet(),
                    outMeshClose.isSet(), outMeshOBJ.isSet()
                );
                sim.writeVDB(fileName.str());
                sim.writeCrackTip(fileName.str()+"_ct");

				t1=omp_get_wtime();
				printf("\t%.4lfs (%s)",t1-t0, fileName.str().c_str());
				t0=t1;
			}

			doFracture = true;
			unsigned int numElems= sim.getElems().size();
			if(k<steps && doFracture){
				printf("\n%% simulating fractures ... ");

				newCracks=0; // just for output
				if( nCracks.getValue() < 0 ) // unlimited seeding
					newCracks = sim.seedCracksAndPropagate();
				else if ( sim.getCracks().size() < nCracks.getValue() ) // limited seeding
					newCracks = sim.seedCracksAndPropagate( nCracks.getValue() - sim.getCracks().size());
				else // no seeding
					sim.propagateCracks();
				
				if(sim.getElems().size() == numElems) doFracture=false; // stop if no new elements have been added (in this the BEM solution will not change, so the next step won't do anything either)
				numElems=sim.getElems().size();

				t1=omp_get_wtime();
				printf("\t%.4lfs - have %d elements",t1-t0,sim.getElems().size());
				if(newCracks>0) printf(", just added %d new fractures", newCracks);
				t0=t1;
			}
		}
        if( outFile.isSet() ){ // write final state
            sim.computeBEM(); // make sure we have a consistent state to start with
            printf("\n%% writing final VTK & VDB output ... ");
            stringstream fileName;
            fileName << outFile.getValue() << "_";
            fileName.fill('0'); fileName.width(4);
            fileName << k;
            sim.writeMesh(fileName.str(),
                outMeshQual.getValue(), !outMeshNoU.isSet(), !outMeshNoCOD.isSet(),
                outMeshClose.isSet(), outMeshOBJ.isSet()
            );
            sim.writeVDB(fileName.str());
            sim.writeCrackTip(fileName.str()+"_ct");

            t1=omp_get_wtime();
            printf("\t%.4lfs (%s)",t1-t0, fileName.str().c_str());
            t0=t1;
        }

        t1=omp_get_wtime();
		printf("\n\n%% all done \n%% total time %.2lf seconds\n", t1-startTime);
        
		if(material!=NULL){ delete material; material==NULL;}
	}catch (TCLAP::ArgException &e)  {// catch any exceptions
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		if(material!=NULL){ delete material; material==NULL;}
		return -1;
	}
	return 0;
}
