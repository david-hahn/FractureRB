#include "types.h"
#include "Reader.h"
#include "VDBLoader.h"

#include <tclap/CmdLine.h>
#include <sstream>

using namespace FractureSim;
using namespace std;

int readMeshVTK(
    string filename, node_map& nodes, elem_map& elems, id_map& regions,
    vector_type& u, vector_type& u_c
);

int removeCompressiveCOD(
	node_map& nodes, elem_map& elems, id_map& regions,
    id_set& cracks, vector_type& u
);

int main(int argc, char* argv[]){
	try {
		printf("%% parsing command line ...\n");
		TCLAP::CmdLine cmd("Post process FractureBEM results to OBJ or STL visual quality meshes", ' ', "using OpenVDB 2.2 backend");
		// input params
		TCLAP::UnlabeledValueArg<string> inFile("input-file", 
			"VTK & VDB input files containing the simulation results", true,
            "default", "name", cmd
		);
		TCLAP::ValueArg<int> fullSteps("s","steps","number of crack propagation steps", false, 0,"steps", cmd);
		TCLAP::ValueArg<int> subSteps ("S","substeps","number of crack propagation substeps", false, -1,"substeps", cmd);
		TCLAP::ValueArg<int> startStep("x","start","start at the given full step (skip previous ones, default 0)", false, 0,"substeps", cmd);
		TCLAP::ValueArg<double> resMesh("r","res-mesh","resolution of BEM mesh as in the simulation, used for nearest triangle search, default is 0.1", false, 0.1,"res-mesh", cmd);

		//output params
        TCLAP::ValueArg<string> outFile("o","out","name of VTK/VDB output files", true, "","out-file", cmd);
		TCLAP::ValueArg<double> outMeshVDBQual("q","vdb-qual","quality parameter for VDB meshing in [0,1] (default is 0.01)", false, 0.01,"vdb-qual", cmd);
		TCLAP::ValueArg<double> outMeshVCGQual("Q","vis-qual","quality parameter for mesh decimation: quadric error tolerance (default is no decimation)", false, -1.0,"max-error", cmd);
		TCLAP::ValueArg<double> outMeshPostQual("P","post-qual","quality parameter for post-interpolation mesh decimation: quadric error tolerance (default is no decimation)", false, -1.0,"max-error", cmd);
		TCLAP::ValueArg<double> outMeshSegment("","segment","if thr>0.0 segment before meshing with handle-removal, default: no segmentation, use thr<=1.0", false, -1.0,"thr", cmd);
		TCLAP::ValueArg<double> outLastSegment("","segment-last","same as segment, but applies only to last (full step) output, default: no segmentation", false, -1.0,"thr", cmd);
		TCLAP::SwitchArg outMeshNoU("", "vis-no-disp", "do NOT apply any displacements when producing hi-res output", cmd);
		TCLAP::SwitchArg outMeshNoCOD("", "vis-no-cod", "do NOT apply crack opening displacements when producing hi-res output", cmd);
		TCLAP::SwitchArg outMeshClose("", "vis-close", "make crack faces coincident for 0 COD in hi-res output", cmd);
		TCLAP::SwitchArg outMeshAllowNegCOD("", "vis-allow-neg-cod", "allow compressive crack opening displacement (may lead to mesh fold-over), otherwise compressive COD will be projected out of the visual mesh", cmd);
		TCLAP::SwitchArg outMeshOBJ("", "vis-obj", "write hi-res output as .OBJ file (otherwise as .STL)", cmd);
		TCLAP::SwitchArg outMeshSegments("", "segment-files", "write one file per segment instead of all segments in a single file", cmd);
		cmd.parse( argc, argv );

        printf("%% reading input files %s_####.vtk/.vdb as full steps", inFile.getValue().c_str());
        if(subSteps.isSet())
            printf(" and %s_####_####.vdb as sub-steps",inFile.getValue().c_str());
        printf("\n%% expecting full steps numbered 0 to %u", fullSteps.getValue());
        if(subSteps.isSet())
            printf(" with sub-steps 0 to %u each", subSteps.getValue());
        printf("\n%% output files will be %s_#####.%s\n%% ...\n",
            outFile.getValue().c_str(), outMeshOBJ.isSet()?"obj":"stl" );
        
        // load coarse mesh (nodes, elems, regions, displacements, crack-displacements)
        // -- for two full steps
        // load vdb grids for one full step
        // load vdb grids for substep-update and apply in each substep
        // interpolate coarse displacements in time on the later mesh
        // produce hi-res mesh output
        
        unsigned int count = startStep.getValue() + startStep.getValue()*(1+subSteps.getValue());
		double seg = outMeshSegment.getValue();
        int ret;
        stringstream file, subfile, outfile;
        
        node_map nodes, nodes2;
        elem_map elems, elems2;
        id_map   regions, regions2;
        vector_type u, u_c, u2, u_c2;
        id_set   cracks;
		/*_quad_** // for quadratic interpolation, we need a third full step
        node_map nodes3;
        elem_map elems3;
        id_map   regions3;
        vector_type u3, u_c3;
		/**/
        
        VDBLoader vdbLoader;
        vdbLoader.addToNearTriGrid(nodes,elems,resMesh.getValue()); // empty init
        
        for(int fullStep=startStep.getValue(); fullStep<=fullSteps.getValue(); ++fullStep){
            // load/copy VTK
            if(fullStep==startStep.getValue()){ // load first step
				/*_quad_** // also preload the second step
                file.str(""); file.clear();
                file << inFile.getValue() << "_";
                file.fill('0'); file.width(4);
                file << (fullStep+1);
				printf("\n%% reading %s.vtk",file.str().c_str());
                ret =readMeshVTK( file.str()+".vtk", nodes2, elems2, regions2, u2, u_c2);
				/**/
                // build file-name
                file.str(""); file.clear();
                file << inFile.getValue() << "_";
                file.fill('0'); file.width(4);
                file << fullStep;
				printf("\n%% reading %s.vtk",file.str().c_str());
                ret =readMeshVTK( file.str()+".vtk", nodes, elems, regions, u, u_c);
            }else{ // copy from already loaded following step
                nodes=nodes2; elems=elems2; regions=regions2; u=u2; u_c=u_c2;
				/*_quad_**
				nodes2=nodes3; elems2=elems3; regions2=regions3; u2=u3; u_c2=u_c3;
				/**/
            }
			if(fullStep==fullSteps.getValue() && outLastSegment.isSet()) seg = outLastSegment.getValue();
            // load VDB
			/*_quad_** // need to rebuild file name, since we mess it up when preloading the 3rd vtk
	        file.str(""); file.clear();
            file << inFile.getValue() << "_";
            file.fill('0'); file.width(4);
            file << fullStep; /**/
			printf("\n%% reading %s.vdb",file.str().c_str());
            ret+=vdbLoader.loadFull(file.str()+".vdb" ,cracks);
            if( ret<0){
                printf("\n%% Error reading full step %s\n", (file.str()+".vtk/.vdb").c_str());
                return -1;
            }
            vdbLoader.addToNearTriGrid(nodes,elems);
			if( !outMeshAllowNegCOD.isSet() && fullStep==startStep.getValue()){
                removeCompressiveCOD(nodes, elems, regions, cracks, u);
                printf("\n%% removing compressive COD from %s.vtk", file.str().c_str());
            }

			
            // write output for the full step
            printf("\n%% %04d      --> %05d ... ", fullStep,count);
            outfile.str(""); outfile.clear();
            outfile << outFile.getValue() << "_";
            outfile.fill('0'); outfile.width(5);
            outfile << (count++); // post-inc!

            vdbLoader.writeVisualMesh(
                nodes,elems,regions,cracks,u,u_c,
                outfile.str(), outMeshVDBQual.getValue(), outMeshVCGQual.getValue(),
                !outMeshNoU.isSet(),!outMeshNoCOD.isSet(),
				outMeshClose.isSet(),outMeshOBJ.isSet(),
				seg, outMeshSegments.isSet(), outMeshPostQual.getValue()
            );
			printf(" ok");

			/*_quad_**
			vdbLoader.addToNearTriGrid(nodes2,elems2);
            if(fullStep<fullSteps.getValue()-1){ // except for the second to last, preload the next full step
                // build file-name
                file.str(""); file.clear();
                file << inFile.getValue() << "_";
                file.fill('0'); file.width(4);
                file << (fullStep+2);
				printf("\n%% reading %s.vtk",file.str().c_str());
                ret =readMeshVTK( file.str()+".vtk", nodes3, elems3, regions3, u3, u_c3);
                if(ret<0) return ret; 
                // prepare for interpolation: pad u,u_c with zeros to size of u3,u_c3
                int oldSize=u.cols(), newSize=u3.cols();
                u.conservativeResizeLike(u3);
                u_c.conservativeResizeLike(u_c3);
                u.tail(newSize-oldSize).setZero();
                u_c.tail(newSize-oldSize).setZero();
				oldSize=u2.cols(), newSize=u3.cols();
                u2.conservativeResizeLike(u3);
                u_c2.conservativeResizeLike(u_c3);
                u2.tail(newSize-oldSize).setZero();
                u_c2.tail(newSize-oldSize).setZero();
            } /*/
            if(fullStep<fullSteps.getValue()){ // except for the last, preload the next full step
                // build file-name
                file.str(""); file.clear();
                file << inFile.getValue() << "_";
                file.fill('0'); file.width(4);
                file << (fullStep+1);
				printf("\n%% reading %s.vtk",file.str().c_str());
                ret =readMeshVTK( file.str()+".vtk", nodes2, elems2, regions2, u2, u_c2);
                if(ret>=0) vdbLoader.addToNearTriGrid(nodes2,elems2);
                // prepare for interpolation: pad u,u_c with zeros to size of u2,u_c2
                int oldSize=u.cols(), newSize=u2.cols();
                u.conservativeResizeLike(u2);
                u_c.conservativeResizeLike(u_c2);
                u.tail(newSize-oldSize).setZero();
                u_c.tail(newSize-oldSize).setZero();
            } /**/
            for(int subStep=0; subStep<=subSteps.getValue() &&
                ret>=0 && fullStep<fullSteps.getValue();
                ++subStep){// do substeps, except for last full step
                subfile.str(""); subfile.clear();
                subfile << inFile.getValue() << "_";
                subfile.fill('0'); subfile.width(4);
                subfile << fullStep << "_";
                subfile.fill('0'); subfile.width(4);
                subfile << subStep;
				//printf("\n%% reading %s.vdb",subfile.str().c_str());
                vdbLoader.loadUpdate(subfile.str()+".vdb",cracks);

                if( !outMeshAllowNegCOD.isSet() && subStep==0 && fullStep<fullSteps.getValue()){
                    removeCompressiveCOD(nodes, elems, regions, cracks, u2);
                    printf("\n%% removing compressive COD from %s.vtk", file.str().c_str());
                }
				
                //testing ... if(count==33) vdbLoader.writeGrids(subfile.str()+"_recovered.vdb");
                
                //interpolate displacements between full steps
				/*_quad_**
				double x= (subStep+1.0)/(subSteps.getValue()+2.0),
					f1 = 0.5*(x-1.0)*(x-2.0),
					f2 = x*(2.0-x),
					f3 = 0.5*x*(x-1.0);
				/*/
                double f = (subStep+1.0)/(subSteps.getValue()+2.0);
				/**/

                // write output for the sub-step
				/*_quad_**
				printf("\n%% %04d_%04d --> %05d (f2=%.3lf) ... ", fullStep,subStep,count,f2); /*/
                printf("\n%% %04d_%04d --> %05d (%.3lf) ... ", fullStep,subStep,count,f); /**/
                outfile.str(""); outfile.clear();
                outfile << outFile.getValue() << "_";
                outfile.fill('0'); outfile.width(5);
                outfile << (count++); // post-inc!
				/*_quad_**
                vdbLoader.writeVisualMesh(
                    nodes2,elems2,regions2,cracks,
					f1*u + f2*u2 + f3*u3, f1*u_c + f2*u_c2 + f3*u_c3, 
                    outfile.str(), outMeshVDBQual.getValue(), outMeshVCGQual.getValue(),
                    !outMeshNoU.isSet(),!outMeshNoCOD.isSet(),
                    outMeshClose.isSet(),outMeshOBJ.isSet(),seg, outMeshSegments.isSet()
                ); /*/
                vdbLoader.writeVisualMesh(
                    nodes2,elems2,regions2,cracks,
                    (1.0-f)*u + f*u2, (1.0-f)*u_c + f*u_c2,
                    outfile.str(), outMeshVDBQual.getValue(), outMeshVCGQual.getValue(),
                    !outMeshNoU.isSet(),!outMeshNoCOD.isSet(),
                    outMeshClose.isSet(),outMeshOBJ.isSet(),seg, true, outMeshPostQual.getValue()
                ); /**/
				printf(" ok");
            }
        }
        
	}catch (TCLAP::ArgException &e)  {// catch any exceptions
		std::cerr << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return -1;
	}catch(std::exception& e){ printf("exception: %s\n", e.what());
	}catch(...){ printf("caught ellipsis at end of main\n");}
	return 0;
}




int readMeshVTK(
    string filename, node_map& nodes, elem_map& elems, id_map& regions,
    vector_type& u, vector_type& u_c
){
    ifstream in(filename.c_str());
    if(!in.is_open()){
        printf("\n Failed to open file %s", filename.c_str());
        return -1;
    }
    nodes.clear(); elems.clear(); regions.clear();

    int r;
    unsigned int nNodes, nElems, elemType;
    std::istringstream strstream;
    string line, token;
    bool done=false;
    while(in.good() && !done){ // read points header
        getline(in, line);
        if( line.compare(0,6,"POINTS")==0 ){
            r=sscanf(line.substr(7).c_str(), "%u", &nNodes);
            if( r==1 ) done=true;
            else{
                printf("\n Error reading number of points in %s", filename.c_str());
                in.close();
                return -1;
            }
        }
    }
	for(unsigned int n=NODE_BASE_INDEX; n< (nNodes+NODE_BASE_INDEX) && in.good(); ++n){ // read points
        getline(in, line);
        nodes[n].assign(3,0.0);
        r=sscanf(line.c_str(), "%lf %lf %lf", &(nodes[n][0]),&(nodes[n][1]),&(nodes[n][2]) );
        if( r!=3 ){
            printf("\n Error reading point %d in %s", n, filename.c_str());
            in.close();
            return -1;
        }
    }
    done=false;
    while(in.good() && !done){ // read cells header
        getline(in, line);
        if( line.compare(0,5,"CELLS")==0 ){
            r=sscanf(line.substr(6).c_str(), "%u", &nElems);
            if( r==1 ) done=true;
            else{
                printf("\n Error reading number of cells in %s", filename.c_str());
                in.close();
                return -1;
            }
        }
    }
    for(unsigned int n=ELEM_BASE_INDEX; n< (nElems+ELEM_BASE_INDEX) && in.good(); ++n){ // read cells
        getline(in, line);
        elems[n].assign(3,0);
        r=sscanf(line.c_str(), "%u %u %u %u", &elemType, &(elems[n][0]),&(elems[n][1]),&(elems[n][2]) );
        if( r!=4 || elemType!=3){
            printf("\n Error reading cell %d in %s", n, filename.c_str());
            in.close();
            return -1;
        }
		elems[n][0]+=NODE_BASE_INDEX;
		elems[n][1]+=NODE_BASE_INDEX;
		elems[n][2]+=NODE_BASE_INDEX;
    }
    done=false;
    while(in.good() && !done){ // read cells types header
        getline(in, line);
        if( line.compare(0,10,"CELL_TYPES")==0 ) done=true;
    }
    for(unsigned int n=ELEM_BASE_INDEX; n< (nElems+ELEM_BASE_INDEX) && in.good(); ++n){ // check cell types
        getline(in, line);
        if(line!="5"){
            printf("\n Cell %d in %s has wrong type", n, filename.c_str());
            in.close();
            return -1;
        }
    }
    done=false;
    while(in.good() && !done){ // read point data header
        getline(in, line);
        if( line.compare(0,10,"POINT_DATA")==0 ) done=true;
    }
    done=false;
    while(in.good() && !done){ // read point data definition
        getline(in, line);
        strstream.clear();
        strstream.str(line);
        getline(strstream, token, ' ');
        if( token=="VECTORS" ||
            token=="SCALARS" ){
            getline(strstream, token, ' '); // token is now field name
            
            if( token=="displacement"){ // read u
                u.resize(3*nNodes);
                for(unsigned int n=0; n<nNodes && in.good(); ++n){
                    getline(in, line);
                    r=sscanf(line.c_str(), "%lf %lf %lf", &(u[3*n  ]),&(u[3*n+1]),&(u[3*n+2]) );
                    if( r!=3 ){
                        printf("\n Error reading displacement of point %d in %s", n, filename.c_str());
                        in.close();
                        return -1;
                    }
                }
            } else
            if( token=="crack_base_displacement"){ // read u_c
                u_c.resize(3*nNodes);
                for(unsigned int n=0; n<nNodes && in.good(); ++n){
                    getline(in, line);
                    r=sscanf(line.c_str(), "%lf %lf %lf", &(u_c[3*n  ]),&(u_c[3*n+1]),&(u_c[3*n+2]) );
                    if( r!=3 ){
                        printf("\n Error reading displacement of point %d in %s", n, filename.c_str());
                        in.close();
                        return -1;
                    }
                }
            }
            else{ //other data field, skip lines
                for(unsigned int n=0; n<nNodes && in.good(); ++n)
                    getline(in, line);
            }
        }
        if( token=="CELL_DATA") done=true; // found cell data header -- end of point data
    }
    while(in.good()){ // read cell data definition
        getline(in, line);
        strstream.clear();
        strstream.str(line);
        getline(strstream, token, ' ');
        if( token=="VECTORS" ||
            token=="SCALARS" ){
            getline(strstream, token, ' '); // token is now field name
            
            if( token=="region"){ // read region
                getline(in, line); // skip line (contains lookup table definition for scalar field)
                for(unsigned int n=ELEM_BASE_INDEX; n< (nElems+ELEM_BASE_INDEX) && in.good(); ++n){
                    getline(in, line);
                    r=sscanf(line.c_str(), "%u", &elemType);
                    if( r!=1){
                        printf("\n Error reading region of cell %d in %s", n, filename.c_str());
                        in.close();
                        return -1;
                    }
                    regions[n]=elemType;
                }
            }
            else{ //other data field, skip lines
                for(unsigned int n=0; n<nElems && in.good(); ++n)
                    getline(in, line);
            }
        }
    }
    in.close();
    return 0;
}

int removeCompressiveCOD(
	node_map& nodes, elem_map& elems, id_map& regions,
    id_set& cracks, vector_type& u
){
	// for all crack surfaces
	// first compute node-normals
	// then project COD to normal
	// and remove normal component from COD if negative

	vect3d_map n; // map node-ID to node-normal

	for(elem_map::iterator it=elems.begin(); it!=elems.end(); ++it){
		if(cracks.count(regions[it->first])!=0){ // if this is a crack element
			// compute face normal
			Eigen::Vector3d // a,b,c are node coordinates in material space; ua,ub,uc are nodal displacements, so world coordinates are a+ua etc.
				a (nodes[it->second[0]][0], nodes[it->second[0]][1], nodes[it->second[0]][2]),
				b (nodes[it->second[1]][0], nodes[it->second[1]][1], nodes[it->second[1]][2]),
				c (nodes[it->second[2]][0], nodes[it->second[2]][1], nodes[it->second[2]][2]),
				fn;
			fn=(a-c).cross(b-c);
			// initialize node normals
			if(n.count(it->second[0])==0) n[it->second[0]]=Eigen::Vector3d(0.0,0.0,0.0);
			if(n.count(it->second[1])==0) n[it->second[1]]=Eigen::Vector3d(0.0,0.0,0.0);
			if(n.count(it->second[2])==0) n[it->second[2]]=Eigen::Vector3d(0.0,0.0,0.0);
			// sum up
			n[it->second[0]]+=fn;
			n[it->second[1]]+=fn;
			n[it->second[2]]+=fn;
		}
	}
	for(vect3d_map::iterator it=n.begin(); it!=n.end(); ++it)
	if( u.size()> (3*(it->first-NODE_BASE_INDEX)+2) ){
		// normalize node normals (it->first is a node-ID, it->second is the node normal) 
		it->second.normalize();
		unsigned int i1=3*(it->first-NODE_BASE_INDEX);
		unsigned int i2=i1+1, i3=i1+2;
		Eigen::Vector3d u_node(u[i1], u[i2], u[i3]);

		if( u_node.dot(it->second) > 0.0 ){ // opening cracks have u oriented opposite the surface normal as we are using outward normals
			u_node -= it->second * u_node.dot(it->second);
			u[i1]=u_node[0];
			u[i2]=u_node[1];
			u[i3]=u_node[2];
		}
	}
	return 0;
}