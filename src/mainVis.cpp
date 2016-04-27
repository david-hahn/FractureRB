#include "types.h"
#include "VDBLoader.h"
#include "FractureRB_config.h"

#include <openvdb/tools/LevelSetMeasure.h>
#include <openvdb/tools/LevelSetRebuild.h>

#include "vcgHelper.h"
#include <vcg/wrap/io_trimesh/export_stl.h>

#include <tclap/CmdLine.h>
#include <string>
#include <sstream>
#include <algorithm>

//#include <openvdb/tools/LevelSetTracker.h> // for testing

using namespace FractureSim;
using namespace std;
using namespace myVCG;

// returns the name of the file written (or empty string on error)
string writeMesh(vdb::FloatGrid::Ptr grid, string file, double vdbQual, bool asObj, bool isPrimary=false);

void writeUpdateCommands(ofstream& out, const string& name, const string& file);

inline bool file_exists(const string& name){
    if( FILE *file = fopen(name.c_str(), "r")){
        fclose(file);
        return true;
    }else{
        return false;
    }   
}

inline std::string& myReplace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
        return str;
    str.replace(start_pos, from.length(), to);
    return str;
}

int main(int argc, char* argv[]){
	try {
		printf("\n%% parsing command line ... ");
		TCLAP::CmdLine cmd("Post process FractureRB results to OBJ or STL visual quality meshes", ' ', "using OpenVDB 2.2 backend");
		// input params
		TCLAP::ValueArg<string> inFile("i","in","prefix for input files", true, "","in", cmd);
		TCLAP::ValueArg<double> inSteps("n","steps","number of input time-steps to read", false, 1000,"steps", cmd);
		TCLAP::UnlabeledMultiArg<string> objNames("object",
			"specify the names of all objects to be processed",
			true, "name", cmd
		);
		TCLAP::SwitchArg vdbRebuild("", "vdb-rebuild", "rebuild level-set representation for huge fragments after segmentation (use if number of fragments in hi-res output does not match sim results exactly", cmd);
		//TCLAP::ValueArg<int> fullSteps("s","steps","number of crack propagation steps", false, 0,"steps", cmd);
		//TCLAP::ValueArg<int> subSteps ("S","substeps","number of crack propagation substeps", false, -1,"substeps", cmd);
		//TCLAP::ValueArg<int> startStep("x","start","start at the given full step (skip previous ones, default 0)", false, 0,"substeps", cmd);
		//TCLAP::ValueArg<double> resMesh("r","res-mesh","resolution of BEM mesh as in the simulation, used for nearest triangle search, default is 0.1", false, 0.1,"res-mesh", cmd);

		//output params
        TCLAP::ValueArg<string> outFile("o","out","postfix for output files", true, "","out", cmd);
		TCLAP::ValueArg<double> outMeshVDBQual("q","vdb-qual","quality parameter for VDB meshing in [0,1] (default is 0.01)", false, 0.01,"vdb-qual", cmd);
		TCLAP::SwitchArg noPrimMesh("", "no-prim-mesh", "suppress mesh output for primary objects", cmd);
		//TCLAP::ValueArg<double> outMeshVCGQual("Q","vis-qual","quality parameter for mesh decimation: quadric error tolerance (default is no decimation)", false, -1.0,"max-error", cmd);
		//TCLAP::ValueArg<double> outMeshPostQual("P","post-qual","quality parameter for post-interpolation mesh decimation: quadric error tolerance (default is no decimation)", false, -1.0,"max-error", cmd);
		//TCLAP::ValueArg<double> outMeshSegment("","segment","if thr>0.0 segment before meshing with handle-removal, default: no segmentation, use thr<=1.0", false, -1.0,"thr", cmd);
		//TCLAP::ValueArg<double> outLastSegment("","segment-last","same as segment, but applies only to last (full step) output, default: no segmentation", false, -1.0,"thr", cmd);
		//TCLAP::SwitchArg outMeshNoU("", "vis-no-disp", "do NOT apply any displacements when producing hi-res output", cmd);
		//TCLAP::SwitchArg outMeshNoCOD("", "vis-no-cod", "do NOT apply crack opening displacements when producing hi-res output", cmd);
		TCLAP::SwitchArg outMeshClose("", "vis-close", "make crack faces coincident for 0 COD in hi-res output", cmd);
		//TCLAP::SwitchArg outMeshAllowNegCOD("", "vis-allow-neg-cod", "allow compressive crack opening displacement (may lead to mesh fold-over), otherwise compressive COD will be projected out of the visual mesh", cmd);
		TCLAP::SwitchArg outMeshOBJ("", "vis-obj", "write hi-res output as .OBJ file (otherwise as .STL)", cmd);
		TCLAP::SwitchArg closeCracks("", "close-cracks", "close cracks such that fragments appear touching in material space", cmd);
		//TCLAP::SwitchArg outMeshSegments("", "segment-files", "write one file per segment instead of all segments in a single file", cmd);
		cmd.parse( argc, argv );

        int ret=0; bool doUpdate;
        stringstream file, outfile, postfix;
		
		vector<string> objects = objNames.getValue();
		vector<double> volumes; volumes.assign(objects.size(), 0.0);
		vector<int> baseIds;    baseIds.assign(objects.size(), -1);
		id_map updateCount, fragmentCount;
		int nPrimaryObjs = objects.size();
		map<string,string> meshForName; // map object names (in the scene) to file names containing their hi-res meshes
        string meshFile;
		
        VDBLoader vdbLoader;
		id_set    cracks;
		map<int, vdb::FloatGrid::Ptr> baseGrids;
        
		for(int id=0; id<objects.size(); ++id){
			printf("\n\n%% processing \"%s\" ... ", objects[id].c_str());
			// first load the base file
			file.str(""); file.clear();
			file << inFile.getValue() << objects[id] << ".vdb";
			outfile.str(""); outfile.clear();
			outfile << inFile.getValue() << objects[id] << outFile.getValue();

			if( file_exists(file.str()) ){
				printf("\n%% base file %s found",file.str().c_str());
				
				cracks.clear();
				ret=vdbLoader.loadFull(file.str(), cracks);
				if( ret<0) printf("\n%% Failed to read base file %s, loader returned %d\n", file.str().c_str(), ret);

				printf(" ... have %d cracks", cracks.size());

				updateCount[id]  =0;
				fragmentCount[id]=0;
				if( id<nPrimaryObjs ){
					volumes[id]  = vdb::tools::levelSetVolume(*vdbLoader.getObjectGrid());
					baseGrids[id]= vdbLoader.getObjectGrid()->deepCopy();
					baseIds[id]  = id; // this is its own base object
				}
				
				if( !noPrimMesh.getValue() && id<nPrimaryObjs ){ // write output for primary mesh
					meshFile=writeMesh(vdbLoader.getObjectGrid(),outfile.str()+"_0",outMeshVDBQual.getValue(),outMeshOBJ.getValue(),true);
					if(!meshFile.empty()) meshForName[objects[id]]=meshFile;
				}
			}
			// now check for updates
			for(int step=0; ret==0 && step<=inSteps.getValue(); ++step){
				doUpdate=false;
				file.str(""); file.clear();
				file << inFile.getValue() << objects[id] << "_" << step << ".vdb";
				if( file_exists(file.str()) ){
					printf("\n%% update file %s found",file.str().c_str());
					ret=vdbLoader.loadUpdate(file.str(), cracks);
					if( ret<0){
						printf("\n%% Failed to read update %s, loader returned %d\n", file.str().c_str(), ret);
					}else{
						printf(" ... have %d cracks", cracks.size());
						
						//implemented analogous to FractureRB_fragments -- could maybe move into common function?
						if( cracks.size() > 0 ){
							vector<vdb::FloatGrid::Ptr> segments = vdbLoader.getSegments(SEG_HDL_THR, false /*use tiles*/);
							if( segments.size() > 1){
								double parentVolume = vdb::tools::levelSetVolume(*vdbLoader.getObjectGrid());
								if( parentVolume<0.0 ) parentVolume*=-1.0;
								double segmentVolume;
								for(unsigned int i=0; i<segments.size(); ++i){
									// decide the proper postfix (either _f#, _F#, or _#) for this fragment
									postfix.str(""); postfix.clear();
									bool ignore=false;
									segmentVolume = vdb::tools::levelSetVolume(*segments[i]);
									if( segmentVolume<0.0 ) segmentVolume*=-1.0;
									
									if( segmentVolume <= FRAGMENT_IGNORE_THR* volumes[id] ){
										// do nothing
										ignore=true;
									}else
									if( segmentVolume <= FRAGMENT_SMALL_THR* volumes[id] ){ //parentVolume ){ // small fragment
										postfix << "_f" << ++fragmentCount[id];
									}else if( segmentVolume >= FRAGMENT_ORI_THR* parentVolume ){// replace the implicit surface of the parent
										postfix << "_" << ++updateCount[id];
										if( vdbRebuild.getValue() ){
											printf("\n%% rebuilding level set on update (%s)",(objects[id]+postfix.str()).c_str() );
											segments[i]=vdb::tools::levelSetRebuild(*segments[i],0.0,vdbLoader.getBandwidth()); // make sure the narrow-band is properly resolved (also removes most interior cracks)
											segments[i]->setName("objectGrid");
											segments[i]->setGridClass(vdb::GRID_LEVEL_SET);
										}
										vdbLoader.setObjectGrid( segments[i] );
										doUpdate=true;
									}else{ // large fragment
										postfix << "_F" << ++fragmentCount[id];
										objects.push_back( objects[id]+postfix.str() );
										volumes.push_back( volumes[id] );
										baseIds.push_back( baseIds[id] );
									}

									if(!ignore){
										if( closeCracks.isSet() ){
											//to make cracks appear (almost) closed: dilate segment grid by halfDiag*voxelSize then intersect with the object grid
											segments[i]=vdb::tools::levelSetRebuild(*segments[i],0.0,vdbLoader.getBandwidth()); // make sure the narrow-band is properly resolved (also removes most interior cracks)
											vdbLoader.intersectWithCrackGrids(segments[i]);
											for(vdb::FloatGrid::ValueOnIter it= segments[i]->beginValueOn(); it.test(); ++it){
												it.setValue(it.getValue()-vdbLoader.getVoxelSize()*halfDiag);
											}
											vdb::tools::csgIntersection(*segments[i], *(baseGrids[baseIds[id]]->deepCopy()));
										}

										meshFile=writeMesh(segments[i],(outfile.str()+postfix.str()),outMeshVDBQual.getValue(),outMeshOBJ.getValue());
										if(!meshFile.empty()) meshForName[objects[id]+postfix.str()]=meshFile;
										
										
//										if(0){//testing 
//											vdb::FloatGrid::Ptr testGrid=segments[i]->deepCopy();
//											for(vdb::FloatGrid::ValueAllIter vit = testGrid->beginValueAll(); vit.test(); ++vit ){
//												vit.setValue(*vit-testGrid->voxelSize()[0]); // dilate surface by a voxel
//											}
//////											// test tracking
////											vdb::tools::LevelSetTracker<vdb::FloatGrid> tracker(*testGrid);
////											tracker.track(); // rebuild the narrow-band
//											testGrid=vdb::tools::levelSetRebuild(*testGrid,0.0,vdbLoader.getBandwidth());
//											for(vdb::FloatGrid::ValueAllIter vit = testGrid->beginValueAll(); vit.test(); ++vit ){
//												vit.setValue(*vit+testGrid->voxelSize()[0]); // erode surface by a voxel
//											}
//											// write file for debugging ...
//											vdb::GridPtrVec grids;
//											grids.push_back(testGrid);
//											vdb::io::File file("debug_tracking.vdb");
//											file.write(grids);
//											file.close();
//										}
									}
								}
								if( doUpdate ){
									postfix.str(""); postfix.clear();
									postfix << "_" << updateCount[id]; // no change to counter here!
									objects[id].append(postfix.str());
									outfile.str(""); outfile.clear();
									outfile << inFile.getValue() << objects[id] << outFile.getValue();
								}
							}
						}
					}
				}
			}
        }
		
		file.str(""); file.clear();
		file << inFile.getValue() << "vis.mel";
		printf("\n writing update script (%s) ... ",file.str().c_str());
		ofstream outMel(file.str().c_str());
		printf("\n\n%% object names --> new mesh file:");
		for(map<string,string>::iterator it=meshForName.begin(); it!=meshForName.end(); ++it){
			printf("\n%% %s --> %s", it->first.c_str(), it->second.c_str());
			//int id = objects.size();
			//vector<string>::iterator oit;
			//if( (oit = find(objects.begin(),objects.end(),it->first))!=objects.end() ) id=objects.end()-oit;
			//printf(" id is %d", id);
			string chkFile(it->second);
			if( (find(objNames.getValue().begin(),objNames.getValue().end(),it->first)==objNames.getValue().end() ) &&
				!file_exists( myReplace(chkFile, outFile.getValue(), "") ) ){
				printf("\n!!! no low-res counterpart found for %s", it->second.c_str());
			}else
			if( outMeshOBJ.getValue()) //ToDo: support STL here
				writeUpdateCommands(outMel, it->first, it->second);
		}
		outMel.close();
	}catch (TCLAP::ArgException &e)  {// catch any exceptions
		std::cerr << std::endl << "error: " << e.error() << " for arg " << e.argId() << std::endl;
		return -1;
	}catch(std::exception& e){ printf("\nexception: %s\n", e.what());
	}catch(...){ printf("\ncaught ellipsis at end of main\n");}
	return 0;
}

void writeUpdateCommands(ofstream& out, const string& name, const string& file){
	//ToDo: support STL here as well
	out << "file -import -type \"OBJ\" -ignoreVersion -mergeNamespacesOnClash true \"" << file << "\";" << endl;
	out << "select -r |Mesh ;" << endl;
	out << "select -add " << name << " ;" << endl;
	out << "parent -r;" << endl;
	out << "string $shp[] = `listRelatives -shapes \"" << name << "\"`;" << endl;
	out << "setAttr ($shp[0]+\".visibility\") 0;" << endl;
}

string writeMesh(vdb::FloatGrid::Ptr grid, string file, double vdbQual, bool asObj, bool isPrimary){
	vdb::tools::VolumeToMesh* mesher;
	vcg::tri::io::ExporterOBJ<VcgMesh> writerObj;
	vcg::tri::io::ExporterSTL<VcgMesh> writerStl;
	VcgMesh mesh;

	mesher = new vdb::tools::VolumeToMesh(0.0,vdbQual);
	(*mesher)(*grid);
	mesh.Clear();
	copyVDBtoVCG(*mesher,mesh);
	delete mesher; mesher=NULL;

	//ToDo: apply crack-closing-displacements as in VDBWrapper_mesh

	if(!isPrimary){
		//COM correction (don't do this for primary objects, since they are allowed to have non-COM coordinate origins)
		double volume=0.0, vol;
		Eigen::Vector3d a,b,c, com;
		com.setZero();
		for(int i=0; i<mesh.face.size(); ++i) if(!mesh.face[i].IsD()){
			mesh.face[i].V(0)->P().ToEigenVector(a);
			mesh.face[i].V(1)->P().ToEigenVector(b);
			mesh.face[i].V(2)->P().ToEigenVector(c);
			vol = a.dot( b.cross( c )) / 6.0;
			volume += vol;
			com += vol*( a+b+c )/4.0; // volume weighted centroid of tet (a,b,c,0)
		}
		com/=volume;
		for(int i=0; i<mesh.vert.size(); ++i) if(!mesh.vert[i].IsD()){
			mesh.vert[i].P().ToEigenVector(a);
			a-=com;
			mesh.vert[i].P().FromEigenVector(a);
		}
	}

	printf("\n%% writing %s%s",file.c_str(),asObj?".obj":".stl");
	if( asObj ){
		vcg::tri::UpdateNormal<VcgMesh>::PerVertex(mesh); // compute area-weighted normals
		vcg::tri::UpdateNormal<VcgMesh>::NormalizePerVertex(mesh); // make them unit-vectors
		int r=writerObj.Save(mesh,(file+".obj").c_str(), vcg::tri::io::Mask::IOM_VERTNORMAL);
		if(r==0) return (file+".obj");
	}else{
		int r=writerStl.Save(mesh,(file+".stl").c_str(),true);
		if(r==0) return (file+".stl");
	}
	return ""; // something went wrong if we reach this
}


//int readMeshVTK(
//    string filename, node_map& nodes, elem_map& elems, id_map& regions,
//    vector_type& u, vector_type& u_c
//){
//    ifstream in(filename.c_str());
//    if(!in.is_open()){
//        printf("\n Failed to open file %s", filename.c_str());
//        return -1;
//    }
//    nodes.clear(); elems.clear(); regions.clear();
//
//    int r;
//    unsigned int nNodes, nElems, elemType;
//    std::istringstream strstream;
//    string line, token;
//    bool done=false;
//    while(in.good() && !done){ // read points header
//        getline(in, line);
//        if( line.compare(0,6,"POINTS")==0 ){
//            r=sscanf(line.substr(7).c_str(), "%u", &nNodes);
//            if( r==1 ) done=true;
//            else{
//                printf("\n Error reading number of points in %s", filename.c_str());
//                in.close();
//                return -1;
//            }
//        }
//    }
//	for(unsigned int n=NODE_BASE_INDEX; n< (nNodes+NODE_BASE_INDEX) && in.good(); ++n){ // read points
//        getline(in, line);
//        nodes[n].assign(3,0.0);
//        r=sscanf(line.c_str(), "%lf %lf %lf", &(nodes[n][0]),&(nodes[n][1]),&(nodes[n][2]) );
//        if( r!=3 ){
//            printf("\n Error reading point %d in %s", n, filename.c_str());
//            in.close();
//            return -1;
//        }
//    }
//    done=false;
//    while(in.good() && !done){ // read cells header
//        getline(in, line);
//        if( line.compare(0,5,"CELLS")==0 ){
//            r=sscanf(line.substr(6).c_str(), "%u", &nElems);
//            if( r==1 ) done=true;
//            else{
//                printf("\n Error reading number of cells in %s", filename.c_str());
//                in.close();
//                return -1;
//            }
//        }
//    }
//    for(unsigned int n=ELEM_BASE_INDEX; n< (nElems+ELEM_BASE_INDEX) && in.good(); ++n){ // read cells
//        getline(in, line);
//        elems[n].assign(3,0);
//        r=sscanf(line.c_str(), "%u %u %u %u", &elemType, &(elems[n][0]),&(elems[n][1]),&(elems[n][2]) );
//        if( r!=4 || elemType!=3){
//            printf("\n Error reading cell %d in %s", n, filename.c_str());
//            in.close();
//            return -1;
//        }
//		elems[n][0]+=NODE_BASE_INDEX;
//		elems[n][1]+=NODE_BASE_INDEX;
//		elems[n][2]+=NODE_BASE_INDEX;
//    }
//    done=false;
//    while(in.good() && !done){ // read cells types header
//        getline(in, line);
//        if( line.compare(0,10,"CELL_TYPES")==0 ) done=true;
//    }
//    for(unsigned int n=ELEM_BASE_INDEX; n< (nElems+ELEM_BASE_INDEX) && in.good(); ++n){ // check cell types
//        getline(in, line);
//        if(line!="5"){
//            printf("\n Cell %d in %s has wrong type", n, filename.c_str());
//            in.close();
//            return -1;
//        }
//    }
//    done=false;
//    while(in.good() && !done){ // read point data header
//        getline(in, line);
//        if( line.compare(0,10,"POINT_DATA")==0 ) done=true;
//    }
//    done=false;
//    while(in.good() && !done){ // read point data definition
//        getline(in, line);
//        strstream.clear();
//        strstream.str(line);
//        getline(strstream, token, ' ');
//        if( token=="VECTORS" ||
//            token=="SCALARS" ){
//            getline(strstream, token, ' '); // token is now field name
//            
//            if( token=="displacement"){ // read u
//                u.resize(3*nNodes);
//                for(unsigned int n=0; n<nNodes && in.good(); ++n){
//                    getline(in, line);
//                    r=sscanf(line.c_str(), "%lf %lf %lf", &(u[3*n  ]),&(u[3*n+1]),&(u[3*n+2]) );
//                    if( r!=3 ){
//                        printf("\n Error reading displacement of point %d in %s", n, filename.c_str());
//                        in.close();
//                        return -1;
//                    }
//                }
//            } else
//            if( token=="crack_base_displacement"){ // read u_c
//                u_c.resize(3*nNodes);
//                for(unsigned int n=0; n<nNodes && in.good(); ++n){
//                    getline(in, line);
//                    r=sscanf(line.c_str(), "%lf %lf %lf", &(u_c[3*n  ]),&(u_c[3*n+1]),&(u_c[3*n+2]) );
//                    if( r!=3 ){
//                        printf("\n Error reading displacement of point %d in %s", n, filename.c_str());
//                        in.close();
//                        return -1;
//                    }
//                }
//            }
//            else{ //other data field, skip lines
//                for(unsigned int n=0; n<nNodes && in.good(); ++n)
//                    getline(in, line);
//            }
//        }
//        if( token=="CELL_DATA") done=true; // found cell data header -- end of point data
//    }
//    while(in.good()){ // read cell data definition
//        getline(in, line);
//        strstream.clear();
//        strstream.str(line);
//        getline(strstream, token, ' ');
//        if( token=="VECTORS" ||
//            token=="SCALARS" ){
//            getline(strstream, token, ' '); // token is now field name
//            
//            if( token=="region"){ // read region
//                getline(in, line); // skip line (contains lookup table definition for scalar field)
//                for(unsigned int n=ELEM_BASE_INDEX; n< (nElems+ELEM_BASE_INDEX) && in.good(); ++n){
//                    getline(in, line);
//                    r=sscanf(line.c_str(), "%u", &elemType);
//                    if( r!=1){
//                        printf("\n Error reading region of cell %d in %s", n, filename.c_str());
//                        in.close();
//                        return -1;
//                    }
//                    regions[n]=elemType;
//                }
//            }
//            else{ //other data field, skip lines
//                for(unsigned int n=0; n<nElems && in.good(); ++n)
//                    getline(in, line);
//            }
//        }
//    }
//    in.close();
//    return 0;
//}
