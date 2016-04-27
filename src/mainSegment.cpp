#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>
#include <openvdb/util/Util.h>
#include <openvdb/tools/LevelSetSphere.h>

#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>
#include <cfloat>
#include <cmath>

#define _DEBUG_SEGMENT_
#include "mySegment.h"

using namespace std;
using namespace openvdb::v2_2_0;

int main(int argc, char* argv[]){ try{
	if(argc<2){
		printf("useage: %s filename [threshold] [notiles]\n",argv[0]);
		printf("filename is the input file in vdb format and must contain grids named\n"
			   " 'objectGrid' and 'crackGrid...'\n"
			   "threshold is used to remove small handles during segmentation\n"
			   " it must be a non-negative float\n"
			   " it is used relative to the voxel-size in the input grid,\n"
			   " so a threshold value of 0.1 will remove handles smaller than 1/10th of a voxel\n"
			   " ideally threshold should be less than 1.0 to avoid artefacts\n"
               " default is 0.0\n"
			   "notiles can be set to 1 to suppress usage of tiles during segmentation\n"
			   " use this only if you're certain that all crack touch another surface\n"
		);
		return -1;
	}

	string filename = string(argv[1]);
	double thresholdInput=0.0;
	if(argc>2) sscanf(argv[2],"%lf",&thresholdInput);
	bool useTiles=true; int tilesInput=0;
	if(argc>3) sscanf(argv[3],"%d",&tilesInput);
	if(tilesInput==1) useTiles=false;
	
	double narrowBandHalfWidth=3;
	initialize();
	GridPtrVec grids;
	io::File fileInput(filename);
	fileInput.open();
	GridBase::Ptr baseGrid;
	baseGrid = fileInput.readGrid("objectGrid");
	FloatGrid::Ptr objectGrid = gridPtrCast<FloatGrid>(baseGrid);
	vector<FloatGrid::Ptr> crackGrids;
	for(io::File::NameIterator nameIter = fileInput.beginName();
		nameIter != fileInput.endName(); ++nameIter){
		if(nameIter.gridName().substr(0,9)=="crackGrid"){
			//printf("reading %s as crack\n",nameIter.gridName().c_str());
			baseGrid = fileInput.readGrid(nameIter.gridName());
			crackGrids.push_back(gridPtrCast<FloatGrid>(baseGrid));
		}
	}
	fileInput.close();

	grids.clear();
	grids.push_back(objectGrid->deepCopy());
	grids.back()->setName("copy of input object");

	// intersect object with fractures
	for(int i=0; i<crackGrids.size(); ++i){

		for(FloatTree::RootNodeType::ChildOnIter cit = objectGrid->tree().beginRootChildren(); cit.test(); ++cit){
			crackGrids[i]->getAccessor().touchLeaf(cit.getCoord()); // make sure the crack-grid covers the object grid
			crackGrids[i]->getAccessor().setValueOn(cit.getCoord());
			// if we find a background value, the sign needs to be flipped
			if( crackGrids[i]->getAccessor().getValue(cit.getCoord()) == crackGrids[i]->background()/* > crackGrids[i]->voxelSize()[0]*/ ){
				crackGrids[i]->getAccessor().setValue(cit.getCoord(), -crackGrids[i]->getAccessor().getValue(cit.getCoord()));
			}
		}
		crackGrids[i]->signedFloodFill();

		tools::csgIntersection(*objectGrid, *crackGrids[i]); // don't use negative background values!
	}
	grids.push_back(objectGrid->deepCopy());
	stringstream name; name << "object after fracture (" << objectGrid->memUsage()/1024.0/1024.0 << " MB)";
	grids.back()->setName(name.str());
	printf("%s\n",name.str().c_str());

	double voxelSize = objectGrid->voxelSize()[0]; // assume voxels are the same size in all directions
	double handleThreshold = thresholdInput*voxelSize; // 0.0 reproduces the surface perfectly, but keeps all handles connected, higher values remove handles but cause extrapolation artefacts -- maybe these can be avoided by accessing the voxels of interest in the original grid (if they are outside)
	printf("objectGrid after pre-process has %.3lf MB\n", objectGrid->memUsage()/1024.0/1024.0);

	// segment ...
	vector<FloatGrid::Ptr> segments = mySegment(*objectGrid,handleThreshold,useTiles);

	printf("objectGrid after segmentation has %.3lf MB\n", objectGrid->memUsage()/1024.0/1024.0);

	// write output
	for(int i=0; i<segments.size(); ++i){
		name.str(""); name.clear();
		name << "segment " << std::setw(7) << std::setfill('0') << (i+1);
		name << " of " << segments.size() << " (" << segments[i]->memUsage()/1024.0/1024.0 << " MB)";
		grids.push_back(segments[i]);
		grids.back()->setName(name.str());
		printf("%s\n",name.str().c_str());
	}
    
	stringstream outfn;
	outfn << filename.substr(0, filename.find_last_of('.'));
	outfn << "_segmented" << filename.substr(filename.find_last_of('.'));
	io::File fileResult(outfn.str());
	fileResult.write(grids);
	fileResult.close();

	}catch(std::exception& e){ printf("exception: %s\n", e.what());
	}catch(...){ printf("caught ellipsis at end of main\n");}
}
