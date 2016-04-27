/* 
 * File:   VDBLoader.cpp
 * Author: David
 * 
 * Created on 16. Juli 2014
 */

#include "VDBLoader.h"

namespace FractureSim{
    using namespace vdb;
    int VDBLoader::loadFull(std::string filename, id_set& cracks){
        //printf("%% loading full step from %s\n", filename.c_str());
        objectGrid.reset();
        resXform.reset();
        crackGrids.clear();
        signedCrackGrids.clear();
        
        unsigned int crackId;
        GridPtrVec grids;
        io::File fileInput(filename);
        fileInput.open();
        GridBase::Ptr baseGrid;

        baseGrid = fileInput.readGrid("objectGrid");
        objectGrid = gridPtrCast<FloatGrid>(baseGrid);

        voxelSize = objectGrid->voxelSize()[0];
        resXform  = objectGrid->transformPtr();

        for(io::File::NameIterator nameIter = fileInput.beginName();
            nameIter != fileInput.endName(); ++nameIter){
            if(nameIter.gridName().substr(0,10)=="crackGrid_"){
                baseGrid = fileInput.readGrid(nameIter.gridName());
                sscanf(nameIter.gridName().substr(10).c_str(), "%u", &crackId);
                crackGrids[crackId]=gridPtrCast<FloatGrid>(baseGrid);
//                printf("%% ... added crack grid %u\n", crackId);
                cracks.insert(crackId);
            }
            if(nameIter.gridName().substr(0,16)=="signedCrackGrid_"){
                baseGrid = fileInput.readGrid(nameIter.gridName());
                sscanf(nameIter.gridName().substr(16).c_str(), "%u", &crackId);
                signedCrackGrids[crackId]=gridPtrCast<FloatGrid>(baseGrid);
//                printf("%% ... added signed crack grid %u\n", crackId);
            }
        }
        fileInput.close();
        return 0;
    }
    
    int VDBLoader::loadUpdate(std::string filename, id_set& cracks){
        //printf("\n%% loading sub step from %s", filename.c_str());

        unsigned int crackId;
        GridPtrVec grids;
        io::File fileInput(filename);
        fileInput.open();
        GridBase::Ptr baseGrid;
		//printf(" ... opened");
        for(io::File::NameIterator nameIter = fileInput.beginName();
            nameIter != fileInput.endName(); ++nameIter){
			//printf("\n%% processing grid %s", nameIter.gridName().c_str());
            if(nameIter.gridName().substr(0,10)=="crackGrid_"){
                baseGrid = fileInput.readGrid(nameIter.gridName());
				//printf(" ... read ok");
                sscanf(nameIter.gridName().substr(10).c_str(), "%u", &crackId);
				//printf(" ID is %u", crackId);
                if(crackGrids.count(crackId)==0){
                    crackGrids[crackId]=gridPtrCast<FloatGrid>(baseGrid);
                    //printf(" ... new crack (%u)", crackId);
                    cracks.insert(crackId);
                }else{ // combine
                    tools::csgIntersection(
                        *crackGrids[crackId],
                        *gridPtrCast<FloatGrid>(baseGrid)
                    );
					//printf(" ... merged crack (%u)", crackId);
                }
            }
            if(nameIter.gridName().substr(0,16)=="signedCrackGrid_"){
                baseGrid = fileInput.readGrid(nameIter.gridName());
                sscanf(nameIter.gridName().substr(16).c_str(), "%u", &crackId);
                if(signedCrackGrids.count(crackId)==0){
                    signedCrackGrids[crackId]=gridPtrCast<FloatGrid>(baseGrid);
				//printf(" ... new signed crack (%u)", crackId);
				}else{ // combine
                    vdbToolsExt::compAbsMin(
                        *signedCrackGrids[crackId],
                        *gridPtrCast<FloatGrid>(baseGrid)
                    );
					//printf(" ... merged signed crack (%u)", crackId);
                }
            }
        }
        fileInput.close();
        return 0;
    }

	int VDBLoader::intersectWithCrackGrids(vdb::FloatGrid::Ptr grid){
		vdb::FloatGrid::Ptr currentCrack;
		for(std::map<unsigned int, vdb::FloatGrid::Ptr>::const_iterator it=crackGrids.begin(); it!=crackGrids.end(); ++it){
			currentCrack = it->second->deepCopy();
			vdb::tools::csgIntersection(*grid, *currentCrack );
		}
	return 0;
	}
}