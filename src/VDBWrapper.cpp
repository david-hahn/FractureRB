/* 
 * File:   VDBWrapper.cpp
 * Author: David
 * 
 * Created on 14. JÃ¤nner 2014, 13:47
 */

#include "VDBWrapper.h"

#include <iostream>
#include <sstream>
#include <cfloat>

#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/ValueTransformer.h>
#include <openvdb/math/Operators.h>
#include <openvdb/tools/LevelSetRebuild.h>

using namespace std;
using namespace vdb;

namespace FractureSim{
	const vdb::Vec3d VDBWrapper::ux(1.0,0.0,0.0), VDBWrapper::uy(0.0,1.0,0.0), VDBWrapper::uz(0.0,0.0,1.0);

	// Functor to convert a crack-UDF to a level-set modelling a thin sheet of "dead" material
	struct crackUDFtoLS{
		const double voxelSize_;
		crackUDFtoLS(const double voxelSize) : voxelSize_(voxelSize) {}
		inline void operator()(const FloatGrid::ValueAllIter& it) const{
			it.setValue(halfDiag*voxelSize_ - std::abs(*it));
		}
	};

    int VDBWrapper::init(
        node_map& nodes, elem_map& elems,
        id_map& regions, id_set& cracks
    ){
        // 1. build the objectGrid from all regions that are not cracks
        //    -- these must form a closed object
        // 2. build the crackGrids from each region that is a crack
        // 3. build the closestTri grid from ALL regions
        vector<Vec3s> pointList;
        vector<Vec4I> objectList; // build list of trianlges on the object's boundary
        map<unsigned int, vector<Vec4I> > crackLists;   // build list of triangles on cracks
        vector<Vec4I> allTriList;
        
        for(node_map::iterator it = nodes.begin(); it != nodes.end(); ++it)
            pointList.push_back(resXform->worldToIndex(
                Vec3s(it->second[0],it->second[1],it->second[2]
            )));
        for(elem_map::iterator it = elems.begin(); it != elems.end(); ++it){
            if( cracks.count(regions[it->first]) >0 ){
                crackLists[regions[it->first]].push_back(Vec4I(
                    it->second[0]-NODE_BASE_INDEX,
                    it->second[1]-NODE_BASE_INDEX,
                    it->second[2]-NODE_BASE_INDEX,
                    util::INVALID_IDX
                ));
                //printf("elem %d into crackLists of region %d\n",it->first,regions[it->first]);
            }else{
                objectList.push_back(Vec4I(
                    it->second[0]-NODE_BASE_INDEX,
                    it->second[1]-NODE_BASE_INDEX,
                    it->second[2]-NODE_BASE_INDEX,
                    util::INVALID_IDX
                ));
                //printf("elem %d into objectList -- region %d\n",it->first,regions[it->first]);
            }
            allTriList.push_back(Vec4I(
                it->second[0]-NODE_BASE_INDEX,
                it->second[1]-NODE_BASE_INDEX,
                it->second[2]-NODE_BASE_INDEX,
                util::INVALID_IDX
            ));
        }

        tools::MeshToVolume<FloatGrid> voxelizer(resXform,tools::GENERATE_PRIM_INDEX_GRID); // GENERATE_PRIM_INDEX_GRID required for signedCrackGrids below
        voxelizer.convertToLevelSet(pointList, objectList, NBHW, NBHW);
        objectGrid = voxelizer.distGridPtr()->deepCopy();
        objectGrid->setName("objectGrid");
        
        crackGrids.clear(); signedCrackGrids.clear(); 
        for(set<unsigned int>::iterator it = cracks.begin(); it != cracks.end(); ++it){
            voxelizer.clear();
            voxelizer.convertToUnsignedDistanceField(pointList, crackLists[*it], NBHW+1); // add 1 voxel because we dilate later

			signedCrackGrids[*it] = voxelizer.distGridPtr()->deepCopy();

			// need to do something like a tile-topology union to activate tiles in the crack-grid such that the object is covered
			//printf("\nroot loop");
			for(FloatTree::RootNodeType::ChildOnIter cit = objectGrid->tree().beginRootChildren(); cit.test(); ++cit){
				voxelizer.distGridPtr()->getAccessor().touchLeaf(cit.getCoord()); // make sure the crack-grid covers the object grid
				voxelizer.distGridPtr()->getAccessor().setValueOn(cit.getCoord());
			}

#ifdef _MSC_VER
            tools::foreach(voxelizer.distGridPtr()->beginValueAll(),crackUDFtoLS(voxelSize) ); // subtract 0.87*voxelSize and flip sign
#else // there's a strange problem when using tools::foreach with MinGW, maybe due to the tbb-lib
			for(FloatGrid::ValueAllIter vit = voxelizer.distGridPtr()->beginValueAll(); vit.test(); ++vit ){
				vit.setValue(-std::abs(*vit) + halfDiag*voxelSize);
			}
#endif
			voxelizer.distGridPtr()->signedFloodFill(); // propagate new sign to tiles with off-voxels

			crackGrids[*it] = voxelizer.distGridPtr()->deepCopy();
			crackGrids[*it]->setGridClass(GRID_LEVEL_SET);
            stringstream tmp; tmp << *it;
            crackGrids[*it]->setName("crackGrid_" + tmp.str());
			//tmp << ".txt";
			//dumpGrid(*voxelizer.distGridPtr(),tmp.str());

			signedCrackGrids[*it]->setGridClass(GRID_LEVEL_SET);
			signedCrackGrids[*it]->setName("signedCrackGrid_" + tmp.str());
			//signedCrackGrids[*it]->setBackground(0.1*voxelSize);
			//ToDo: ? make this a functor and use tools::foreach ?
			Coord c; Vec3s v,n; unsigned int tri;
			Int32Grid::ConstAccessor triAcc = voxelizer.indexGridPtr()->getConstAccessor();
			for(FloatGrid::ValueOnIter vit = signedCrackGrids[*it]->beginValueOn(); vit.test(); ++vit){
				c = vit.getCoord();
				tri = triAcc.getValue(c);
				v = (( pointList[crackLists[*it][tri][0]] + pointList[crackLists[*it][tri][1]] + pointList[crackLists[*it][tri][2]] )/3.0) - c.asVec3s();
				n = ( pointList[crackLists[*it][tri][1]] - pointList[crackLists[*it][tri][0]] ).cross( pointList[crackLists[*it][tri][2]] - pointList[crackLists[*it][tri][0]] );
				if( v.dot(n) < 0.0 ) // see if (c - triangle midpoint ) is positive or negative projected onto the normal of the triangle
					vit.setValue(-vit.getValue());
			}
        }
        voxelizer.clear();
        
//        tools::MeshToVolume<FloatGrid> voxelizer2(resXform,tools::GENERATE_PRIM_INDEX_GRID);
//        voxelizer2.convertToUnsignedDistanceField(pointList, allTriList, NBHW+1); // add 1 voxel because we dilate later
//#ifdef _MSC_VER
//		tools::foreach(voxelizer2.distGridPtr()->beginValueAll(), crackUDFtoLS(voxelSize) ); // subtract 0.87*voxelSize
//#else // there's a strange problem when using tools::foreach with MinGW, maybe due to the tbb-lib
//			for(FloatGrid::ValueAllIter it = voxelizer2.distGridPtr()->beginValueAll(); it.test(); ++it ){
//				it.setValue(-std::abs(*it) + halfDiag*voxelSize);
//			}
//#endif
//		//voxelizer2.distGridPtr()->setBackground(-voxelizer2.distGridPtr()->background()); // we flipped the sign on all values, also flip the background
//		voxelizer2.distGridPtr()->signedFloodFill(); // propagate new sign to tiles with off-voxels
//		//allRegionsGrid = tools::levelSetRebuild(*voxelizer2.distGridPtr(),0.0,NBHW,NBHW);
//		allSurfacesGrid = voxelizer2.distGridPtr()->deepCopy();
//		allSurfacesGrid->setGridClass(GRID_LEVEL_SET);
//        closestTri      = voxelizer2.indexGridPtr()->deepCopy();
//		allSurfacesGrid->setName("allSurfacesGrid");
//        closestTri->setName("closestTri");
//        voxelizer2.clear();

        return 0;
    }

	int VDBWrapper::addCrack(
		node_map& nodes, elem_map& newElems,
		id_map& newRegions, unsigned int ndBaseIndex
	){
		// build VDB compatible point and polygon lists for each region
		// voxelize to UDF
		// convert to thin level-set
		// union(intersect) into existing crack-grids
        vector<Vec3d> pointList;
        map<unsigned int, vector<Vec4I> > crackLists; // build list of triangles on cracks
        for(node_map::iterator it = nodes.begin(); it != nodes.end(); ++it)
            pointList.push_back(resXform->worldToIndex(
                Vec3d(it->second[0],it->second[1],it->second[2]
            )));
        for(elem_map::iterator it = newElems.begin(); it != newElems.end(); ++it){
            crackLists[newRegions[it->first]].push_back(Vec4I(
                it->second[0]-ndBaseIndex,
                it->second[1]-ndBaseIndex,
                it->second[2]-ndBaseIndex,
                util::INVALID_IDX
            ));
        }
		return addCrack(pointList,crackLists,true);
	}

	int VDBWrapper::addCrack(
		std::vector<vdb::Vec3d>& pointList,
		std::map<unsigned int, std::vector<vdb::Vec4I> >& crackLists,
		bool haveIndexCoords, bool updateSignedCrack, std::string updateFilename
	){
        GridPtrVec updateGrids;
		std::vector<vdb::Vec3s> pointListFloat(pointList.size());
		for(unsigned int i=0; i<pointList.size(); ++i){
			//if(!haveIndexCoords){ //need to transform world coords to index space
			//	pointList[i]=resXform->worldToIndex(pointList[i]);
			//}
			//pointListFloat[i]=Vec3s(pointList[i]);
			if(!haveIndexCoords) //need to transform world coords to index space
				pointListFloat[i]=Vec3s(resXform->worldToIndex(pointList[i]));
			else
				pointListFloat[i]=Vec3s(pointList[i]);
		}
        tools::MeshToVolume<FloatGrid> voxelizer(resXform, tools::GENERATE_PRIM_INDEX_GRID); // GENERATE_PRIM_INDEX_GRID required for signedCrackGrids below
		FloatGrid::Ptr tmpGrid;
        for(map<unsigned int, vector<Vec4I> >::iterator it = crackLists.begin(); it != crackLists.end(); ++it){
			voxelizer.clear();
            voxelizer.convertToUnsignedDistanceField(pointListFloat, it->second, NBHW+1); // add 1 voxel because we dilate later

			if(updateSignedCrack) tmpGrid = voxelizer.distGridPtr()->deepCopy();

			for(FloatTree::RootNodeType::ChildOnIter cit = objectGrid->tree().beginRootChildren(); cit.test(); ++cit ){
				//printf("\n%% VDB addCrack touching leaf in root (%d, %d, %d)\n",cit.getCoord()[0],cit.getCoord()[1],cit.getCoord()[2]);
                voxelizer.distGridPtr()->getAccessor().touchLeaf( cit.getCoord());
				voxelizer.distGridPtr()->getAccessor().setValueOn(cit.getCoord());
			}
#ifdef _MSC_VER
            tools::foreach(voxelizer.distGridPtr()->beginValueAll(),crackUDFtoLS(voxelSize) ); // subtract 0.87*voxelSize and flip sign
#else // there's a strange problem when using tools::foreach with MinGW, maybe due to the tbb-lib
			for(FloatGrid::ValueAllIter vit = voxelizer.distGridPtr()->beginValueAll(); vit.test(); ++vit ){
				vit.setValue(-std::abs(*vit) + halfDiag*voxelSize);
			}
#endif
			voxelizer.distGridPtr()->signedFloodFill(); // propagate new sign to tiles with off-voxels
            if(!updateFilename.empty()){
                updateGrids.push_back(voxelizer.distGridPtr()->deepCopy());
                stringstream tmp; tmp << "crackGrid_" << it->first;
                updateGrids.back()->setName(tmp.str());
				updateGrids.back()->setGridClass(GRID_LEVEL_SET);
            }
			if(crackGrids.count(it->first)==0){ // new crack
                crackGrids[it->first]=voxelizer.distGridPtr()->deepCopy();
                crackGrids[it->first]->setGridClass(GRID_LEVEL_SET);
                stringstream tmp; tmp << "crackGrid_" << it->first;
                crackGrids[it->first]->setName(tmp.str());
            }else{ // propagated crack
                tools::csgIntersection(*crackGrids[it->first], *voxelizer.distGridPtr());
//				for(FloatTree::RootNodeType::ChildOnIter cit = objectGrid->tree().beginRootChildren(); cit.test(); ++cit ){
//					crackGrids[it->first]->getAccessor().touchLeaf( cit.getCoord());
//					crackGrids[it->first]->getAccessor().setValueOn(cit.getCoord());
//				}
			}


			if(updateSignedCrack){
				Coord c; Vec3d cv,v,n,n2; unsigned int tri;
				Int32Grid::ConstAccessor triAcc = voxelizer.indexGridPtr()->getConstAccessor();
				for(FloatGrid::ValueOnIter vit = tmpGrid->beginValueOn(); vit.test(); ++vit){
					c = vit.getCoord();
					tri = triAcc.getValue(c);
					/**/ //use the original (input) values in (possibly) world space
					// in some cases this avoids reducing normal-lengths below numeric precission
					// but we still get some errors in the sign when cracks propagate slowly (maybe due to smoothing?)
					if(!haveIndexCoords) cv=resXform->indexToWorld(c); else cv=c.asVec3d();
					v = pointList[it->second[tri][0]] - cv;
					n = ( pointList[it->second[tri][1]] - pointList[it->second[tri][0]] ).cross( pointList[it->second[tri][2]] - pointList[it->second[tri][0]] );
					if(it->second[tri][3] != util::INVALID_IDX){
						n2 = ( pointList[it->second[tri][2]] - pointList[it->second[tri][0]] ).cross( pointList[it->second[tri][3]] - pointList[it->second[tri][0]] );
						if(n2.lengthSqr()>n.lengthSqr()) n=n2;
					}
					//if(n.lengthSqr() < 1e-50) printf("\n possibly degenerate element in crack update (n.lengthSqr=%.3le ; %.1le)", n.lengthSqr(), n2.lengthSqr());
					/*/ //using the float values in index space
					v = pointListFloat[it->second[tri][0]] - c.asVec3s();//(( pointListFloat[it->second[tri][0]] + pointListFloat[it->second[tri][1]] + pointListFloat[it->second[tri][2]] )/3.0) - c.asVec3s();
					n = ( pointListFloat[it->second[tri][1]] - pointListFloat[it->second[tri][0]] ).cross( pointListFloat[it->second[tri][2]] - pointListFloat[it->second[tri][0]] );
					if(it->second[tri][3] != util::INVALID_IDX){
						n2 = ( pointList[it->second[tri][2]] - pointList[it->second[tri][0]] ).cross( pointList[it->second[tri][3]] - pointList[it->second[tri][0]] );
						if(n2.lengthSqr()>n.lengthSqr()) n=n2;
					}
					if(n.lengthSqr() < DBL_EPSILON) printf("\n possibly degenerate element in crack update (n.lengthSqr=%.3le)", n.lengthSqr());
					/**/
					if( n.lengthSqr() > DBL_EPSILON /*!=0.0*/ && v.dot(n) < 0.0 ) // check if (c - triangle midpoint ) is positive or negative projected onto the normal of the triangle
						vit.setValue(-vit.getValue());
				}
                if(!updateFilename.empty()){
                    updateGrids.push_back(tmpGrid->deepCopy());
                    stringstream tmp; tmp << "signedCrackGrid_" << it->first;
                    updateGrids.back()->setName(tmp.str());
					updateGrids.back()->setGridClass(GRID_LEVEL_SET);
                }
				if( signedCrackGrids.count(it->first)==0){ // new crack
					signedCrackGrids[it->first]=tmpGrid;
					signedCrackGrids[it->first]->setGridClass(GRID_LEVEL_SET);
					stringstream tmp; tmp << "signedCrackGrid_" << it->first;
					signedCrackGrids[it->first]->setName(tmp.str());
				}else{ // propagated crack
					vdbToolsExt::compAbsMin(*signedCrackGrids[it->first], *tmpGrid);
				}
			}
		}

        if(!updateFilename.empty()){
            io::File file(updateFilename+".vdb");
            file.write(updateGrids);
            file.close();
        }
        
		return 0;
	}

	double VDBWrapper::getValueAndGradient(const vdb::Vec3d& p, vdb::Vec3d& g, unsigned int self){
        tools::GridSampler<FloatGrid::Accessor, tools::BoxSampler> sampler(
            objectGrid->getAccessor(), objectGrid->transform()
        );
        double val=sampler.wsSample(p);
		g[0] = 2.0*( sampler.wsSample(p+0.25*voxelSize*ux) - sampler.wsSample(p-0.25*voxelSize*ux) )/voxelSize;
		g[1] = 2.0*( sampler.wsSample(p+0.25*voxelSize*uy) - sampler.wsSample(p-0.25*voxelSize*uy) )/voxelSize;
		g[2] = 2.0*( sampler.wsSample(p+0.25*voxelSize*uz) - sampler.wsSample(p-0.25*voxelSize*uz) )/voxelSize;

        //also check crack-grids and eval gradient in grid with highest value
        for(crIt it=crackGrids.begin(); it!=crackGrids.end(); ++it){
            if(it->first != self){ //do not treat self-grid here
                tools::GridSampler<FloatGrid::Accessor, tools::BoxSampler>
                    sampler(it->second->getAccessor(), *resXform);
                if(sampler.wsSample(p) > val){
                    val=sampler.wsSample(p);
                    g[0] = 2.0*( sampler.wsSample(p+0.25*voxelSize*ux) - sampler.wsSample(p-0.25*voxelSize*ux) )/voxelSize;
                    g[1] = 2.0*( sampler.wsSample(p+0.25*voxelSize*uy) - sampler.wsSample(p-0.25*voxelSize*uy) )/voxelSize;
                    g[2] = 2.0*( sampler.wsSample(p+0.25*voxelSize*uz) - sampler.wsSample(p-0.25*voxelSize*uz) )/voxelSize;
                }
            }
        }
        
        return val;
	}

    bool VDBWrapper::isInside(vdb::Vec3d p) const{
        tools::GridSampler<FloatGrid::Accessor, tools::BoxSampler> sampler(
            objectGrid->getAccessor(), objectGrid->transform()
        );
        return (sampler.wsSample(p) < -FLT_EPSILON);
    }
    bool VDBWrapper::isInside(Eigen::Vector3d p) const{
        return isInside(Vec3d(p[0],p[1],p[2]));
    }
    bool VDBWrapper::isInside(double x, double y, double z) const{
        return isInside(Vec3d(x,y,z));
    }

	int VDBWrapper::getClosestTriangleToLine(Vec3d p, Vec3d q) const{
		printf("OLD VERSION - DON'T USE (VDBWrapper::getClosestTriangleToLine)");
		// old version depending on allSurfacesGrid and closestTri grid which are no longer used
		/** 
		// search along the line (p,q) starting at p
        // for the first local minimum in allSurfacesGrid and return the
        // triangle closest to that voxel (or -1 if none has been found)
		double lenSqr = (q-p).lengthSqr(); // squared length of the line pq
		double d = 0; // distance along the line away from p
		Vec3d line = (q-p); line.normalize();
		Coord c_new, c = resXform->worldToIndexCellCentered(p); //start at p
		Int32Grid::ConstAccessor acc = closestTri->getConstAccessor();
        FloatGrid::ConstAccessor valAcc = allSurfacesGrid->getConstAccessor();
        double val = allSurfacesGrid->background();
        bool found=false;
//        printf("%%finding closest tri ...\n");
		while( (d*d)<=lenSqr && !found ){
//            printf("%% ... d=%.2le, c=(%d,%d,%d), v=%.2le\n",d,c[0],c[1],c[2],val);
			d+=halfDiag*voxelSize;
			c_new =resXform->worldToIndexCellCentered(p + d*line);
            // acc.isValueOn(c) might be obsolete here, as both grids should have the same active-states (generated by the same voxelizer)
            if(acc.isValueOn(c) && valAcc.getValue(c_new) > val){
//                printf("\n%% ... found, next value would be %.2le\n",valAcc.getValue(c_new));
//                Vec3d vox = resXform->indexToWorld(c);
//                printf("voxel=[%.3le %.3le %.3le];\n",vox[0],vox[1],vox[2]);
                found=true;
            }else{
                val=valAcc.getValue(c_new);
                c=c_new;
            }
		}
		if( found ){
			return acc.getValue(c)+ELEM_BASE_INDEX;
		}
		/**/
		return -1;
	}
	int VDBWrapper::getClosestTriangleToLine(Eigen::Vector3d p, Eigen::Vector3d q) const{
		return getClosestTriangleToLine(Vec3d(p[0],p[1],p[2]),Vec3d(q[0],q[1],q[2]));
	}
	int VDBWrapper::getClosestTriangleToLine(
		double p_x, double p_y, double p_z,
		double q_x, double q_y, double q_z) const
	{
		return getClosestTriangleToLine(Vec3d(p_x,p_y,p_z),Vec3d(q_x,q_y,q_z));
	}

	CRACK_STATE VDBWrapper::checkCrackState(Vec3d a, Vec3d old_a, int region){
		if( !isInside(a) ) return INACTIVE;

		if(!noSI) if(lineCrackSelfIntersection(a,old_a, region)) return INACTIVE;

		if( lineCracksIntersection(a,old_a, region)) return INACTIVE;

		return ACTIVE;
	}

	bool VDBWrapper::lineCrackSelfIntersection(const Vec3d& a, const Vec3d& b, int grid){
		if(crackGrids.count(grid)==0) return false; // invalid grid specified
		//check for self-intersection of a crack
		tools::GridSampler<FloatGrid::Accessor, tools::BoxSampler> sampler(
			crackGrids[grid]->getAccessor(), *resXform
		);
		double d,val,tmp,dir; // d is distance in world space
		double valAtB = sampler.wsSample(b);
		Vec3d  p;
		bool found = false;
		Vec3d line = a-b;
		double lenSqr = line.lengthSqr();
		line.normalize(); // line is now a unit vector, the length is stored in lenSqr

		d= std::max(0.0, valAtB+FLT_EPSILON); // should always satisfy tmp<(-FLT_EPSILON) condition
		while( d*d < lenSqr && !found){
			p=b+d*line; // move from (b) towards (a)
			// estimate level-set value, assume that it is locally maximal at (b), ie. (b) is already represented in the level-set
			// also assume that the level-set is based on a distance function
			tmp = valAtB - d *0.707; // d*sin(45deg) allows for non-perpendicular path
			//if( val > (tmp+halfDiag*voxelSize) )
			//	printf("\nSI found: d %.4lf, at b %.4lf, est. %.4lf, is %.4lf, diff %.2le",d,valAtB, tmp, val, tmp-val);
			val = sampler.wsSample(p);
			if( tmp< (-FLT_EPSILON) && val > (-halfDiag*voxelSize) ){ // near a crack
				//printf("%s",tmp>=0.0?"+":"-");
				dir = sampler.wsSample( p + 0.25*voxelSize*line ) 
					- sampler.wsSample( p ); // use only forward difference, as (b) might hold a local max.
				if( dir > FLT_EPSILON ){ // trigger only when moving towards a surface, as the self-surface will be part of the grid, and the crack-tip goes away from it
					found=true;
					//printf("\nSI found: d %.4lf, at b %.4lf, at p %.4lf, tmp %.4lf, dir %.4le",d,valAtB, val, tmp, dir);
				}
			}
			d+=0.5*halfDiag*voxelSize; // advance
		}
		return found;
	}

	bool VDBWrapper::lineCracksIntersection(
		const Vec3d& a, const Vec3d& b,
		int self//const vector<tools::GridSampler<FloatGrid::Accessor, tools::BoxSampler> >& samplers
	){
		double d=0.0,val,dir,tmp; // d is distance in world space
		Vec3d  p;
        int whichGrid;
		bool found = false;
		Vec3d line = a-b;
		double lenSqr = line.lengthSqr();
		line.normalize(); // line is now a unit vector, the length is stored in lenSqr

		while( d*d < lenSqr && !found){
			p=a-d*line; // start at a (d==0) then move towards b
            val=-1.0; //printf("\n***");
//            for(int i=0; i<samplers.size(); ++i){
//                tmp=samplers[i].wsSample(p);
//                if(tmp >= val){
//                    val=tmp; whichGrid=i;
//                    //printf("\n v=%.3lf (%d), bg=%.3lf",val,i,objectGrid->background());
//                }
//                if(tmp>0.0) return true;
//            }
            for(crIt it=crackGrids.begin(); it!=crackGrids.end(); ++it){
                if(it->first != self){ //do not treat self-intersections here
                    tools::GridSampler<FloatGrid::Accessor, tools::BoxSampler>
                        sampler(it->second->getAccessor(), *resXform);
                    if(sampler.wsSample(p)>0.0) return true;
                }
            }
            d+=0.5*halfDiag*voxelSize;
//			if( val > 0.0){
////				dir = samplers[whichGrid].wsSample( p + 0.25*voxelSize*line ) 
////					- samplers[whichGrid].wsSample( p - 0.25*voxelSize*line );
////				if( dir < -FLT_EPSILON ){ //try to trigger once the line goes away from the surface (has moved trough, not just towards)
//					found=true;
////				}else{
////					d+=0.4*halfDiag*voxelSize;
////				}
//			}else{
//				d+=0.5*halfDiag*voxelSize; // advance until we find a value that is positive
//			}
		}
		return found;
	}

	int VDBWrapper::writeGrids(std::string filename){
		//ToDo: check whether any of these methods throws (file open failed etc.) and handle exceptions
		GridPtrVec grids;
        for(crIt it=crackGrids.begin(); it!=crackGrids.end(); ++it){
            grids.push_back(it->second);
        }
        for(crIt it=signedCrackGrids.begin(); it!=signedCrackGrids.end(); ++it){
            grids.push_back(it->second);
        }
        grids.push_back(objectGrid);
        io::File file(filename+".vdb");
        file.write(grids);
        file.close();
		return 0;
	}
    
    VDBWrapper::VDBWrapper(double voxelSize_, double nbHWidth, bool noSI_)
        : NBHW(nbHWidth)
    {
        initialize(); // initialize OpenVDB
        resXform = math::Transform::createLinearTransform(voxelSize_);
		voxelSize=voxelSize_;
		noSI=noSI_;
		nearTriGridInitialized=false;
    }

    VDBWrapper::VDBWrapper(const VDBWrapper& ori) : NBHW(ori.NBHW) {
        resXform = ori.resXform->copy();
		voxelSize = ori.voxelSize;
        objectGrid = ori.objectGrid->deepCopy();
		noSI=ori.noSI;
        //closestTri = ori.closestTri->deepCopy();
		//allSurfacesGrid = ori.allSurfacesGrid->deepCopy();
		nearTriGridInitialized=ori.nearTriGridInitialized;
		if(nearTriGridInitialized) nearTriGrid=ori.nearTriGrid->deepCopy();
		crackGrids.clear();
        crackGrids.insert(ori.crackGrids.begin(), ori.crackGrids.end());
        signedCrackGrids.clear();
        signedCrackGrids.insert(ori.signedCrackGrids.begin(), ori.signedCrackGrids.end());
    }

    VDBWrapper::~VDBWrapper() { // using smart-pointers (boost) ;-)
    }

	//void dumpGrid(FloatGrid& g, string file){
	//template<class T> void dumpGrid(Grid<T>& g, string file){
	//	ofstream out(file);
	//	if(out.is_open()){
	//		for(typename Grid<T>::ValueOnIter it=g.beginValueOn(); it.test(); ++it){
	//			out << "(" << it.getCoord()[0] << ", " << it.getCoord()[1] << ", " << it.getCoord()[2] << ") = " << it.getValue() << endl;
	//		}
	//		out.close();
	//	}
	//}

	// NEW FOR FractureRB
	vdb::FloatGrid::Ptr VDBWrapper::maskCracks(vdb::FloatGrid::Ptr surface){
		vdb::FloatGrid::Ptr result=crackGrids.begin()->second->deepCopy();
		vdb::FloatGrid::Ptr surfCpy=surface->deepCopy();
		result->setName("masked_cracks");
		// first copy all the cracks together into one grid
        crIt it=crackGrids.begin(); ++it;
		for(; it!=crackGrids.end(); ++it){
			vdb::tools::csgIntersection(*result, *(it->second)->deepCopy());
        }
		for(FloatGrid::ValueAllIter vit = result->beginValueAll(); vit.test(); ++vit ){
			vit.setValue(-*vit); // flip the sign
		}
		//for(FloatGrid::ValueAllIter vit = surfCpy->beginValueAll(); vit.test(); ++vit ){
		//	vit.setValue(*vit+voxelSize); // erode the surface for masking
		//}
		// now intersect the surface with the grid containing all cracks
		vdb::tools::csgIntersection(*result, *surfCpy);
//		result->prune();
//		result->treePtr()->pruneInactive();
		// needs #include <openvdb/tools/LevelSetRebuild.h>
		result=vdb::tools::levelSetRebuild(*result,0.0,1);

		//// write file for debugging ...
		//GridPtrVec grids;
  //      grids.push_back(surface);
		//grids.push_back(surfCpy);
		//grids.push_back(result);
  //      io::File file("debug_masking.vdb");
  //      file.write(grids);
  //      file.close();

		return result;
	}

//******************************************************************************
//******************************************************************************
//******************************************************************************
// old stuff ...

	CRACK_STATE VDBWrapper::checkCrackState(Vec3d a, Vec3d b, Vec3d old_a, Vec3d old_b, const CRACK_STATE oldState, int region){
        printf("\nOLD VERSION - DON'T USE (VDBWrapper::checkCrackState(Vec3d,Vec3d,Vec3d,Vec3d,const CRACK_STATE,int))\n");
        //before doing complicated stuff make sure at least one of (a,b) is inside the object
		bool found_a = !isInside(a);
		bool found_b = !isInside(b);
		if( found_a && found_b ) return INACTIVE;

        vector<tools::GridSampler<FloatGrid::Accessor, tools::BoxSampler> > samplers; // samplers for all crack-grids except the one containing 'region'
		int selfGrid=-1;
        for(crIt it=crackGrids.begin(); it!=crackGrids.end(); ++it){ // build list of accessors and samplers of grids that will be tested for intersections
            //printf("\n%d --> %p",it->first, it->second.get());
            if(it->first != region){ //ignore self-intersections
                samplers.push_back( tools::GridSampler<FloatGrid::Accessor, tools::BoxSampler>( 
                        it->second->getAccessor(), *resXform
                ));
            }else{
				selfGrid=it->first;
			}
        }
//		found_a = found_a || lineCrackSelfIntersection(a,old_a, selfGrid);
//		found_b = found_b || lineCrackSelfIntersection(b,old_b, selfGrid);
//        if( lineCrackSelfIntersection(a,old_a, selfGrid) || // stricter version for self-intersections
//            lineCrackSelfIntersection(b,old_b, selfGrid) ) return INACTIVE;
        
        // these methods have a new interface now ...
//		found_a = found_a || lineCracksIntersection(a,old_a, samplers);
//		found_b = found_b || lineCracksIntersection(b,old_b, samplers);
		if( (found_a && found_b) ||
			(found_a && oldState==UPDATE_A) ||
			(found_b && oldState==UPDATE_B) // one node has intersected earlier, the other now --> arrest
		) return INACTIVE;
        if(found_a) return ACTIVE_B;
        if(found_b) return ACTIVE_A;
        // no change, just switch back from UPDATE_* to ACTIVE_*
        if(oldState==UPDATE)   return ACTIVE;
        if(oldState==UPDATE_A) return ACTIVE_A;
        if(oldState==UPDATE_B) return ACTIVE_B;
        return oldState;
	}

}
