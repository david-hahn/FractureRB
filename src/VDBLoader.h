/* 
 * File:   VDBLoader.h
 * Author: David
 *
 * Created on 16. Juli 2014
 */

#ifndef VDBLOADER_H
#define	VDBLOADER_H

#include "VDBWrapper.h"

#include <string>

namespace FractureSim{
    class VDBLoader : protected VDBWrapper {
    public:
        int loadFull(std::string file, id_set& cracks);
        int loadUpdate(std::string file, id_set& cracks);
		int intersectWithCrackGrids(vdb::FloatGrid::Ptr grid); // NEW FOR FractureRB

        inline int writeVisualMesh(
			node_map& nodes, elem_map& elems, id_map& regions, id_set& cracks,
			const vector_type& u, const vector_type& u_c, std::string filename,
            double adaptiveVDB, double adaptiveQuadric, bool visDisplace=true, bool visCOD=true,
            bool visClose=false, bool visOBJ=false, double segment=-1.0, bool writePerSegment=false, double postDecimation=-1.0
		){
            return VDBWrapper::writeVisualMesh(
                nodes, elems, regions, cracks, u, u_c, filename, adaptiveVDB, adaptiveQuadric,
				visDisplace, visCOD, visClose, visOBJ, segment, writePerSegment, postDecimation
            );
        }
        inline int writeGrids(std::string filename){
            return VDBWrapper::writeGrids(filename);
        }

        inline int addToNearTriGrid(node_map& nodes, elem_map& elems, double res=1.0){
            return VDBWrapper::addToNearTriGrid(nodes, elems, res);
        }

		// NEW FOR FractureRB
		inline std::vector<vdb::FloatGrid::Ptr> getSegments(double handleThreshold=0.0, bool useTiles=true){
			return VDBWrapper::getSegments(handleThreshold, useTiles);
		}
		// NEW FOR FractureRB
		inline vdb::FloatGrid::Ptr getObjectGrid(){
			return VDBWrapper::getObjectGrid();
		}
		// NEW FOR FractureRB
		inline void setObjectGrid(vdb::FloatGrid::Ptr grid){
			VDBWrapper::setObjectGrid(grid);
		}
		// NEW FOR FractureRB
		inline double getBandwidth(){
			return VDBWrapper::getBandwidth();
		}
		// NEW FOR FractureRB
		inline double getVoxelSize(){
			return VDBWrapper::getVoxelSize();
		}
    };
}
#endif	/* VDBLOADER_H */

