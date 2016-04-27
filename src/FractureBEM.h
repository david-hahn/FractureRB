/* 
 * File:   FractureBEM.h
 * Author: David
 *
 * Created on 15. Jänner 2014, 13:49
 */

#ifndef FRACTUREBEM_H
#define	FRACTUREBEM_H

#include "types.h"
#include <boost/scoped_ptr.hpp>
#include <string>

namespace FractureSim{
	// fwd decl
	class HyENAWrapper;
	class PostProcessor;
	class VDBWrapper;
	class FractureModel;
	class MaterialModel;
    class SubsampledCrackTip;
    
    class FractureBEM {
    public:
		FractureBEM(double bemMeshSize=1.0, FRACTURE_METHOD m=FULL_BEM);
        //FractureBEM(const FractureBEM& ori);
		virtual ~FractureBEM();
        
        // byRef getters ...
        inline node_map& getNodes(){ return nodes; }
        inline elem_map& getElems(){ return elems; }
		inline id_map& getRegions(){ return regions; }
        inline id_set& getCracks(){ return cracks; }
        inline elem_map& getCrackTips(){ return crackTips; }
        inline elem_map& getCrackTipParents(){ return crackTipParents; }
        inline state_map& getCrackTipStates(){ return crackTipStates; }
		// NEW FOR FractureRB:
		inline id_set& getFracturedNodes(){ return fracturedNodes; }
		inline id_set& getFracturedElems(){ return fracturedElems; }
		inline void recomputeCrackAreas(){ crackAreas.clear(); updateCrackAreas(elems); }
		inline FRACTURE_METHOD getMethod(){ return method; }
		inline void setMethodEstimated(){ method=INITIAL_BEM; }
		
		unsigned int getActiveCount();
        
        /* Initialize the BEM solver using a homogeneous (E,nu) material
         * Call after filling initial geometry (i.e. nodes, elems, regions, crackTips, crackTipParents)
         * Make sure crackTips (line elements!) and crackTipParents (containing tris)
         * hold data consistent with bndCnds.
         * Boundary conditions given as specified in HyENAWrapper::buildBoundaryCnd
         */
        int initBEM(MaterialModel& material, bool is2D, const std::vector<std::string>& bndCnds);
		int initBEM(
            double youngsModulus, double poissonsRatio, bool is2D,
            const std::vector<std::string>& bndCnds
        );
        
        /* ToDo: write comment here
         */
        int initVDB(double voxelSize=1.0, bool noSI=false, double nbHWidth=3.0);
        
        /* ToDo: write comment here
		 * ToDo: implement setters such that the caller can specify a fracture model when propagating cracks
		 */
		int initFractureModel(double strength=1.0, double toughness=1.0, double density=1.0, double compressive=3.0);
		int initFractureModel(MaterialModel& material, bool substepOutput=false, double samplingDensity=0.8);

        /* Solve the BEM deformation based on the current matrices and boundary conditions
		 * also updates the nodal stress intensity factors based on the BEM result
         */
		int computeBEM();

		/* Update the BEM matrices to include recently added cracks, then solve
		 * use this method after adding cracks using the startCrack(...) method
		 * calls computeBEM()
		 */
		int computeAddedCracks(); // call this after manually adding crack geometry

        /* Allow cracks to propagate for one time-step, based on the material model
		 * and the most recently computed nodal stress intensities - see computeBEM()
         */
		int propagateCracks();

		/* Create new cracks based on the material model and principal surface stresses
		 * the number of additional cracks can be limited using the maxSeed parameter
		 * once seeding is complete, propagateCracks() is called
         */
		int seedCracksAndPropagate(int maxSeed=-1);
        
        /* All crack tips whose parent has a region-ID listed in regions
         * are set to inactive state - won't propagate anymore.
         */
        void setCracksInactive(const id_set& regions);

		/* Add a completely new edge crack, starting from the given element
		 * in the plane given by its normal vector
		 * for fracture elements, we can specify to start on the "negative" side of the element
		 */
		int startCrack(unsigned int elem, Eigen::Vector3d normal, bool negSide=false);
		int startCrack(Eigen::Vector3d point, Eigen::Vector3d planeNormal, Eigen::Vector3d tipNormal);

		//ToDo: write comment here
		int remesh(int nTris, double adaptive=0.1, double offsetVoxels=0.0);
		int writeMesh(std::string filename, double visualQuality=2.0, // visualQuality >1 no output, in [0,1] use VDB adaptive meshing, <0 use quadric decimation
            bool visDisplace=true, bool visCOD=true, bool visClose=false, bool visOBJ=false
        );
		int writeVDB(std::string filename, bool updateSeeds=true); //NEW FOR FractureRB: option to suppress updates to the last written file when seeding cracks
		int writeCrackTip(std::string filename);
		int innerEvalTest(std::string filename);
		
		void setCPVersion(unsigned int v); // use 0, 1 or 2 to select the crack propagation algorithm (see FractureModel.h for details)

		/* NEW FOR FractureRB:
		 * find the closest triangle using the near triangle grid
		 * returns -1 if no triangle is found in the vicinity of p
		 */
		unsigned int findClosestSurfaceTri(const Eigen::Vector3d& p);
		/* NEW FOR FractureRB:
		 * allow access to the BEM results
		 */
		vector_type& getDisplacements();
		vector_type& getTractions();

		/* NEW FOR FractureRB:
		 * read the current timestep from the fracture model
		 */
		double getTimeStep();

		/* NEW FOR FractureRB:
		 * update the boundary data (but not the type of the boundary condition)
		 */
		int updateBoundaryData(vect3d_map& newBcData);
		
		/* NEW FOR FractureRB:
		 * allow access to the implicit surface and hi-res crack-front
		 */
		VDBWrapper& getLevelSet();
		SubsampledCrackTip& getHiResCrackTip();

		/* NEW FOR FractureRB:
		 * report the target edge length in the BEM mesh
		 */
		inline double getTargetMeshSize(){ return crackMeshSize; }
		
		/* print all sorts of stuff to stdout
		 * like the current mesh, data arrays, etc.
		 * warning: output will be modified as needed
		 */
		void dumpDebugData();
	protected:
        // explicit geometry representation (BEM mesh)
        node_map nodes; // map node-ID --> coordinates
        elem_map elems; // map element-ID --> node-IDs (triangles)
        id_map regions; // map element-ID --> region-ID
        // fracture handling (at BEM mesh resolution)
        id_set cracks; //set of region-IDs which are cracks (as opposed to boundaries of the object)
        elem_map crackTips; // line elements along the crack tip
        elem_map crackTipParents; // map a crack tip (line) element to its containing triangle(s)
        state_map crackTipStates; // map a crack tip (line) element to it's state
		vect3d_map nodeSIFs; // stress intensity factors at nodes on the crack-tip
		vect3d_map crackTipFaceNormals; // region normals at crack tip nodes
		vect3d_map crackTipTangents; // crack-tip tangents at crack tip nodes
		value_map crackAreas; // map region-id of a crack to it's mesh surface area (only used for SIF estimates, not for BEM/DCT based SIF computations)
		double crackMeshSize;
		FRACTURE_METHOD method;

		id_set fracturedNodes; // <- for node based seeding: avoid starting a crack from the any one node more than once
		id_set fracturedElems; // <- for elem based seeding: avoid starting a crack from the any one element more than once

        // object pointers
        boost::scoped_ptr<HyENAWrapper>       bemSolver; // manages HyENA BEM solving
        boost::scoped_ptr<PostProcessor>      postPro;   // takes care of crack propagation and VTK output
        boost::scoped_ptr<VDBWrapper>         levelSet;  // maintains VDB grids representing the geometry implicitly
		boost::scoped_ptr<FractureModel>      fractureModel; // implements the fracture criterion
		boost::scoped_ptr<MaterialModel>      materialModel; // enables variable fracture toughness
		boost::scoped_ptr<SubsampledCrackTip> fineCrackTip; // detailed sampling of the crack tip at level-set resolution
        // state flags
		bool bemInitialized, vdbInitialized, fractureInitialized, bemFirst;
		unsigned int firstCrackElement; // we require that all elements with IDs < firstCrackElement are regular surface elements and all IDs >= firstCrackElement are crack elements

		// helper functions
		bool checkMesh();
		bool checkDistanceNd(unsigned int node, const id_set& nodes); // return true iff the node is further than mesh-resolution from all nodes in the set
		bool checkDistanceEl(unsigned int elem, const id_set& nodes); // return true iff all nodes of the given element are further than mesh-resolution from all given nodes
        
		// add the areas of the elements in newElems to the corresponding crackAreas
		// requires that the coordinates and regions of the elements are already stored
		void updateCrackAreas(const elem_map& newElems);

        //// using template to avoid inclusion of openvdb header
        //template<typename Vec3d>
        //int intersectElemLine(Vec3d a, Vec3d b, unsigned int tri, Vec3d& p);
		std::string lastVDBfile; //store last written vdb output file, overwrite this after crack seeding, if new cracks were added.
    };
}

#endif	/* FRACTUREBEM_H */
