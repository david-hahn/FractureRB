/* 
 * File:   VDBWrapper.h
 * Author: David
 *
 * Created on 14. Jänner 2014, 13:47
 */

#ifndef VDBWRAPPER_H
#define	VDBWRAPPER_H

#include "types.h"

#include <openvdb/openvdb.h>
#include <openvdb/tools/Composite.h>

#include <Eigen/Dense>
#include <string>

namespace vdb = openvdb::v2_2_0;

typedef std::map<unsigned int, vdb::FloatGrid::Ptr>::iterator crIt;

namespace vdbToolsExt{ // for signedCrackGrids we need to merge two grids taking the value closest to zero at each voxel
	using namespace vdb;
	using namespace vdb::tools;
	// this is a modified version of <openvdb/tools/Composite.h> compMin(...) method
	template<typename GridOrTreeT>
	OPENVDB_STATIC_SPECIALIZATION inline void
	compAbsMin(GridOrTreeT& aTree, GridOrTreeT& bTree)
	{
		typedef TreeAdapter<GridOrTreeT>    Adapter;
		typedef typename Adapter::TreeType  TreeT;
		typedef typename TreeT::ValueType   ValueT;
		struct Local {
			static inline void op(CombineArgs<ValueT>& args) {
				args.setResult(
					( std::abs(args.a()) < std::abs(args.b()) ) ? args.a() : args.b()
				);
			}
		};
		Adapter::tree(aTree).combineExtended(Adapter::tree(bTree), Local::op, /*prune=*/false);
	}
}

namespace myVCG{
	class VcgMesh; // fwd decl
}

namespace FractureSim{
	const double halfDiag = 0.867; // roughly the half space diagonal of the unit cube ~~ sqrt(3)/2



    class VDBWrapper {
    public:
		const static vdb::Vec3d ux,uy,uz; // unit vectors in x,y,z axis
        VDBWrapper(double voxelSize=1.0, double nbHWidth=3.0, bool noSI=false);
        VDBWrapper(const VDBWrapper& ori);
        virtual ~VDBWrapper();
        
        /* Inizialize the wrapped grids from a boundary-mesh
         * elems maps an    element-ID to 3 node-IDs (triangle)
         * nodes maps a     node-ID    to 3 coordinates (xyz)
         * regions maps an element-ID to a region-ID
         * cracks is a set of region-IDs indicating which regions are cracks
         * all regions not in cracks are treated as boundaries of the object
         */
        int init(
            node_map& nodes, elem_map& elems,
            id_map& regions, id_set& cracks
        );

		/* Update the level-set representation of cracks
		 * to include the given geometry.
		 */
		int addCrack(
			node_map& nodes, elem_map& newElems,
			id_map& newRegions, unsigned int ndBaseIndex = NODE_BASE_INDEX
		);
		int addCrack(
			std::vector<vdb::Vec3d>& pointList,
			std::map<unsigned int, std::vector<vdb::Vec4I> >& crackLists,
			bool haveIndexCoords=false, // pointList may contain world or index space coordinates
			bool updateSignedCrack=true, // also update the signed crack grids from the normals of the new elements
            std::string updateFilename=""
		);

        /* Test whether a point (in world coords (as opposed to index coords))
         * is inside the object.
         */
        bool isInside(vdb::Vec3d p) const;
        bool isInside(Eigen::Vector3d p) const;
        bool isInside(double x, double y, double z) const;

		/* Search along a line for the triangle closest to it
		 * starting at p
		 */
		int getClosestTriangleToLine(vdb::Vec3d p, vdb::Vec3d q) const;
		int getClosestTriangleToLine(Eigen::Vector3d p, Eigen::Vector3d q) const;
		int getClosestTriangleToLine(double p_x, double p_y, double p_z, double q_x, double q_y, double q_z) const;

		/* Check whether a crack segment (a,b) which has propagated from (old_a,old_b)
		 * should be allowed to propagate any further
		 * returns ACTIVE if ok, INACTIVE if not, ACTIVE_A if node a is valid but b is not
		 * or ACTIVE_B if node b is valid but a is not.
         * If region is specified, crack-grids having that region-ID are ignored
		 */
		CRACK_STATE checkCrackState(vdb::Vec3d a, vdb::Vec3d b, vdb::Vec3d old_a, vdb::Vec3d old_b, const CRACK_STATE oldState, int region=-1);
        /* Simpler version: check only along one line, rather than both sides of a quad
         * returns ACTIVE if the line has not intersected a surface, INACTIVE otherwise
         */
        CRACK_STATE checkCrackState(vdb::Vec3d a, vdb::Vec3d old_a, int region=-1);

		int writeGrids(std::string filename);
		inline double getVoxelSize(){ return voxelSize; }
		inline double getBandwidth(){ return NBHW; }

		double getValueAndGradient(const vdb::Vec3d& p, vdb::Vec3d& g, unsigned int self);
		inline double getValueAndGradient(const Eigen::Vector3d& p, Eigen::Vector3d& g, unsigned int self){
			vdb::Vec3d vg;
			double v=getValueAndGradient(vdb::Vec3d(p[0],p[1],p[2]), vg, self);
			g[0]=vg[0]; g[1]=vg[1]; g[2]=vg[2];
			return v;
		}

		// The following public methods are implemented in VDBWrapper_mesh.cpp ...

		/* Create a tri-mesh from the level-set (object grid)
		 * and decimate to nTris, default is no change (current voxel-size) when nTris<0
		 * adaptivity parameter is between 0 (no adaptivity) and 1 (max adaptivity)
         * is passed on to the VDB-meshing routine
		 * the final mesh can be offset along the vertex normals by offsetVoxels to make sure
		 * that most of the detailed surface is inside the coarse mesh.
		 */
		int mesh(node_map& nodes, elem_map& elems, int nTris=-1, double adaptive=0.1, double offsetVoxels=0.0);

		/* Add elements to the near triangle grid
		 * the resolution of the grid (res) must be specified at the first call
		 * to initialize the grid, later on it will be ignored
		 */
		int addToNearTriGrid(node_map& nodes, elem_map& elems, double res=1.0);

		/* NEW FOR FractureRB:
		 * find the closest triangle using the near triangle grid
		 * the specified nodes and elems should match the union of all
		 * for which addToNearTriGrid has been called previously
		 * returns -1 if no triangle is found in the vicinity of p
		 * triangles that map to a region in cracks will be ignored
		 */
		unsigned int findClosestSurfaceTri(const Eigen::Vector3d& p,
			node_map& nodes, elem_map& elems, id_map& regions, id_set& cracks
		);

		/* Output a high-detail tri-mesh mapped to world-space by interpolating
		 * displacements in the coarse BEM mesh
		 * u   is a vector containing surface displacements and crack-opening displacements
		 * u_c is a vector containing crack-base displacements
         * if adaptive>0 it is passed on to OpenVDB's volumeToMesh function
         * if adaptive<0 a quadric decimation is run on the full-res mesh until
         * the quadric error reaches (-adaptive)
		 * if segment is >=0.0 segmentation is performed with the given value as threshold for handle-removal
		 * it is recommended that segment<=1.0 (in voxel units)
		 */
		int writeVisualMesh(
			node_map& nodes, elem_map& elems, id_map& regions, id_set& cracks,
			const vector_type& u, const vector_type& u_c, std::string filename,
            double adaptiveVDB=0.01, double adaptiveQuadric=-1.0, bool visDisplace=true, bool visCOD=true,
            bool visClose=false, bool visOBJ=false, double segment=-1.0, bool writePerSegment=false, double postDecimation=-1.0
		);
//        int writeVisualMesh2(
//            std::string filename, HyENAWrapper* bem ,double adaptive=0.0
//		);
		
		/* NEW FOR FractureRB:
		 * run the segmentation and return raw results
		 */
		std::vector<vdb::FloatGrid::Ptr> getSegments(double handleThreshold=0.0, bool useTiles=false) const;
		
		/* NEW FOR FractureRB:
		 * allow access to the objectGrid (stores the surface of the object without any fractures)
		 */
		vdb::FloatGrid::Ptr getObjectGrid() {return objectGrid;}
		void setObjectGrid(vdb::FloatGrid::Ptr grid) { objectGrid = grid; }
		vdb::math::Transform::Ptr getTransform() {return resXform; }

		/* NEW FOR FractureRB:
		 * report whether crack self-intersections are ignored (returns true) or handled as normal surface intersections (returns false)
		 */
		double ignoreCrackSelfIntersection() {return noSI;}
		
		/* NEW FOR FractureRB:
		 * deep-copy this object's crack grids to the given maps
		 * newCrackIDs allows to assign new region IDs for the target maps
		 * input maps will be cleared!
		 */
		void copyCrackGrids(
			std::map<unsigned int, vdb::FloatGrid::Ptr>& crackGrids_,
			std::map<unsigned int, vdb::FloatGrid::Ptr>& signedCrackGrids_,
			id_map newCrackIDs , vdb::FloatGrid::Ptr surface
		){
			crackGrids_.clear();
			signedCrackGrids_.clear();
			try{
				//for(crIt it=crackGrids.begin(); it!=crackGrids.end(); ++it){
				//	crackGrids_[ ((newCrackIDs.count(it->first)>0)?(newCrackIDs[it->first]):(it->first)) ]
				//		=it->second->deepCopy();
				//}
				//for(crIt it=signedCrackGrids.begin(); it!=signedCrackGrids.end(); ++it){
				//	signedCrackGrids_[ ((newCrackIDs.count(it->first)>0)?(newCrackIDs[it->first]):(it->first)) ]
				//		=it->second->deepCopy();
				//}
				unsigned int id; bool test;
				vdb::FloatGrid::Accessor surf_acc = surface->getAccessor();
				for(crIt it=crackGrids.begin(); it!=crackGrids.end(); ++it){
					id = ((newCrackIDs.count(it->first)>0)?(newCrackIDs[it->first]):(it->first));
					test = false;
					// test if the current crack grid has any active voxel inside of the given surface
					for(vdb::FloatGrid::ValueOnIter vit = it->second->beginValueOn(); vit.test() && !test; ++vit ){
						if( surf_acc.getValue( vit.getCoord() ) <= 0.0 ) test=true;
					}

					if( test ){
						crackGrids_[id]       = it->second->deepCopy();
						signedCrackGrids_[id] = signedCrackGrids[it->first]->deepCopy();
					}//else{ printf("\n%% ignoring crack %d->%d no inside voxels",it->first, id); }
				}
			}catch(std::exception& ex){
				printf("\n!!! %s in copyCrackGrids()", ex.what());
				crackGrids_.clear();
				signedCrackGrids_.clear();
				return;
			}catch(...){ printf("\n!!! error in copyCrackGrids()");}
			//printf(" (copied %d cracks) ", crackGrids_.size());
		}
		
		/* NEW FOR FractureRB:
		 * deep-copy this object's crack grids to another VDBWrapper object
		 * the crack grids of the input object will be cleared!
		 * newCrackIDs allows to assign new region IDs for the target maps
		 */
		void copyCrackGrids( VDBWrapper& other , id_map newCrackIDs , vdb::FloatGrid::Ptr surface ){
			copyCrackGrids( other.crackGrids, other.signedCrackGrids, newCrackIDs, surface );
		}

		/* NEW FOR FractureRB:
		 * write all crack IDs stored in this object into the output param. set
		 * output set will be cleared first
		 */
		void getCrackIDs( std::set<unsigned int>& crackRegions ){
			crackRegions.clear();
			for(crIt it=crackGrids.begin(); it!=crackGrids.end(); ++it)
				crackRegions.insert( it->first );
		}
		
		/* NEW FOR FractureRB:
		 * build a new grid which contains all cracks stored in this object
		 * and intersect them with the given surface
		 */
		vdb::FloatGrid::Ptr maskCracks(vdb::FloatGrid::Ptr surface);
		
		/* NEW FOR FractureRB:
		 * get a set of nearby triangles for a given point p
		 * output parameter (tris) will be cleared before writing
		 * returns a bbox specifying the grid box from which the
		 * triangle set has been read
		 */
		vdb::BBoxd getNearTris(const vdb::Vec3d& p, id_set& tris);

    protected:
		bool noSI;
		// objectGrid is a (proper) (narrow-band) level-set used for in-/outside test
        vdb::FloatGrid::Ptr objectGrid;

		template<typename ValueType> class VdbValueAdapter{
		public:
			ValueType& get() const {return data;}
			void set(ValueType& value) const {data=value;}

			VdbValueAdapter(int zero=0){}
			VdbValueAdapter(const VdbValueAdapter& ori){data=ori.data;}

			operator bool() const{return false;}
			VdbValueAdapter& operator=(const VdbValueAdapter& ori){data=ori.data; return *this;}
			bool operator==(const VdbValueAdapter& rhs) const {return data ==rhs.data; }
			bool operator!=(const VdbValueAdapter& rhs) const {return data !=rhs.data; }
			bool operator< (const VdbValueAdapter& rhs) const {return data < rhs.data; }
			bool operator> (const VdbValueAdapter& rhs) const {return data > rhs.data; }
			template<class T> const VdbValueAdapter& operator+(const T&) const { return *this; }
			template<class T> const VdbValueAdapter& operator-(const T&) const { return *this; }
		private:
			mutable ValueType data;
		};
        typedef VdbValueAdapter<id_set> IdSetAdapter;
		typedef vdb::Grid<vdb::tree::Tree4<IdSetAdapter,5,4,3>::Type > IdSetGrid;
		IdSetGrid::Ptr nearTriGrid;
        vdb::math::Transform::Ptr nearTriXform;
		bool nearTriGridInitialized;
        
        // each crackGrid is a (narrow-band) unsigned distance function to a
        // bounded 2-manifold (single-sided) crack surface
        std::map<unsigned int, vdb::FloatGrid::Ptr> crackGrids;

		std::map<unsigned int, vdb::FloatGrid::Ptr> signedCrackGrids; // as above, but keep information about which is the positive/negative side of the crack, do not change dist.fcn. values after voxelization
		//ToDo: if signedCrackGrids work as intended (primarily for post-processing atm.) re-write all functions using crackGrids to work on signedCrackGrids instead and remove crackGrids.

        // transform for all grids (voxel resolution)
        vdb::math::Transform::Ptr resXform;
		double voxelSize;

        // narrow band half width
        const double NBHW; // usually 3 [voxels]

		// test for intersection of the line (a,b)
		bool lineCracksIntersection(
			const vdb::Vec3d& a, const vdb::Vec3d& b,
			int self
		);
		// test for intersection of the line (a,b) with the crack-grid specified in 'grid',
		// treat the grid as the representation of the crack which is propagated (self-intersection test)
		bool lineCrackSelfIntersection(const vdb::Vec3d& a, const vdb::Vec3d& b, int grid);

		char interpDisplacement(
			Eigen::Vector3d& u_p, const Eigen::Vector3d p,const Eigen::Vector3d n, const vector_type& u,
			unsigned int tri, double s, double t, node_map& nodes, elem_map& elems, int crackRegion=-1, double bgval=1.0
		);

		int outputSubMesh( int subMeshId, myVCG::VcgMesh& mesh,
			node_map& nodes, elem_map& elems, id_map& regions, id_set& cracks,
			const vector_type& u, const vector_type& u_c, std::string filename,
            double adaptive, bool visDisplace, bool visCOD, bool visClose, bool visOBJ, double postDecimation
		);
    };
}
#endif	/* VDBWRAPPER_H */

