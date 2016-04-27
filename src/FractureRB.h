/* 
 * File:   FractureRB.h
 * Author: David
 *
 * Created on 08. Jan 2015
 */

#ifndef FRACTURERB_H
#define	FRACTURERB_H

#include "types.h"
#include "ColliderData.h"

#include <string>

//fwd decl
class btRigidBody;
class btCollisionShape;
class btStridingMeshInterface;
class btConvexHullShape;
class btBoxShape;

namespace FractureSim{
	
	//fwd decl
	class FractureBEM;
	class MaterialModel;
	class VDBWrapper;
	
    class FractureRB { // this class encapsulates a Bullet rigid body and a FractureBEM object
    public:
		FractureRB();
		FractureRB(btRigidBody* rb, bool ownRB=true); //if ownRB==true, takes ownership of rigid body - will be deleted by destructor
		virtual ~FractureRB();
		
		int initFractureSim(
			double crackMeshSize,
			double voxelSize,
			double young,
			double poisson,
			double density,
			double strength,
			double toughness,
			double compFactor = 1.0,
			int    remesh = 0,
			bool   ignoreCrackSelfIntersect = false,
			std::string matSpec = "default",
			int maxCracks = -1,
			FRACTURE_METHOD method = FULL_BEM,
			double samplingDensity=1.0
		);
		
		inline btRigidBody* getRB(){return rb;}
		inline MaterialModel* getMaterial(){return mat;}
		
		/* computes E/(1-nu^2) where E is Young's modulus and nu is Poisson's ratio
		 */
		double getEovNuSq();

		/* For a given point p, find the closest surface element in the BEM mesh
		 * writes the element ID into the elem parameter
		 * returns the local mean curvature
		 */
		double getCurvatureAndElementID(const Eigen::Vector3d& p, unsigned int& elem);
		
		/* Clear all previously added contact points
		 */
		void clearContacts();

		/* Add a contact point which will contribute a boundary condition for the fracture simulation
		 * p is the point of impact in this rigid body's local coordinate system (alternatively a triangle-ID can be specified)
		 * d is the direction of the impact (also in local coords), should point towards the inside for collisions
		 * the simulation will run for as many time-steps as required by the contact with the longest duration
		 */
		int  addContact(Eigen::Vector3d& p, Eigen::Vector3d& d, double impulse, double duration=1.0);
		int  addContact(unsigned int tri, Eigen::Vector3d& d, double impulse, double duration=1.0);

		/* Check if the current set of contact points is sufficient to run a fracture simulation
		 * returns 0 if the longest contact duration is too short for the time-step size
		 * otherwise returns the total magnitude of deformational forces
		 * the caller should then decide whether the force magnitude warrants running the simulation
		 */
		double getTotalContactForce();

		/* Run the fracture simulation using the previously added contacts as boundary conditions
		 * rbTimeCode is used to distinguish output files from different calls to this method
		 * returns the number of new BEM elements (0 ... no fractures) or negative values for errors
		 */
		int  runFractureSim(double maxTime, int rbTimeCode=0);

		/* Segment the implicit surface to find separate fragments
		 * create breakable rigid bodies for large fragments
		 * and standard rigid bodies for small ones
		 * the output parameters largeFragments and smallFragments will be cleared
		 * if one fragment is almost as big as this object,
		 * update the implicit surface of this object with the huge fragment
		 * returns a positive value if the implicit surface was updated (hence this object should be kept in the RB world)
		 * the value indicates the number of the mesh-file to be drawn from now on for this rigid body
		 * returns 0 if this object should be removed from the RB world and negative values for errors
		 */
		int splitFragments(std::set<FractureRB*>& largeFragments, std::set<btRigidBody*>& smallFragments);

		/* sets the prefix for all output files
		 */
		inline void setOutputDir(std::string newOutDir){outDir = newOutDir;}
		inline std::string getOutputDir(){ return outDir; }
		
		/* notify the BulletWrapper that this objects has a new collision shape
		 * that should be updated, the wrapper must first remove this objects RB
		 * from the dynamics world, then call doCollisionUpdate() and finally
		 * add this object's RB back into the dynamics world
		 */
		inline bool pendingCollisionUpdate(){ return collUpdate.shape!=NULL; }
		/* if pendingCollisionUpdate()==true this method actually changes the
		 * collision shape of the encapsulated rigid body, which MUST NOT be
		 * part of a Bullet dynamics world at the time.
		 */
		void doCollisionUpdate();

		/* set the mass and inertia of the rigid body
		 * based on the volume of the implicit surface
		 * and the density of the material
		 * returns the implicit surface volume (in world space)
		 */
		double updateMass();
		
		void setEstSIFsThreshold( int useEstSIFs_ ){ useEstSIFs=useEstSIFs_; }

		inline const Eigen::Vector3d& getCOMshift(){ return vdbCOM; }
		
	protected:
		bool haveMesh, haveParams, haveBEM; // basic state checking
		int maxCracks; // seeding cutoff for new fractures (-1 is unlimited)
		int useEstSIFs; // switch to estimated SIFs if the number of elements exceeds this value at the start of a time-step (-1 == never switch)
		unsigned int remeshTarget; // number of surface elements for BEM mesh generation
		double samplingDensity; // density of crack-front sampling (1.0 means roughly 1 sample per voxel)
		std::string outDir; // directory and file prefix for output files
		btRigidBody* rb; // Bullet rigid body
		bool ownRB; // store if this object should delete the RB in the destructor
		FractureBEM* fractSim; // BEM fracture sim handler
		MaterialModel* mat; // material descriptor
		vect3d_map elemCtr; // store centroid of surface elements
		value_map elemArea, nodeCurv; // store area of surface elements and mean curvature magnitude on surface nodes
		unsigned int fragmentCounter; // store how many fragments this object has created (keep consistent numbering of fragments)
		unsigned int updateCounter; // store how many times this object's implicit surface has been updated (keep consistent numbering of output files)
		double originalVolume; // store volume of the unfractured object
		std::vector<std::string> bndCnds;

		double updateVolume; // store the new volume when updating the collision shape
		Eigen::Vector3d vdbCOM, updateCOM; // updateCOM stores the shift of the center of mass due to a collision shape update, vdbCOM stores the offset of the collision mesh COM to the VDB COM (we want to avoid re-sampling the level-set to account for COM shifts of fragments)
		ColliderData collInUse, collUpdate; // store data required to properly delete stuff when changing the collision shape
		
		double contactDuration;
		vect3d_map contactTractions;

		template<typename GridPtrType>
		int initFractureSim(
			ColliderData& collMesh,
			GridPtrType vdbImplicitSurface,
			const Eigen::Vector3d& com,
			double crackMeshSize,
			double voxelSize, // should be consistent with the vdbImplicitSurface parameter
			const MaterialModel& matMdl,
			int    remesh,
			double originalVolume,
			FractureBEM* parent,
			int maxCracks = -1,
			double samplingDensity=1.0
		);
		int copyInteriorCrackMesh(FractureBEM* source, FractureBEM* target, id_map newCrackIDs );

		// these methods build a BEM mesh from the RB's collision shape
		int getMeshFromRB();
		int getMeshFromShape(btCollisionShape* shape, node_map& nodes, elem_map& elems);
		int getMeshFromBtInterface(btStridingMeshInterface* meshIntf, node_map& nodes, elem_map& elems);
		int getMeshFromConvHull(btConvexHullShape* convhull, node_map& nodes, elem_map& elems);
		int getMeshFromBox(btBoxShape* box, node_map& nodes, elem_map& elems);
		int initVDBfromSphere(double voxelSize, bool ignoreCrackSelfIntersect); // special case: for btSphereShape, init the VDB first
		
		void precompMeshData();
		void getBalancedContactTractions(vect3d_map& q);
	};
}

#endif	/* FRACTURERB_H */