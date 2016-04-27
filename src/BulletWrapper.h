/* 
 * File:   BulletWrapper.h
 * Author: David
 *
 * Created on 08. Jan 2015
 */

#ifndef BULLETWRAPPER_H
#define	BULLETWRAPPER_H

#include "types.h"

#include <string>
#include <iosfwd>

//fwd decl
class btBroadphaseInterface;
class btDefaultCollisionConfiguration;
class btCollisionDispatcher;
class btSequentialImpulseConstraintSolver;
class btDiscreteDynamicsWorld;
class btVector3;
class btBulletWorldImporter;
class btTransform;

namespace FractureSim{
	//fwd decl
	class FractureRB;
	
	// helper-fcn
	Eigen::Vector3d toEigen(const btVector3& in);
	btVector3 fromEigen(const Eigen::Vector3d& in);

    class BulletWrapper { // this class encapsulates a Bullet dynamics world
    public:
		enum OUT_FORMAT{
			OUT_MEL, // write keyframe output in Autodesk Maya MEL script format
			OUT_BPY // same for Blender-Python
			// ToDo? support other formats (could be Maya Python, Blender Python, Alembic, FBX, ...)
		};
		
		// constructs a new BulletWrapper with an empty dynamics world and default configuration
		BulletWrapper( bool useDefaultSolver=false );
		// deletes all contained Bullet objects
		virtual ~BulletWrapper();
		
		// load a scene from a .bullet file into the dynamics world
		// the fracture parameters for each object in the scene are specified
		// in a csv-formatted parameter file
		// impulse splitting threshold is specified relative to the smallest voxel size of any breakable rigid body
		// fracture simulations will be switched to fast estimated SIF evaluations
		// if the exceed the number of elements given in useEstSIFs (default -1 == never)
		int initSceneFromFile(
			std::string bulletfile,
			std::string paramfile,
			double splitImpulse=3.0,
			int useEstSIFs=-1 );
		
		// take one timestep of the rigid body simulation
		// the duration of the timestep is dt
		// all collisions involving breakable rigid bodies where
		// the collision implulse exceeds impulseThreshold
		// OR the total deformational force exceeds forceThreshold are
		// forwarded to the fracture simulation
		int stepSimulation(
			std::ostream* out=NULL, OUT_FORMAT outFmt = OUT_MEL, double dt=1.0/250.0,
			double impulseThreshold=DBL_EPSILON, double forceThreshold=DBL_EPSILON
		);

		void setOutputDir(std::string outDir);
		
	protected:
		btBroadphaseInterface* broadphase;
		btDefaultCollisionConfiguration* collisionConfiguration;
		btCollisionDispatcher* dispatcher;
		btSequentialImpulseConstraintSolver* solver;
		btDiscreteDynamicsWorld* dynamicsWorld;
		btBulletWorldImporter* importer;
		
		int doneTimesteps; // count how many times stepSimulation was called
		int nextRB; // keep a unique numbering on rigid bodies
		std::string outDir; // directory and file prefix for output files
		std::map<int, FractureRB*> breakableRBs;

		std::map<unsigned int, double> timings;
		id_map counters;
		
		void loadParamFile(
			std::map<std::string, std::vector<std::string> >& params,
			std::string filename
		);
		void storeMotionData(
			std::map<int, btTransform>& pos,
			std::map<int, Eigen::Vector3d>& linvel,
			std::map<int, Eigen::Vector3d>& angvel
		);
		void restoreMotionData(
			std::map<int, btTransform>& pos,
			std::map<int, Eigen::Vector3d>& linvel,
			std::map<int, Eigen::Vector3d>& angvel
		);
	};
}

#endif	/* BULLETWRAPPER_H */