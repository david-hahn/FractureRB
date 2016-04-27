/* 
 * File:   BulletWrapper.cpp
 * Author: David
 * 
 * Created on 08. Jan 2015
 */

#include "BulletWrapper.h"
#include "FractureRB.h"
#include "FractureBEM.h"

#include <ostream>
#include <fstream>
#include <sstream>
#include <omp.h> // just for omp_get_wtime

#include <btBulletDynamicsCommon.h>
#include <BulletWorldImporter/btBulletWorldImporter.h>
#include <BulletCollision/Gimpact/btGImpactCollisionAlgorithm.h>

#include <BulletDynamics/MLCPSolvers/btMLCPSolver.h>
#include <BulletDynamics/MLCPSolvers/btDantzigSolver.h>
//#include <BulletDynamics/MLCPSolvers/btLemkeSolver.h> //only in Bullet versions later than 2.82

using namespace std;

namespace FractureSim{
	enum TIMERS{
		T_START, // store start time
		T_INI, // scene loading
		T_RBD, // rigid body dynamics
		T_FRS // fracture simulation
	};
	
	BulletWrapper::BulletWrapper(bool useDefaultSolver) : outDir("") {
		broadphase = new btDbvtBroadphase();
		collisionConfiguration = new btDefaultCollisionConfiguration();
		dispatcher = new btCollisionDispatcher(collisionConfiguration);
		if( useDefaultSolver)
			solver = new btSequentialImpulseConstraintSolver();
		else
			solver = new btMLCPSolver( new btDantzigSolver() );
		dynamicsWorld = new btDiscreteDynamicsWorld(dispatcher, broadphase, solver, collisionConfiguration);
		btGImpactCollisionAlgorithm::registerAlgorithm(dispatcher);
		timings[T_INI]=0.0;
		timings[T_RBD]=0.0; counters[T_RBD]=0;
		timings[T_FRS]=0.0; counters[T_FRS]=0;
		timings[T_START]=omp_get_wtime();
	}

	BulletWrapper::~BulletWrapper(){
		printf("\n%% cleaning up Bullet data ... ");
		
//		printf("del names\n");
		for(int i=0; i<dynamicsWorld->getNumCollisionObjects(); ++i){
			string* name = (string*)(dynamicsWorld->getCollisionObjectArray()[i]->getUserPointer());
			if( name ){
//				printf("%d (%d) -> %s\n",i,dynamicsWorld->getCollisionObjectArray()[i]->getUserIndex(),name->c_str() );
				delete name;
			}
			dynamicsWorld->getCollisionObjectArray()[i]->setUserPointer(NULL);
		}
		
//		printf("del FRBs: %d\n", breakableRBs.size());
		for(std::map<int, FractureRB*>::iterator it = breakableRBs.begin(); it != breakableRBs.end(); ++it){
//			printf("rm %d\n",it->first);
			dynamicsWorld->removeRigidBody(it->second->getRB());
//			printf("del\n");
			delete it->second; // this will delete the contained RB if FRB has ownership, otherwise the RB belongs to the importer
		}
		breakableRBs.clear();
		
//		printf("clear importer\n");
		importer->deleteAllData(); //removes all objects created by importer from the world, then deletes them
//		printf("del importer\n");
		delete importer;

		if( dynamicsWorld->getNumCollisionObjects() )
			printf("\n%% deleting %d remaining objects from dynamics world",dynamicsWorld->getNumCollisionObjects());
		for(int i=dynamicsWorld->getNumCollisionObjects()-1; i>=0; --i){
			btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray()[i];
			if( obj ){
//				printf("remove\n");
				dynamicsWorld->removeCollisionObject(obj);
				btRigidBody* rb = dynamic_cast<btRigidBody*>(obj);
//				printf("motion state\n");
				if( rb && rb->getMotionState() ) delete rb->getMotionState();
//				printf("collision shape\n");
				if( obj->getCollisionShape() ) delete obj->getCollisionShape();
//				printf("del obj\n");
				delete obj;
			}
		}

		ColliderData::clearStoredData(false); // we've already delete the colliison shape for all the rigid bodies in the loop above

//		printf("del dynamics world and stuff\n");
		delete dynamicsWorld;
		delete solver;
		delete dispatcher;
		delete collisionConfiguration;
		delete broadphase;
		
		double t_all=omp_get_wtime()-timings[T_START];
		printf("\n%% Collective runtimes:");
		printf("\n%% ... %8.2lfs (%3.0lf%%)\t scene import and initialization",timings[T_INI],100.0*timings[T_INI]/t_all);
		printf("\n%% ... %8.2lfs (%3.0lf%%)\t rigid body dynamics\t avg. %.4lfs (%u)",timings[T_RBD],100.0*timings[T_RBD]/t_all,timings[T_RBD]/counters[T_RBD],counters[T_RBD]);
		printf("\n%% ... %8.2lfs (%3.0lf%%)\t fracture simulation\t avg. %.4lfs (%u)",timings[T_FRS],100.0*timings[T_FRS]/t_all,timings[T_FRS]/counters[T_FRS],counters[T_FRS]);
		printf("\n%% ... %8.2lfs       \t others",t_all-timings[T_INI]-timings[T_RBD]-timings[T_FRS]);
		printf("\n%% total %6.2lfs\n",t_all);
		printf("\n%% BulletWrapper destroyed successfully\n");
	}

	int BulletWrapper::initSceneFromFile(
		std::string filename,
		std::string paramfile,
		double splitImpulse, int useEstSIFs
	){
		double t0=omp_get_wtime();
		double minVoxelSize=DBL_MAX;
		map<string, vector<string> > params;
		loadParamFile(params, paramfile);
		doneTimesteps = 0;
		
		importer = new btBulletWorldImporter(dynamicsWorld);
		importer->setVerboseMode(0);
		importer->loadFile(filename.c_str());
		
		// convert BvhTri...meshes to GImpact meshes
		for(int i=0; i<dynamicsWorld->getNumCollisionObjects(); ++i){
			btRigidBody* rb = dynamic_cast<btRigidBody*>(dynamicsWorld->getCollisionObjectArray()[i]);
			if( rb ){
				btCollisionShape* shape = rb->getCollisionShape();
				btCompoundShape* cmpd;
				while(cmpd=dynamic_cast<btCompoundShape*>(shape)) shape=cmpd->getChildShape(0);
				if( !(rb->isStaticOrKinematicObject()) && dynamic_cast<btBvhTriangleMeshShape*>( shape )  ){
					btBvhTriangleMeshShape* bmesh = dynamic_cast<btBvhTriangleMeshShape*>( shape );
					btGImpactMeshShape* gimesh = new btGImpactMeshShape(bmesh->getMeshInterface()); //ToDo: delete this mesh object in the destructor
					gimesh->updateBound();
					dynamicsWorld->removeRigidBody( rb );
					rb->setCollisionShape(gimesh);
					dynamicsWorld->addRigidBody( rb );
					printf("\n%% Changed collision shape to GImpact mesh on %s", importer->getNameForPointer(rb));
				}
			}
		}

		for(int i=0; i<dynamicsWorld->getNumCollisionObjects(); ++i){
			btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray()[i];
			obj->setUserIndex(i); // user index is set to the position in the array
			obj->setUserPointer(new string( importer->getNameForPointer(obj) )); // will be used later to store the object name of breakable rigid bodies
			printf("\n%% object %d: \"%s\" ",i, importer->getNameForPointer(obj));

			btRigidBody* rb = dynamic_cast<btRigidBody*>(obj);
			if( rb ){
				// read general rb params (name, initial velocity and rotation)
				if(	params.count(importer->getNameForPointer(obj)) >0 ){
					vector<string>& paramList = params[importer->getNameForPointer(obj)];
					//param list format:
					// % collision-shape name;object name;ini. vel. x;init. vel. y;init. vel. z;init. rot. x;init. rot. y;init. rot. z;init. sleeping;mesh size;voxel size;Young's modulus;Possion ratio;density;strength;toughness;compressive factor;remesh target;ignore crack self-intersections;material model;max cracks per BEM run;fracture method;sampling density
					// name
					if( paramList.size()>0 ){ // switch the name to the one given in the param file
						if( obj->getUserPointer() ) delete obj->getUserPointer();
						obj->setUserPointer(new string(paramList[0]));
					}
					printf("--> \"%s\"", ((string*)obj->getUserPointer())->c_str() );
				
					// initial velocity & rotation
					if( paramList.size()>7 ){
						double tmp; int check, tmI;
						btVector3 v(0.,0.,0.), w(0.,0.,0.);
						check=sscanf(paramList[1].c_str(),"%lf",&tmp); if(check==1) v[0]=tmp;
						check=sscanf(paramList[2].c_str(),"%lf",&tmp); if(check==1) v[1]=tmp;
						check=sscanf(paramList[3].c_str(),"%lf",&tmp); if(check==1) v[2]=tmp;
						check=sscanf(paramList[4].c_str(),"%lf",&tmp); if(check==1) w[0]=tmp;
						check=sscanf(paramList[5].c_str(),"%lf",&tmp); if(check==1) w[1]=tmp;
						check=sscanf(paramList[6].c_str(),"%lf",&tmp); if(check==1) w[2]=tmp;
						check=sscanf(paramList[7].c_str(),"%d",&tmI); if(check==1 && tmI==1) rb->setActivationState(ISLAND_SLEEPING);
						if( check==1 && tmI==2 ){ // force to kinematic object
							btVector3 z; z.setZero();
							rb->setMassProps(0.0, z);
							rb->setCollisionFlags( rb->CF_KINEMATIC_OBJECT );
							rb->setActivationState( DISABLE_DEACTIVATION );
						}
						rb->setLinearVelocity (v);
						rb->setAngularVelocity(w);
					}
					// build breakable rigid body ...
					if( paramList.size()>18 && !rb->isStaticOrKinematicObject() ){
						double tmp, meshSize=1.0, voxelSize=0.1, E=1.0, nu=0.3, rho=1.0, Sc=1.0, Kc=1.0, cf=1.0, smplDens=1.0;
						bool noSI=false; int remesh=0, tmI, check, maxCracks=-1;
						string material = "default";
						FRACTURE_METHOD method = FULL_BEM;
						check=sscanf(paramList[8 ].c_str(),"%lf",&tmp); if(check==1) meshSize=tmp;
						check=sscanf(paramList[9 ].c_str(),"%lf",&tmp); if(check==1) voxelSize=tmp;
						check=sscanf(paramList[10].c_str(),"%lf",&tmp); if(check==1) E=tmp;
						check=sscanf(paramList[11].c_str(),"%lf",&tmp); if(check==1) nu=tmp;
						check=sscanf(paramList[12].c_str(),"%lf",&tmp); if(check==1) rho=tmp;
						check=sscanf(paramList[13].c_str(),"%lf",&tmp); if(check==1) Sc=tmp;
						check=sscanf(paramList[14].c_str(),"%lf",&tmp); if(check==1) Kc=tmp;
						check=sscanf(paramList[15].c_str(),"%lf",&tmp); if(check==1) cf=tmp;
						check=sscanf(paramList[16].c_str(),"%d", &tmI); if(check==1) remesh=tmI;
						check=sscanf(paramList[17].c_str(),"%lf",&tmp); if(check==1 && tmp>0.0) noSI=true;
						material =   paramList[18];
						// some optional params ...
						if( paramList.size()>19 ){
							check=sscanf(paramList[19].c_str(),"%d", &tmI); if(check==1) maxCracks=tmI;
						}
						if( paramList.size()>20 ){
							check=sscanf(paramList[20].c_str(),"%d", &tmI); if(check==1) method=(FRACTURE_METHOD) tmI;
						}
						if( paramList.size()>21 ){
							check=sscanf(paramList[21].c_str(),"%lf", &tmp); if(check==1) smplDens=tmp;
						}

						if(minVoxelSize > voxelSize) minVoxelSize = voxelSize; // store the smallest voxel-size in the scene
						printf("\n%% ... building breakable rigid body:"
							   "\n%% ... crack mesh size %.3lg / voxel size %.3lg = resolution ratio %.1lf"
							   "\n%% ... Young's modulus %.3lg, Poisson's ratio %.3lf, density %.3lg"
							   "\n%% ... tensile strength %.3lg, tensile toughness %.3lg, compressive factor %.3lf",
							   meshSize, voxelSize, meshSize/voxelSize, E,nu,rho, Sc,Kc,cf
						);/**/
						breakableRBs[i] = new FractureRB( rb ,false); //FractureRB class should not delete this rb since it's the importer's job to do so
						breakableRBs[i]->setOutputDir(outDir);
						if( useEstSIFs >=0 ) breakableRBs[i]->setEstSIFsThreshold( useEstSIFs );
						breakableRBs[i]->initFractureSim(
							meshSize,voxelSize,E,nu,rho,Sc,Kc,cf,remesh,noSI,material,maxCracks,method,smplDens
						);
					}
				}
				printf("\n%% ... is %ssleeping", (rb->getActivationState()==ACTIVE_TAG)?"not ":"");
				//printf("\n%% ... sleeping thresholds (%.3lf, %.3lf)", rb->getLinearSleepingThreshold(), rb->getAngularSleepingThreshold());
				if(!(rb->isStaticObject()) ) printf("\n%% ... mass is %.3le", 1.0/rb->getInvMass());
			}
			printf("\n%% ... is %s%s%s, shape type %d", !(obj->isStaticOrKinematicObject())?"active":"", obj->isKinematicObject()?"kinematic":"", obj->isStaticObject()?"static":"", obj->getCollisionShape()->getShapeType());
			printf("\n%% ... friction %.3lf, restitution %.3lf", obj->getFriction(), obj->getRestitution());
			//printf("%% -- collision shape %s is %s %s %s\n", obj->getCollisionShape()->getName(), obj->getCollisionShape()->isConvex()?"convex":"", obj->getCollisionShape()->isCompound()?"compound":"", obj->getCollisionShape()->isConcave()?"concave":"");
		}
		nextRB = dynamicsWorld->getNumCollisionObjects(); // store where to continue numbering rigid bodies -> goes into UserIndex
		btContactSolverInfo& info = dynamicsWorld->getSolverInfo();
		if(splitImpulse>0.0){
			info.m_splitImpulse = 1; //enable split impulse feature
			info.m_splitImpulsePenetrationThreshold = -splitImpulse*minVoxelSize;
			printf("\n%% split impulse threshold set to %.3lg (%.3lg * min. voxel size)",splitImpulse*minVoxelSize,splitImpulse);
		}else{
			info.m_splitImpulse = 0;
			printf("\n%% impulse splitting disabled");
		}
		timings[T_INI]+=(omp_get_wtime()-t0);
		return 0;
	}
	
	
	int BulletWrapper::stepSimulation(
		ostream* outPtr, OUT_FORMAT outFmt,
		double dt, double impulseThreshold, double forceThreshold
	){
		double t0=omp_get_wtime(), t1;
		// remember the velocities before taking the timestep, we'll need those for the collision handling
		map<int, btTransform> oldPosition;
		map<int, Eigen::Vector3d> oldLinVelocity, oldAngVelocity;
		storeMotionData(oldPosition, oldLinVelocity, oldAngVelocity);
		
		// limited support for kinematic objects (seems we can't export keyframes to the .bullet scene)
		// but we let the user specify constant linear and angular velocities and update these in every timestep
		for(int i=0; i<dynamicsWorld->getNumCollisionObjects(); ++i){
			btRigidBody* rb = dynamic_cast<btRigidBody*>(dynamicsWorld->getCollisionObjectArray()[i]);
			if( rb && rb->isKinematicObject() ){
				btTransform pt;
				rb->predictIntegratedTransform(dt, pt);
				rb->setCenterOfMassTransform( pt );
			}
		}
		dynamicsWorld->stepSimulation(dt, 1, dt); // force exactly 1 sub-step of the given duration
		++doneTimesteps; ++counters[T_RBD];
		printf("%%");
//
//		// produce output for current state
//		if( outPtr!=NULL ){
//			if( outFmt == OUT_MEL ){
//				ostream& melfile = (*outPtr); // for convenience
//				melfile << "currentTime " << doneTimesteps << " ;" << endl;
//				//printf("%% generating MEL output ...\n");
//				for(int i=0; i<dynamicsWorld->getNumCollisionObjects(); ++i){
//					btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray().at(i);
//					if(!obj->isStaticObject() ){
//						string& name = (*(string*)obj->getUserPointer());
//						melfile << "select -r " << name << " ;" << endl;
//						btScalar m[16]; m[15]=1.0;
//						obj->getWorldTransform().getOpenGLMatrix(m);
//						melfile << "xform -ws -m ";
//						for(int k=0; k<16; ++k) melfile << m[k] << " ";
//						melfile << " ;" << endl;
//						melfile << "setKeyframe \"" << name << ".tx\";" << endl;
//						melfile << "setKeyframe \"" << name << ".ty\";" << endl;
//						melfile << "setKeyframe \"" << name << ".tz\";" << endl;
//						melfile << "setKeyframe \"" << name << ".rx\";" << endl;
//						melfile << "setKeyframe \"" << name << ".ry\";" << endl;
//						melfile << "setKeyframe \"" << name << ".rz\";" << endl;
//					}
//				}
//			}else
//			if( outFmt == OUT_BPY ){
//				ostream& pyfile = (*outPtr);
//				pyfile << "bpy.context.scene.frame_set(" << doneTimesteps << ")" << endl;
//				for(int i=0; i<dynamicsWorld->getNumCollisionObjects(); ++i){
//					btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray().at(i);
//					if(!obj->isStaticObject() ){
//						string& name = (*(string*)obj->getUserPointer());
//						pyfile << "bpy.context.scene.objects.active = bpy.data.objects[\"" << name << "\"]" << endl;
//						btScalar m[16]; m[15]=1.0;
//						obj->getWorldTransform().getOpenGLMatrix(m);
//						pyfile << "bpy.context.object.matrix_local = Matrix(((";
//						//1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)))
//						for(int k=0; k<4; ++k) for(int l=0; l<4; ++l) pyfile << ((k>0 && l==0)?",(":"") << m[4*l+k] << ((l==3)?")":","); //       
//						pyfile << "))" << endl;
//						pyfile << "x=bpy.context.active_object.keyframe_insert(\"location\")" << endl;
//						pyfile << "x=bpy.context.active_object.keyframe_insert(\"rotation_euler\")" << endl;
//						pyfile << "x=bpy.context.active_object.keyframe_insert(\"scale\")" << endl;
//					}
//				}
//			}
//		}
		
		
		// handle collisions and call the fracture simulation when appropriate
		// could do this in a callback -- should not matter as long as we don't allow bullet to sub-step
		// http://www.bulletphysics.org/mediawiki-1.5.8/index.php/Simulation_Tick_Callbacks
		int count=0; // count notable collisions for return value
		// collect data about collisions in these maps, the key is the object id
		map<int, vector<unsigned int> >    collElems;
		map<int, vector<Eigen::Vector3d> > collDirections;
		map<int, vector<double> >          collImpulses, collDurations;
		
		// collect collision info from contact manifolds
		int numManifolds = dynamicsWorld->getDispatcher()->getNumManifolds();
		for(int i=0; i < numManifolds; i++){
			btPersistentManifold* contactManifold = dynamicsWorld->getDispatcher()->getManifoldByIndexInternal(i);
			const btRigidBody* rbA = dynamic_cast<const btRigidBody*>(contactManifold->getBody0());
			const btRigidBody* rbB = dynamic_cast<const btRigidBody*>(contactManifold->getBody1());
			//printf("%% %d: contact between %s and %s:\n", doneTimesteps, ((string*)rbA->getUserPointer())->c_str(), ((string*)rbB->getUserPointer())->c_str());
			
			for(int j=0; j < contactManifold->getNumContacts(); j++){
				btManifoldPoint& pt = contactManifold->getContactPoint(j);
				btVector3 dir = pt.m_normalWorldOnB;//pt.getPositionWorldOnB() - pt.getPositionWorldOnA();
				//dir.normalize(); // should be normalized already
				
				
				//get contact duration using Hertzian contact (sphere-sphere)
				// - at least one object must be a breakable rigid body
				// - if the other object is not breakable use it's bounding sphere radius (or distance from COM to contact pt?)
				// - for breakable RBs we can find an element-ID closest to the contact point we store this to register the collision with the BEM solver
				// - we'll also get a local mean curvature from the BEM mesh to estimate the radius of the model sphere
				// - the contact duration depends on:
				// -- the kinetic energy (in a center-of-mass frame) MV^2/2 with M=m1+m2 and V=v1-v2
				// -- the elasticity of both objects (it's ok if one is rigid) 1/E = (1-nu1^2)/E1 + (1-nu2^2)/E2
				// -- the radii of the sphere's in contact (it's ok if one is a plane) 1/R = 1/r1 + 1/r2
				// -- then we have t = 2.94*X/V, X=( 15/16*M*V^2/(E*sqrt(R)) )^(2/5), see Chap. 5 Problem 3 in Contact Mechanics and Friction, Springer 2010 http://www.springer.com/materials/mechanics/book/978-3-642-10802-0
				// --- http://www.springer.com/cda/content/document/cda_downloaddocument/9783642108020-c1.pdf (Jan.2015)
				// -- the pre-factor c=2.94*(15/16)^(2/5) can be computed as (wolfram alpha) "integrate[ 1/sqrt(1-x^(5/2)), {x, 0,1}] * 2 * (15/16)^(2/5)" = 2.8682656991953313822367...
				// -- with this we can reduce to the form found in the papers by Glondu et al. (2013) and Koschier et al. (2014)
				// -- t_c = c*( M^2 / ( E^2 * R * V ) )^(1/5)

				double e=0.0, m=0.0, r=0.0, v=0.0, t_c=0.0;
				unsigned int elemA, elemB;
				Eigen::Vector3d pOnA=toEigen(pt.m_localPointA), pOnB=toEigen(pt.m_localPointB), vOnA, vOnB;

				if( breakableRBs.count( rbA->getUserIndex() )){
					r += std::abs( breakableRBs[rbA->getUserIndex()]->
						 getCurvatureAndElementID(pOnA, elemA) );
					e += 1.0 / breakableRBs[rbA->getUserIndex()]->getEovNuSq();
				}
				if( breakableRBs.count( rbB->getUserIndex() )){
					r += std::abs( breakableRBs[rbB->getUserIndex()]->
						 getCurvatureAndElementID(pOnB, elemB) );
					e += 1.0 / breakableRBs[rbB->getUserIndex()]->getEovNuSq();
				}
				// if one of the objects is not breakable we treat is as a locally planar rigid body
				// consquently e is not changed because it's rigid
				// and r is not changed because we treat is as planar -- so just make everything that's sharp breakable in the scene!
				m += rbA->getInvMass(); // note that inv.mass is set to 0 for static objects by the Bullet system
				m += rbB->getInvMass();
				// the effective radius, elasticity and mass are all computed from reciprocal summation
				r=1.0/r;
				e=1.0/e;
				m=1.0/m;
				
				// next we need the velocity of approach -- we need to take this from the previous RB timestep, since the collision has already been processed once we get here
				// following btRigidBody::getVelocityInLocalPoint but for the pre-collision velocities
				vOnA = oldLinVelocity[rbA->getUserIndex()] + oldAngVelocity[rbA->getUserIndex()].cross(pOnA);
				vOnB = oldLinVelocity[rbB->getUserIndex()] + oldAngVelocity[rbB->getUserIndex()].cross(pOnB);
				// project velocity difference onto contact normal direction
				v = std::abs( ( vOnA - vOnB ).dot( toEigen(dir).normalized() ) );
				// now we can compute the contact duration
				t_c = std::pow( m*m / ( e*e * r * v ) , 0.2 )  *2.8682656991953313822367; // see Chap. 5 Problem 3 in Contact Mechanics and Friction, Springer 2010 to compute the constant
				//printf("%% ... duration %.3les: eff. m=%.1le, r=%.3lf, E=%.1le, v=%.1le\n",t_c,m,r,e,v);

				if( breakableRBs.count( rbA->getUserIndex() )){
					collElems[rbA->getUserIndex()].push_back( elemA );
					collDirections[rbA->getUserIndex()].push_back( toEigen( rbA->getWorldTransform().inverse().getBasis()*   dir  ));
					collImpulses[rbA->getUserIndex()].push_back( pt.getAppliedImpulse() );
					collDurations[rbA->getUserIndex()].push_back( t_c );
				}
				if( breakableRBs.count( rbB->getUserIndex() )){
					collElems[rbB->getUserIndex()].push_back( elemB );
					collDirections[rbB->getUserIndex()].push_back( toEigen( rbB->getWorldTransform().inverse().getBasis()* (-dir) ));
					collImpulses[rbB->getUserIndex()].push_back( pt.getAppliedImpulse() );
					collDurations[rbB->getUserIndex()].push_back( t_c );
				}
			}
		}

		set<FractureRB*>  largeFragments;
		set<btRigidBody*> smallFragments;
		set<int> runTheseSims;

		// run fracture sim for each breakable RB exceeding the impulse threshold
		for(map<int, vector<unsigned int> >::iterator it = collElems.begin(); it != collElems.end(); ++it ){
			if(breakableRBs.count(it->first)) {
				double totalImpulse=0.0, totalForce;
				for(int i=0; i < it->second.size(); ++i) totalImpulse += collImpulses[it->first][i];
				//if( totalImpulse >= impulseThreshold ){
					//printf("\n%% %d: collisions for %s, total impulse is %.4lg ...\n",doneTimesteps, ((string*)breakableRBs[it->first]->getRB()->getUserPointer())->c_str(),totalImpulse);
					//printf("%% collision points for object %d (%s):\n",it->first,((string*)breakableRBs[it->first]->getRB()->getUserPointer())->c_str() );
					//printf("%% collision points for object %d:\n",it->first );
				breakableRBs[it->first]->clearContacts(); // clear all old contact points
				for(int i=0; i < it->second.size(); ++i){ // add the current contact points
					breakableRBs[it->first]->addContact(it->second[i], collDirections[it->first][i],collImpulses[it->first][i],collDurations[it->first][i]);
				}
				totalForce=breakableRBs[it->first]->getTotalContactForce();
				if((totalImpulse >= impulseThreshold || totalForce >= forceThreshold) && totalForce > 0.0){
					printf("\n%% %d: I %-10.4lg  F %-10.4lg (%s)",doneTimesteps,totalImpulse,totalForce, ((string*)breakableRBs[it->first]->getRB()->getUserPointer())->c_str());
					runTheseSims.insert( it->first );
				}
			}
		}
		// restore the motion data of the previous time-step since fragments will be created such that
		// they match the parent's position and velocity, we'll re-compute the time-step once all the fragments are in place
		if(!runTheseSims.empty()){
			printf("\n%% ************************************************************");
			restoreMotionData(oldPosition, oldLinVelocity, oldAngVelocity);
		}
		t1=omp_get_wtime();
		timings[T_RBD]+=(t1-t0);
		t0=t1;
		for(set<int>::iterator it = runTheseSims.begin(); it != runTheseSims.end(); ++it){
			int check = breakableRBs[*it]->runFractureSim(dt,doneTimesteps); //ToDo: make all fracture sims run in parallel (?)
			if( check > 0 ){ // if check is positive, the fracture sim has added new cracks ...
				// ... so we look for new fragments
				check = breakableRBs[*it]->splitFragments(largeFragments, smallFragments);

				// add all the large fragments to the dynamics world
				for(set<FractureRB*>::iterator rbit=largeFragments.begin(); rbit!=largeFragments.end(); ++rbit){
					(*rbit)->getRB()->setUserIndex( nextRB );
					dynamicsWorld->addRigidBody( (*rbit)->getRB() );
					breakableRBs[nextRB]=(*rbit);
					++nextRB;

					// write a MEL import command for the new fragment
					if( outPtr!=NULL ){
						if( outFmt == OUT_MEL ){
							ostream& melfile = (*outPtr);
							string& name = *((string*)(*rbit)->getRB()->getUserPointer());
							melfile << "file -import -type \"OBJ\" -ignoreVersion -mergeNamespacesOnClash true \""
									<< breakableRBs[*it]->getOutputDir()
									<< name << ".obj\";" << endl;
							melfile << "rename \"Mesh\" \"" << name << "\";" << endl;
							//set visibility keys such that the fragment becomes visible in the next RB timestep
							melfile << "setAttr \"" << name << ".visibility\" 0;" << endl; // invisible before
							melfile << "setKeyframe \"" << name << ".v\";"  << endl;
							melfile << "currentTime " << (doneTimesteps+1) << ";" << endl;
							melfile << "setAttr \"" << name << ".visibility\" 1;" << endl; // visible after switch
							melfile << "setKeyframe \"" << name << ".v\";"  << endl;
							melfile << "currentTime " << doneTimesteps << ";" << endl;
						} else
						if( outFmt == OUT_BPY ){
							ostream& pyfile = (*outPtr);
							string& name = *((string*)(*rbit)->getRB()->getUserPointer());
							pyfile << "bpy.ops.import_scene.obj(filepath=\""
								   << breakableRBs[*it]->getOutputDir()
								   << name << ".obj\")" << endl;
							pyfile << "bpy.context.selected_objects[0].name = \"" << name << "\"" << endl;
							pyfile << "bpy.context.scene.objects.active = bpy.data.objects[\"" << name << "\"]" << endl;
							pyfile << "bpy.context.active_object.hide=True" << endl; // set NOT visible
							pyfile << "x=bpy.context.active_object.keyframe_insert(\"hide\")" << endl;
							pyfile << "bpy.context.active_object.hide=False" << endl; // set visible after current timestep
							pyfile << "x=bpy.context.active_object.keyframe_insert(\"hide\",frame=" << (doneTimesteps+1) << ")" << endl;
						}
					}
				}

				// add all the small fragments to the dynamics world
				for(set<btRigidBody*>::iterator rbit=smallFragments.begin(); rbit!=smallFragments.end(); ++rbit){
					(*rbit)->setUserIndex( nextRB );
					dynamicsWorld->addRigidBody( (*rbit) );
					++nextRB;

					// write a MEL import command for the new fragment
					if( outPtr!=NULL ){
						if( outFmt == OUT_MEL ){
							ostream& melfile = (*outPtr);
							string& name = *((string*)(*rbit)->getUserPointer());
							melfile << "file -import -type \"OBJ\" -ignoreVersion -mergeNamespacesOnClash true \""
									<< breakableRBs[*it]->getOutputDir()
									<< name << ".obj\";" << endl;
							melfile << "rename \"Mesh\" \"" << name << "\";" << endl;
							//set visibility keys such that the fragment becomes visible in the next RB timestep
							melfile << "setAttr \"" << name << ".visibility\" 0;" << endl; // invisible before
							melfile << "setKeyframe \"" << name << ".v\";"  << endl;
							melfile << "currentTime " << (doneTimesteps+1) << ";" << endl;
							melfile << "setAttr \"" << name << ".visibility\" 1;" << endl; // visible after switch
							melfile << "setKeyframe \"" << name << ".v\";"  << endl;
							melfile << "currentTime " << doneTimesteps << ";" << endl;
						} else
						if( outFmt == OUT_BPY ){
							ostream& pyfile = (*outPtr);
							string& name = *((string*)(*rbit)->getUserPointer());
							pyfile << "bpy.ops.import_scene.obj(filepath=\""
								   << breakableRBs[*it]->getOutputDir()
								   << name << ".obj\")" << endl; //newly loaded object is selected but not active
							pyfile << "bpy.context.selected_objects[0].name = \"" << name << "\"" << endl;
							pyfile << "bpy.context.scene.objects.active = bpy.data.objects[\"" << name << "\"]" << endl;
							pyfile << "bpy.context.active_object.hide=True" << endl; // set NOT visible
							pyfile << "x=bpy.context.active_object.keyframe_insert(\"hide\")" << endl;
							pyfile << "bpy.context.active_object.hide=False" << endl; // set visible after current timestep
							pyfile << "x=bpy.context.active_object.keyframe_insert(\"hide\",frame=" << (doneTimesteps+1) << ")" << endl;
						}
					}
				}

				// handle changes to the parent object
				if( check > 0 ){ // update the parent object
					// the collision shape might have changed, so we need to remove and re-add the object to the dynamics world
					if( breakableRBs[*it]->pendingCollisionUpdate() ){
						//printf("\n%% removing RB from world ... ");
						//printf("\n%% pre-update COM position:  %.4f %.4f %.4f (%s)",
						//	breakableRBs[*it]->getRB()->getCenterOfMassPosition()[0],
						//	breakableRBs[*it]->getRB()->getCenterOfMassPosition()[1],
						//	breakableRBs[*it]->getRB()->getCenterOfMassPosition()[2],
						//	((string*)(breakableRBs[*it]->getRB())->getUserPointer())->c_str());
						dynamicsWorld->removeRigidBody( breakableRBs[*it]->getRB() );
						breakableRBs[*it]->doCollisionUpdate();
						//printf("\n%% adding RB to world ... ");
						dynamicsWorld->addRigidBody( breakableRBs[*it]->getRB() );
						//printf("\n%% post-update COM position: %.4f %.4f %.4f (%s)",
						//	breakableRBs[*it]->getRB()->getCenterOfMassPosition()[0],
						//	breakableRBs[*it]->getRB()->getCenterOfMassPosition()[1],
						//	breakableRBs[*it]->getRB()->getCenterOfMassPosition()[2],
						//	((string*)(breakableRBs[*it]->getRB())->getUserPointer())->c_str());
						//printf("\n%% update complete\n");
					}
							
					//update the visual mesh of the parent, as it will have new (but possibly incomplete) cracks ...
					string& name = *((string*)breakableRBs[*it]->getRB()->getUserPointer());
					if( outPtr!=NULL ){
						if( outFmt == OUT_MEL ){
							ostream& melfile = (*outPtr);
							// switch the old mesh invisible
							melfile << "setAttr \"" << name << ".visibility\" 1;" << endl; // visible before
							melfile << "setKeyframe \"" << name << ".v\";"  << endl;
							melfile << "currentTime " << (doneTimesteps+1) << ";" << endl;
							melfile << "setAttr \"" << name << ".visibility\" 0;" << endl; // invisible after switch
							melfile << "setKeyframe \"" << name << ".v\";"  << endl;
							melfile << "currentTime " << doneTimesteps << ";" << endl;
							// load the new mesh and set visibility keys
							melfile << "file -import -type \"OBJ\" -ignoreVersion -mergeNamespacesOnClash true \""
									<< breakableRBs[*it]->getOutputDir()
									<< name << "_" << check << ".obj\";" << endl;
						} else
						if( outFmt == OUT_BPY ){
							ostream& pyfile = (*outPtr);
							pyfile << "bpy.context.scene.objects.active = bpy.data.objects[\"" << name << "\"]" << endl;
							pyfile << "bpy.context.active_object.hide=False" << endl; // set visible
							pyfile << "x=bpy.context.active_object.keyframe_insert(\"hide\")" << endl;
							pyfile << "bpy.context.active_object.hide=True" << endl; // set NOT visible after current timestep
							pyfile << "x=bpy.context.active_object.keyframe_insert(\"hide\",frame=" << (doneTimesteps+1) << ")" << endl;
							pyfile << "bpy.ops.import_scene.obj(filepath=\""
								   << breakableRBs[*it]->getOutputDir()
								   << name << "_" << check << ".obj\")" << endl; //newly loaded object is selected but not active
						}
					}
					// rename the current object to match the new mesh name
					stringstream checkstr; checkstr<<check;
					name.append("_").append(checkstr.str());
					if( outPtr!=NULL ){
						if( outFmt == OUT_MEL ){
							ostream& melfile = (*outPtr);
							// rename the new mesh
							melfile << "rename \"Mesh\" \"" << name << "\";" << endl;
							//set visibility keys such that the fragment becomes visible in the next RB timestep
							melfile << "setAttr \"" << name << ".visibility\" 0;" << endl; // invisible before
							melfile << "setKeyframe \"" << name << ".v\";"  << endl;
							melfile << "currentTime " << (doneTimesteps+1) << ";" << endl;
							melfile << "setAttr \"" << name << ".visibility\" 1;" << endl; // visible after switch
							melfile << "setKeyframe \"" << name << ".v\";"  << endl;
							melfile << "currentTime " << doneTimesteps << ";" << endl;
						} else
						if( outFmt == OUT_BPY ){
							ostream& pyfile = (*outPtr);
							pyfile << "bpy.context.selected_objects[0].name = \"" << name << "\"" << endl;
							pyfile << "bpy.context.scene.objects.active = bpy.data.objects[\"" << name << "\"]" << endl;
							pyfile << "bpy.context.active_object.hide=True" << endl; // set NOT visible
							pyfile << "x=bpy.context.active_object.keyframe_insert(\"hide\")" << endl;
							pyfile << "bpy.context.active_object.hide=False" << endl; // set visible after current timestep
							pyfile << "x=bpy.context.active_object.keyframe_insert(\"hide\",frame=" << (doneTimesteps+1) << ")" << endl;
						}
					}
				}else if( check==0 ){ // remove the parent object
					//set visibility keys such that the fragment becomes visible in the next RB timestep
					// write a MEL command to hide the old object
					if( outPtr!=NULL ){
						string& name = *((string*)breakableRBs[*it]->getRB()->getUserPointer());
						if( outFmt == OUT_MEL ){
							ostream& melfile = (*outPtr);
							melfile << "setAttr \"" << name << ".visibility\" 1;" << endl; // visible before
							melfile << "setKeyframe \"" << name << ".v\";"  << endl;
							melfile << "currentTime " << (doneTimesteps+1) << ";" << endl;
							melfile << "setAttr \"" << name << ".visibility\" 0;" << endl; // invisible after switch
							melfile << "setKeyframe \"" << name << ".v\";"  << endl;
							melfile << "currentTime " << doneTimesteps << ";" << endl;
						} else
						if( outFmt == OUT_BPY ){
							ostream& pyfile = (*outPtr);
							pyfile << "bpy.context.scene.objects.active = bpy.data.objects[\"" << name << "\"]" << endl;
							pyfile << "bpy.context.active_object.hide=False" << endl; // set visible
							pyfile << "x=bpy.context.active_object.keyframe_insert(\"hide\")" << endl;
							pyfile << "bpy.context.active_object.hide=True" << endl; // set NOT visible after current timestep
							pyfile << "x=bpy.context.active_object.keyframe_insert(\"hide\",frame=" << (doneTimesteps+1) << ")" << endl;

						}
					}
					dynamicsWorld->removeRigidBody( breakableRBs[*it]->getRB() );
					delete breakableRBs[*it];
					breakableRBs.erase( *it );
				}else{ // something went wrong - print a warning
					//printf("%% !!! splitFragments returned %d for %s!\n",check,((string*)breakableRBs[*it]->getRB()->getUserPointer())->c_str());
				}
			}else{
				// check==0 means we didn't add any elements, but there was no error
				if(check<0) printf("\n%% !!! runFractureSim returned %d for %s!",check,((string*)breakableRBs[*it]->getRB()->getUserPointer())->c_str());
			}
			++count;
		}
		t1=omp_get_wtime();
		timings[T_FRS]+=(t1-t0); counters[T_FRS]+=runTheseSims.size();
		t0=t1;

		if(!runTheseSims.empty()){
			printf("\n%% ************************************************************ %d\n", doneTimesteps);
			for(int i=0; i<dynamicsWorld->getNumCollisionObjects(); ++i){
				btRigidBody* rb = dynamic_cast<btRigidBody*>(dynamicsWorld->getCollisionObjectArray()[i]);
				if( rb && rb->isKinematicObject() ){
					btTransform pt;
					rb->predictIntegratedTransform(dt, pt);
					rb->setCenterOfMassTransform( pt );
				}
			}
			dynamicsWorld->stepSimulation(dt, 1, dt); // re-do the timestep once fragments have been handled
			++counters[T_RBD];
		}


		// produce output for current state
		if( outPtr!=NULL ){
			if( outFmt == OUT_MEL ){
				ostream& melfile = (*outPtr); // for convenience
				melfile << "currentTime " << doneTimesteps << " ;" << endl;
				//printf("%% generating MEL output ...\n");
				for(int i=0; i<dynamicsWorld->getNumCollisionObjects(); ++i){
					btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray().at(i);
					if(!obj->isStaticObject() ){
						string& name = (*(string*)obj->getUserPointer());
						melfile << "select -r " << name << " ;" << endl;
						btScalar m[16]; m[15]=1.0;
						obj->getWorldTransform().getOpenGLMatrix(m);
						melfile << "xform -ws -m ";
						for(int k=0; k<16; ++k) melfile << m[k] << " ";
						melfile << " ;" << endl;
						melfile << "setKeyframe \"" << name << ".tx\";" << endl;
						melfile << "setKeyframe \"" << name << ".ty\";" << endl;
						melfile << "setKeyframe \"" << name << ".tz\";" << endl;
						melfile << "setKeyframe \"" << name << ".rx\";" << endl;
						melfile << "setKeyframe \"" << name << ".ry\";" << endl;
						melfile << "setKeyframe \"" << name << ".rz\";" << endl;
					}
				}
				melfile << "select -cl;" << endl; // for convenience clear the selection at the end of operations to speed up Maya's preview during script exec.
			}else
			if( outFmt == OUT_BPY ){
				ostream& pyfile = (*outPtr);
				pyfile << "bpy.context.scene.frame_set(" << doneTimesteps << ")" << endl;
				for(int i=0; i<dynamicsWorld->getNumCollisionObjects(); ++i){
					btCollisionObject* obj = dynamicsWorld->getCollisionObjectArray().at(i);
					if(!obj->isStaticObject() ){
						string& name = (*(string*)obj->getUserPointer());
						pyfile << "bpy.context.scene.objects.active = bpy.data.objects[\"" << name << "\"]" << endl;
						btScalar m[16]; m[15]=1.0;
						obj->getWorldTransform().getOpenGLMatrix(m);
						pyfile << "bpy.context.object.matrix_local = Matrix(((";
						//1,0,0,0),(0,1,0,0),(0,0,1,0),(0,0,0,1)))
						for(int k=0; k<4; ++k) for(int l=0; l<4; ++l) pyfile << ((k>0 && l==0)?",(":"") << m[4*l+k] << ((l==3)?")":","); //       
						pyfile << "))" << endl;
						pyfile << "x=bpy.context.active_object.keyframe_insert(\"location\")" << endl;
						pyfile << "x=bpy.context.active_object.keyframe_insert(\"rotation_euler\")" << endl;
						pyfile << "x=bpy.context.active_object.keyframe_insert(\"scale\")" << endl;
					}
				}
			}
		}
//		if( outPtr!=NULL ){
//			if( outFmt == OUT_MEL ){
//				ostream& melfile = (*outPtr);
//				melfile << "select -cl;" << endl; // for convenience clear the selection at the end of operations to speed up Maya's preview during script exec.
//			}
//		}
		t1=omp_get_wtime();
		timings[T_RBD]+=(t1-t0);
		t0=t1;
		fflush(stdout); outPtr->flush();
		return count;
	}
	
	
	void BulletWrapper::storeMotionData(
		map<int, btTransform>& pos,
		map<int, Eigen::Vector3d>& linvel,
		map<int, Eigen::Vector3d>& angvel
	){
		for(int i=0; i<dynamicsWorld->getNumCollisionObjects(); ++i){
			btRigidBody* rb = dynamic_cast<btRigidBody*>(dynamicsWorld->getCollisionObjectArray()[i]);
			if(rb!=NULL){
				pos   [rb->getUserIndex()]=rb->getCenterOfMassTransform();
				linvel[rb->getUserIndex()]=toEigen(rb->getLinearVelocity());
				angvel[rb->getUserIndex()]=toEigen(rb->getAngularVelocity());
			}
		}
	}

	void BulletWrapper::restoreMotionData(
		map<int, btTransform>& pos,
		map<int, Eigen::Vector3d>& linvel,
		map<int, Eigen::Vector3d>& angvel
	){
		for(int i=0; i<dynamicsWorld->getNumCollisionObjects(); ++i){
			btRigidBody* rb = dynamic_cast<btRigidBody*>(dynamicsWorld->getCollisionObjectArray()[i]);
			if(rb!=NULL){
				rb->setAngularVelocity( fromEigen(angvel[rb->getUserIndex()]) );
				rb->setLinearVelocity ( fromEigen(linvel[rb->getUserIndex()]) );
				rb->setCenterOfMassTransform( pos[rb->getUserIndex()] );
				rb->clearForces();
			}
		}
	}

	void BulletWrapper::setOutputDir(std::string newOutDir){
		outDir = newOutDir;
		for(std::map<int, FractureRB*>::iterator it = breakableRBs.begin(); it != breakableRBs.end(); ++it){
			it->second->setOutputDir(outDir);
		}
	}

	Eigen::Vector3d toEigen(const btVector3& in){
		return Eigen::Vector3d(in[0],in[1],in[2]);
	}
	btVector3 fromEigen(const Eigen::Vector3d& in){
		return btVector3(in[0],in[1],in[2]);
	}


	// string-trim helpers from http://www.cplusplus.com/faq/sequences/strings/trim/
	std::string& trim_right_inplace( std::string& s,const std::string& delimiters = " \f\n\r\t\v" ){
		return s.erase( s.find_last_not_of( delimiters ) + 1 );
	}

	std::string& trim_left_inplace( std::string& s, const std::string& delimiters = " \f\n\r\t\v" ){
		return s.erase( 0, s.find_first_not_of( delimiters ) );
	}

	std::string& trim_inplace( std::string& s, const std::string& delimiters = " \f\n\r\t\v" ){
		return trim_left_inplace( trim_right_inplace( s, delimiters ), delimiters );
	}

	void BulletWrapper::loadParamFile(
		map<string, vector<string> >& params,
		string filename
	){
		params.clear();
		ifstream in(filename.c_str());
		string line, token, name;
		bool first,comment;
		while(in.good()){
			getline(in,line);
			if(line.empty() || line.at(0)=='%'){ // comment -- ignore
			}else{
				stringstream tokenize(line);
				first=true;
				comment=false;
				while(tokenize.good() && !comment){
					getline(tokenize,token,';');
					trim_inplace(token);
					if(first){
						name=token;
						first=false;
					}else{
						if(!token.empty()) if(token.at(0)=='%') comment=true; //allow commenting out rest of line
						if(!comment) params[name].push_back(token);
						//printf("added token \"%s\" to %s\n",token.c_str(), name.c_str());
					}
				}
			}
		}
		in.close();
	}

}
