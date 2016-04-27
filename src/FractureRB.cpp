/* 
 * File:   FractureRB.cpp
 * Author: David
 * 
 * Created on 08. Jan 2015
 */

#include "FractureRB.h"
#include "FractureBEM.h"
#include "MaterialModel.h"

#include <omp.h> // just for omp_get_wtime
#include "PostProcessor.h" // just for VTK debug output
#include "SubsampledCrackTip.h" // just for SIF debug output
//#include "Reader.h" // just for printMap debug output

#include <btBulletDynamicsCommon.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
#include <LinearMath/btConvexHullComputer.h>

using namespace std;

namespace FractureSim{
	FractureRB::FractureRB()
		: rb(NULL), fractSim(NULL), mat(NULL), outDir(""), fragmentCounter(0),
		  updateCounter(0), haveMesh(false), haveParams(false), maxCracks(-1),
		  ownRB(false), haveBEM(false), remeshTarget(0)
	{
		vdbCOM.setZero();
	}
	FractureRB::FractureRB(btRigidBody* rb_, bool ownRB_)
		: rb(rb_), fractSim(NULL), mat(NULL), outDir(""), fragmentCounter(0),
		  updateCounter(0), haveMesh(false), haveParams(false), maxCracks(-1),
		  ownRB(ownRB_), haveBEM(false), remeshTarget(0)
	{
		vdbCOM.setZero();
		useEstSIFs=-1;
	}
	FractureRB::~FractureRB(){
		if( ownRB && rb  ){
			if( rb->getUserPointer() )    delete rb->getUserPointer();
			if( rb->getMotionState() )    delete rb->getMotionState();
			if( rb->getCollisionShape() ){
				delete rb->getCollisionShape();
				collInUse.shape=NULL; // if this shape was in use we've just deleted it
			}
			delete rb;
		}
		collInUse.deleteAll();
		collUpdate.deleteAll();
		if( fractSim ) delete fractSim;
		if( mat )      delete mat;
//		printf("done\n");
	}
	
	int FractureRB::initFractureSim(
		double crackMeshSize,
		double voxelSize,
		double young,
		double poisson,
		double density,
		double strength,
		double toughness,
		double compFactor,
		int    remesh,
		bool   ignoreCrackSelfIntersect,
		string matSpec,
		int maxCracks_,
		FRACTURE_METHOD method,
		double samplingDensity_
	){
		double t0=omp_get_wtime(), t1;
		fractSim = new FractureBEM(crackMeshSize,method);
		fractSim->setCPVersion(2);

		int check = getMeshFromRB();
		if( check==0 ){
			haveMesh=true;
			printf("\n%% ... initializing implicit surface");
			fractSim->initVDB(voxelSize,ignoreCrackSelfIntersect);
		}else
		if( check==-3 ){ // it is a btSphereShape
			// init VDB from levelsetsphere then remesh!
			check=initVDBfromSphere(voxelSize,ignoreCrackSelfIntersect);
			if( check==0 ) haveMesh=true; // slightly misusing the flag, should probably rename it to haveGeometry or something
			else printf("\n%% init from sphere returned %d",check);
		}else printf("\n%% !!! getMeshFromRB returned %d",check);

		remeshTarget = fractSim->getElems().size();
		//remesh from VDB if requested
		if( remesh > 3 && haveMesh ){ // remesh from level-set to requested number of elements
			printf("\n%% ... building BEM mesh from VDB surface");
			fractSim->remesh(remesh);
			remeshTarget = remesh;
		}
		printf("\n%% ... have %d elements in the BEM mesh",fractSim->getElems().size());
		mat = createMaterialModel(  matSpec,
			young, poisson, density, strength, toughness, compFactor
		);

		samplingDensity=samplingDensity_;
		fractSim->initFractureModel(*mat,false,samplingDensity); // requires initVDB before calling, also initializes nearTriGrid
		haveParams=true;

		fractSim->getRegions().clear();
		for(elem_map::iterator it = fractSim->getElems().begin();
			it != fractSim->getElems().end(); ++it
		){ // this assumes that there are only surface elements in the mesh -- might need to handle pre-existing cracks eventually
			fractSim->getRegions()[it->first]=it->first; // each element is treated as it's own region when placing boundary conditions
		}
		maxCracks=maxCracks_;

		if(haveMesh){
			precompMeshData();

			bndCnds.clear();
			for(id_map::iterator it = fractSim->getRegions().begin(); it != fractSim->getRegions().end(); ++it){
				stringstream bc;
				bc << it->second << "(0,0,0)"; // zero-traction BC, we'll fill in the values later ...
				bndCnds.push_back( bc.str() );
			}
			// we've delayed this until the first fracture simluation is run
			//printf("\n%% ... initializing BEM solver");
			//if( useEstSIFs >=0 && fractSim->getElems().size() > useEstSIFs ) fractSim->setMethodEstimated();
			//fractSim->initBEM(*mat,false,bndCnds);
			//haveBEM=true;

			if( haveParams ){
				// in order to get consistent masses we need to re-compute them at the beginning of the sim
				// we can't trust the the scene file will have reasonable values for masses!
				originalVolume=updateMass();
				
				stringstream outfile;
				outfile << outDir << *(string*)rb->getUserPointer();
				fractSim->writeVDB(outfile.str(),false);

				//debug output
				fractSim->getHiResCrackTip().setSIFsOutputFile(outfile.str()+"_sifs");
			}
		}
		t1=omp_get_wtime();
		printf("\n%% ... built breakable rigid body (%.4lfs)",t1-t0);
		return 0;
	}


	void FractureRB::clearContacts(){
		contactTractions.clear();
		contactDuration=0.0;
	}

	int  FractureRB::addContact(Eigen::Vector3d& p, Eigen::Vector3d& d, double impulse, double duration){
		unsigned int tri = fractSim->findClosestSurfaceTri(p);
		//printf("%% ... %s contact (%.3lg %.3lg %.3lg), tri %d\n",((string*)rb->getUserPointer())->c_str(),p[0],p[1],p[2],tri);
		return addContact(tri,d,impulse,duration);
	}

	int FractureRB::addContact(unsigned int tri, Eigen::Vector3d& d, double impulse, double duration){
		if( tri<ELEM_BASE_INDEX || tri >= fractSim->getElems().size() ) return -1;
		Eigen::Vector3d traction(d);
		traction.normalize(); // d should be normalized by caller, but just to be sure ...
		traction *= impulse/duration;
		traction /= elemArea[tri];

		if( contactTractions.count(tri)==0 ) contactTractions[tri].setZero();
		contactTractions[tri] += traction;
		
		if(duration > contactDuration) contactDuration = duration; // store the max. contact duration for the sim-run
		//printf("%% ... contact on tri %6d - traction is (%.3lg %.3lg %.3lg)\n",tri, traction[0],traction[1],traction[2]);
		return 0;
	}

	double FractureRB::getTotalContactForce(){
		int steps = (int)((contactDuration/fractSim->getTimeStep())+0.5); // round
		if( steps == 0 ) return 0.0;

		double totalForce=0.0;
		vect3d_map q;
		getBalancedContactTractions(q);
		for(vect3d_map::iterator it=q.begin(); it!=q.end(); ++it){
			totalForce += it->second.norm() * elemArea[it->first];
		}
		return totalForce;
	}

	int  FractureRB::runFractureSim(double maxTime, int rbTimeCode){
		double t0=omp_get_wtime(), t1, t2, tl; // runtime
		double t_max=maxTime, t_step; // time-steps
		int addedCracks=0, chk;
		if(contactDuration < t_max) t_max = contactDuration;
		t_step = fractSim->getTimeStep();
		int steps = (int)((t_max/t_step)+0.5); // round
		printf("\n%% fracture simulation for %s_%d ... %d time steps", ((string*)rb->getUserPointer())->c_str(), rbTimeCode, steps);
		if( steps == 0) printf("\n%% contact duration or max time too short (%.3lgs)",t_max);
		if( haveMesh && haveParams && steps > 0){

			if(!haveBEM){ // delayed initialization of BEM solver on first collision
				t1=omp_get_wtime();
				printf("\n%% ... initializing BEM solver");
				if( useEstSIFs >=0 && fractSim->getElems().size() > useEstSIFs ) fractSim->setMethodEstimated();
				try{
					fractSim->initBEM(*mat,false,bndCnds);
					haveBEM=true;
					t2=omp_get_wtime();
					printf("\n%%           BEM initialized ...\t%.4lfs",t2-t1); t1=t2;
				}catch(...){
					printf(" ERROR !!! ");
					haveParams=false; haveMesh=false;
				}
			}

			stringstream outfile;
			outfile << outDir << *(string*)rb->getUserPointer() << "_" << rbTimeCode;

			vect3d_map q;
			t1=omp_get_wtime();
			printf("\n%% setting new boundary data ...");
			getBalancedContactTractions(q);
			fractSim->updateBoundaryData(q);
			t2=omp_get_wtime(); printf("\t%.4lfs",t2-t1); t1=t2;

			int lastNumElems=fractSim->getElems().size(), numElems;
			int startElems=lastNumElems;
			
			bool done=false;
			for(int k=0; k<steps && !done; ++k){ // do fracture steps
				stringstream kstr; kstr << "_" << (k+1);
				printf("\n%% fracturing (%d/%d) ... ",k+1,steps);
				t1=omp_get_wtime(); tl=t1;
				if( useEstSIFs >=0 && lastNumElems > useEstSIFs ) fractSim->setMethodEstimated();
				printf("\n%% computing SIFs    ...");
				fractSim->computeBEM();
				if(k==0){ // only initial output
					t2=omp_get_wtime(); printf("\t%.4lfs",t2-t1); t1=t2;
					printf("\n%% writing output    ...");
					fractSim->writeMesh(outfile.str()+"_0");
				}
				t2=omp_get_wtime(); printf("\t%.4lfs",t2-t1); t1=t2;

				if( maxCracks < 0){
					chk=fractSim->seedCracksAndPropagate(); // unlimited seeding
					if( chk>0 ) addedCracks+=chk;
				}else if( maxCracks > /**/ addedCracks /*/fractSim->getCracks().size()/**/ ){
					chk=fractSim->seedCracksAndPropagate( maxCracks - /**/ addedCracks /*/fractSim->getCracks().size()/**/ ); // seed to limit
					if( chk>0 ) addedCracks+=chk;
				}else{
					fractSim->propagateCracks(); // no seeding
				}
				
				numElems=fractSim->getElems().size();
				if( numElems == lastNumElems ) done=true;
				printf("\t%.4lfs (have %d elements and %d cracks)",omp_get_wtime()-tl, numElems, fractSim->getCracks().size());
				lastNumElems = numElems;
				//printf("\n%% writing output ... ");
				//fractSim->writeMesh(outfile.str()+kstr.str());
				//t2=omp_get_wtime(); printf(" %.4lfs",t2-t1); t1=t2;
			}
			fractSim->writeVDB(outfile.str(),false);

			//// check net force and torque in the resulting traction field
			//vector_type& result_q = bem->getTractions();
			//Eigen::Vector3d force(0.0,0.0,0.0), torque(0.0,0.0,0.0), tmp;
			//for(elem_map::iterator it = bem->getElems().begin();
			//	it != bem->getElems().end(); ++it
			//){
			//	tmp[0] = result_q[(it->first-ELEM_BASE_INDEX)*3  ];
			//	tmp[1] = result_q[(it->first-ELEM_BASE_INDEX)*3+1];
			//	tmp[2] = result_q[(it->first-ELEM_BASE_INDEX)*3+2];
			//	force  += tmp*elemArea[it->first];
			//	torque += tmp.cross(elemCtr[it->first])*elemArea[it->first];
			//}
			//printf("%% ... BEM result force %.3g (%.3lg %.3lg %.3lg)\n%% ... BEM result torque %.3g (%.3lg %.3lg %.3lg)\n",
			//	force.norm(),force[0],force[1],force[2], torque.norm(),torque[0],torque[1],torque[2]);

			printf("\n%% fracture sim done in \t%.4lfs",omp_get_wtime()-t0);
			return numElems-startElems;
		}
		return -1;
	}


	// balance tractions --> goal is to have sum(forces)=0 and sum(torques)=0 over the surface
	// reads from member variable contactTractions, writes into output param q
	void FractureRB::getBalancedContactTractions(vect3d_map& q){
		q.clear();
		Eigen::Vector3d force(0.0,0.0,0.0), torque(0.0,0.0,0.0);

		// initialze q for all surface (i.e. non-crack) elements
		// also compute the net force and torque due to the contactTractions
		for(elem_map::iterator it = fractSim->getElems().begin();
			it != fractSim->getElems().end(); ++it
		)if(fractSim->getCracks().count(fractSim->getRegions()[it->first])==0){
			if( contactTractions.count(it->first) ){
				q[it->first]=contactTractions[it->first];
				force  += contactTractions[it->first]*elemArea[it->first];
				torque += contactTractions[it->first].cross(elemCtr[it->first])*elemArea[it->first];
			}else{
				q[it->first].setZero();
			}
		} //from now on iterate over q to get all elements that are not cracks
		
		// we want to find a traction field that is 
		// -- as close as possible to the contactTractions (in L2-norm)
		// -- has 0 net force and 
		// -- has 0 net torque (both integrated over the entire surface)
		// we find this by setting q = q_in - q_c (i.e. contactTractions - correction)
		// then do a constrained quadratic minimization on the correction
		// the function we want to minimize is the squared L2-norm of the correction q_c'*q_c
		// this results in a linear system which we solve via a Schur complement system
		// since our target function is so simple, and we have only 2 constraints
		// with 3 dimensions each (so 6 in all), we only need to invert a 6x6 matrix.
		
		// build the 6x6 Schur complement system S=B(A^-1)B' = -b, with b=[ force ; torque ]
		// actually A is an identity since we use the squared L2-norm so we ignore it from now on
		// but we could weight e.g. by element areas just as easily
		// as long a A is diagonal (with non-zero entries) we're good and S is symmetric
		matrix_type S(6,6); S.setZero(); double tmp;
		for(vect3d_map::iterator it = q.begin(); it != q.end(); ++it ){
			// run once through the elements and compute the entries of S
			// S(0,0) is (force-x, force-x), S(1,1) is (force-y, force-y) etc.
			// S(0,1)=S(0,2)=S(1,2) = 0 -- forces do not mix
			// S(0,3)=S(1,4)=S(2,5) = 0 -- forces never affect torques in the same axis
			// S(1,3) is (force-y, torque-x), S(2,3) is (force-z,torque-x) etc.
			// we'll build S block-wise S=[ Sff  Sft ; Sft', Stt ]
			// actually we'll first only build the upper triangle of S and then transpose-copy them down

			// ok let's start with the forces, Sff = Bf*Bf' where Bf consists of
			// blocks that are a*I_3x3 with a the element area and I the identity
			// consequently Sff = sum(a^2)*I_3x3 i.e. S(0,0)=S(1,1)=S(2,2)=sum(a^2) the sum of the element areas squared
			tmp = elemArea[it->first]; tmp*=tmp; // tmp is the squared area of the element
			S(0,0)+=tmp; S(1,1)+=tmp; S(2,2)+=tmp;
			
			// now do the torques, Stt = Bt*Bt' where Bt consists of blocks
			// that are a*M with M being the cross-product from the right operator
			// with the elements centroid: M*b = (b x elemCtr[i])
			// i.e. Stt=sum(a^2*M*M') with M=[0 mz -my ; -mz 0 mx ; my -mx 0], mx,my,mz being the components of elemCtr[i]
			// and  M*M' = [my^2+mz^2  -mx*my  -mx*mz ; -mx*my  mx^2+mz^2  -my*mz ; -mx*mz  -my*mz  mx^2+my^2]
			// first do the diagonal entries
			S(3,3) += tmp*( elemCtr[it->first][1]*elemCtr[it->first][1] + elemCtr[it->first][2]*elemCtr[it->first][2] );
			S(4,4) += tmp*( elemCtr[it->first][0]*elemCtr[it->first][0] + elemCtr[it->first][2]*elemCtr[it->first][2] );
			S(5,5) += tmp*( elemCtr[it->first][0]*elemCtr[it->first][0] + elemCtr[it->first][1]*elemCtr[it->first][1] );
			// now the upper triangle
			S(3,4) -= tmp*( elemCtr[it->first][0]*elemCtr[it->first][1] );
			S(3,5) -= tmp*( elemCtr[it->first][0]*elemCtr[it->first][2] );
			S(4,5) -= tmp*( elemCtr[it->first][1]*elemCtr[it->first][2] );
			
			// finally the force-torque interaction block Sft = Bf*Bt'
			// Sft = a*I*a*M' = a^2*M'
			// again we need only the upper triangle of this
			S(0,4) += tmp*elemCtr[it->first][2];
			S(0,5) -= tmp*elemCtr[it->first][1]; // note the minus!
			S(1,5) += tmp*elemCtr[it->first][0];
		}
		// fill in the lower triangle of Stt
		S(4,3)=S(3,4); S(5,3)=S(3,5); S(5,4)=S(4,5);
		// now the lower triangle of Sft (this is a skew-symmetric block!)
		S(1,3)=-S(0,4); S(2,3)=-S(0,5); S(2,4)=-S(1,5);
		// and finally the 3x3 lower block of Sft'
		S(3,1)=S(1,3); S(3,2)=S(2,3); S(4,2)=S(2,4);
		S(4,0)=S(0,4); S(5,0)=S(0,5); S(5,1)=S(1,5);
		
		vector_type v(6);
		v[0]= force[0]; v[1]= force[1]; v[2]= force[2]; // note: no minus here!
		v[3]=torque[0]; v[4]=torque[1]; v[5]=torque[2];
		S.llt().solveInPlace(v);
		//Eigen::LLT<matrix_type> lltSolver(S); // Mathematica could pre-solve this analytically ... but it's a terribly long formula
		//lltSolver.solveInPlace(v);
		
		// now that we have the Schur complement solution
		// the traction correction is given by A*q_c = -B'*v
		// (note that we've already omitted a minus above so we also ignore this one)
		// again we'll ignore A as it is the identity in the case of the squared L2-norm,
		// but obviously we could use any positive diagonal weight very easily
		
		// compute and apply q_c
		Eigen::Vector3d q_ce; // per element traction correction
		for(vect3d_map::iterator it = q.begin(); it != q.end(); ++it ){
			// as above, B consists of a 3x3 block Bf and another one Bt
			// so B'*v = [ Bf' Bt' ]*v = Bf'*v[0-2] + Bt'*v[3-5]
			// the result is the traction correction for this element

			// let's start with the force block, remember Bf has blocks that are a*I_3x3 per element
			q_ce[0] = elemArea[it->first] * v[0];
			q_ce[1] = elemArea[it->first] * v[1];
			q_ce[2] = elemArea[it->first] * v[2];

			// now we'll add the contribution of the torque block a*M
			// with M=[0 mz -my ; -mz 0 mx ; my -mx 0], mx,my,mz being the components of elemCtr[i]
			// this is skew-symmetric, so we just need to swap the signs to get the transpose
			q_ce[0] += elemArea[it->first] * (-elemCtr[it->first][2]*v[4] + elemCtr[it->first][1]*v[5]);
			q_ce[1] += elemArea[it->first] * ( elemCtr[it->first][2]*v[3] - elemCtr[it->first][0]*v[5]);
			q_ce[2] += elemArea[it->first] * (-elemCtr[it->first][1]*v[3] + elemCtr[it->first][0]*v[4]);

			it->second -= q_ce; // apply the correction
		}
		

		//// testing ...
		//printf("%% ... original force %.3g (%.3lg %.3lg %.3lg)\n%% ... original torque %.3g (%.3lg %.3lg %.3lg)\n",
		//	force.norm(),force[0],force[1],force[2], torque.norm(),torque[0],torque[1],torque[2]);
		//force.setZero(); torque.setZero();
		//for(vect3d_map::iterator it = q.begin(); it != q.end(); ++it){
		//	force  += it->second*elemArea[it->first];
		//	torque += it->second.cross(elemCtr[it->first])*elemArea[it->first];
		//}
		//printf("%% ... remaining force %.3g (%.3lg %.3lg %.3lg)\n%% ... remaining torque %.3g (%.3lg %.3lg %.3lg)\n",
		//	force.norm(),force[0],force[1],force[2], torque.norm(),torque[0],torque[1],torque[2]);


		// finally distribute remaining forces from collision elements, since these will receive Diri BCs
		// we're not doing this anymore and solve the Neumann problem with a regularizer instead
		//if(0) for(vect3d_map::iterator it_c = contactTractions.begin(); it_c != contactTractions.end(); ++it_c){				 
		//	//printf("%% ... ... remaining traction on collision (%.3lg %.3lg %.3lg)\n",q[it->first][0],q[it->first][1],q[it->first][2]);
		//	force = -q[it_c->first] / (bem->getElems().size()-contactTractions.size()); // use force vector as buffer -- it's still a traction
		//	for(vect3d_map::iterator it = q.begin(); it != q.end(); ++it){
		//		if( contactTractions.count(it->first)==0 ){ // this is not a collision element
		//			it->second += force;
		//		}
		//	}
		//	q[it_c->first].setZero();
		//}
	}

	int FractureRB::getMeshFromRB(){
		if(rb==NULL) return -1;
		return getMeshFromShape(rb->getCollisionShape(),fractSim->getNodes(),fractSim->getElems());
	}

	int FractureRB::getMeshFromShape(btCollisionShape* shape, node_map& nodes, elem_map& elems){
		if(dynamic_cast<btGImpactMeshShape*>(shape)){
			btGImpactMeshShape* mesh = dynamic_cast<btGImpactMeshShape*>(shape);
			return getMeshFromBtInterface(mesh->getMeshInterface(),nodes,elems);
		} else
		if(dynamic_cast<btBvhTriangleMeshShape*>(shape)){
			btBvhTriangleMeshShape* mesh = dynamic_cast<btBvhTriangleMeshShape*>(shape);
			return getMeshFromBtInterface(mesh->getMeshInterface(),nodes,elems);
		} else
		if(dynamic_cast<btBoxShape*>(shape)){
			btBoxShape* box = dynamic_cast<btBoxShape*>(shape);
			return getMeshFromBox(box, nodes, elems);
		} else
		if(dynamic_cast<btConvexHullShape*>(shape)){
			btConvexHullShape* convhull = dynamic_cast<btConvexHullShape*>(shape);
			return getMeshFromConvHull(convhull,nodes,elems);
		} else
		if(dynamic_cast<btSphereShape*>(shape)){
			collInUse.shape = shape; // using collInUse as temporary storage, will be cleared in initVDBfromSphere(...)
			return -3;
		} else
		if(dynamic_cast<btCompoundShape*>(shape)){
			btCompoundShape* cmpd = dynamic_cast<btCompoundShape*>(shape);
			if( cmpd->getNumChildShapes() >1 ){
				printf("\n%% !!! compound shapes with multiple children not implemented");
				return -4;
			}else return getMeshFromShape(cmpd->getChildShape(0),nodes,elems);
			//ToDo: Bullet can have offset coordinate origin (center of mass) for parts of compound shapes
			//      check that we're not messing stuff up here ...
			//      things sometimes go wrong if the "center of mass" is not set to (0,0,0) on mesh-shapes created by Maya
		}
		// if(dynamic_cast<btConcaveShape*>(shape))
			// printf("%% -- cast to btConcaveShape succeeded\n");
		// if(dynamic_cast<btConvexShape*>(shape))
			// printf("%% -- cast to btConvexShape succeeded\n");
		// if(dynamic_cast<btConvexInternalShape*>(shape))
			// printf("%% -- cast to btConvexInternalShape succeeded\n");
		// if(dynamic_cast<btConvexTriangleMeshShape*>(shape))
			// printf("%% -- cast to btConvexTriangleMeshShape succeeded\n");
		return -1; // default return (nothing found)
	}

	// read a mesh from a Bullet interface
	int FractureRB::getMeshFromBtInterface(btStridingMeshInterface* meshIntf, node_map& nodes, elem_map& elems){
		unsigned int nextNode = NODE_BASE_INDEX;
		unsigned int nextElem = ELEM_BASE_INDEX;
		nodes.clear(); elems.clear();

		for(int j=0; j<meshIntf->getNumSubParts(); ++j){
			void* vertbase; unsigned int* idxbase;
			int numverts,numtris,stride,idxstride;
			PHY_ScalarType type,idxtype;
			meshIntf->getLockedReadOnlyVertexIndexBase(
				(const unsigned char**)&vertbase,numverts,type,stride,
				(const unsigned char**)&idxbase,idxstride,numtris,idxtype,j
			);

			if( type==PHY_FLOAT ){
				float* verts = (float*) vertbase;
				int vertstride = stride/sizeof(float);
				for(int k=0; k<numverts; ++k){
					//printf("%.8lg %.8lg %.8lg\n",verts[vertstride*k  ],verts[vertstride*k+1],verts[vertstride*k+2]);
					nodes[nextNode].assign(3,0.0);
					nodes[nextNode][0] = verts[vertstride*k  ];
					nodes[nextNode][1] = verts[vertstride*k+1];
					nodes[nextNode][2] = verts[vertstride*k+2];
					++nextNode;
				}
			}else if( type==PHY_DOUBLE ){
				double* verts = (double*) vertbase;
				int vertstride = stride/sizeof(double);
				for(int k=0; k<numverts; ++k){
					//printf("%.8lg %.8lg %.8lg\n",verts[vertstride*k  ],verts[vertstride*k+1],verts[vertstride*k+2]);
					nodes[nextNode].assign(3,0.0);
					nodes[nextNode][0] = verts[vertstride*k  ];
					nodes[nextNode][1] = verts[vertstride*k+1];
					nodes[nextNode][2] = verts[vertstride*k+2];
					++nextNode;
				}
			}else printf("\n%% !!! unsupported vertex type %d in FractureRB::getMeshFromBtInterface", type);
			if( idxtype==PHY_INTEGER ){
				for(int k=0; k<numtris; ++k){
					//printf("%u %u %u\n",idxbase[(idxstride/sizeof(int))*k  ],idxbase[(idxstride/sizeof(int))*k+1],idxbase[(idxstride/sizeof(int))*k+2]);
					elems[nextElem].assign(3,0);
					elems[nextElem][0] = idxbase[(idxstride/sizeof(int))*k  ] + NODE_BASE_INDEX;
					elems[nextElem][1] = idxbase[(idxstride/sizeof(int))*k+1] + NODE_BASE_INDEX;
					elems[nextElem][2] = idxbase[(idxstride/sizeof(int))*k+2] + NODE_BASE_INDEX;
					++nextElem;
				}
			}else printf("\n%% !!! unsupported index type %d in FractureRB::getMeshFromBtInterface", idxtype);
			meshIntf->unLockReadOnlyVertexBase(j);
		}
		return 0;
	}

	// triangulate a convex hull using btConvexHullComputer
	int FractureRB::getMeshFromConvHull(btConvexHullShape* convhull, node_map& nodes, elem_map& elems){
		unsigned int nextNode = NODE_BASE_INDEX;
		unsigned int nextElem = ELEM_BASE_INDEX;
		nodes.clear(); elems.clear();

		btConvexHullComputer* ch = new btConvexHullComputer();
		double* points = new double[3*convhull->getNumPoints()];
		for(int i=0; i < convhull->getNumPoints(); ++i){
			points[3*i  ]=convhull->getUnscaledPoints()[i][0];
			points[3*i+1]=convhull->getUnscaledPoints()[i][1];
			points[3*i+2]=convhull->getUnscaledPoints()[i][2];
		}
		ch->compute(points,3*sizeof(double),convhull->getNumPoints(),0.0,0.0);
		for(int i=0; i < ch->vertices.size(); ++i){
			nodes[nextNode].assign(3,0.0);
			nodes[nextNode][0] = ch->vertices[i][0];
			nodes[nextNode][1] = ch->vertices[i][1];
			nodes[nextNode][2] = ch->vertices[i][2];
			++nextNode;
		}
		for(int i=0; i < ch->faces.size(); ++i){
			// faces contains an index into edges giving any one edge of a face
			// each edge has pointers to the next edge ccw around a face or cw around a vertex
			const btConvexHullComputer::Edge* currentEdge = &(ch->edges[ch->faces[i]]);
			int firstVertex = currentEdge->getSourceVertex();
			//printf("\n\n%% building triangulation for face %d of %d, source vertex is %d\n",i+1, ch.faces.size(),firstVertex );
			currentEdge = currentEdge->getNextEdgeOfFace(); // skip first edge (adjacent to first vertex)
			while( currentEdge->getTargetVertex() != firstVertex ){ // also skip last edge (also adjacent to first vertex)
				elems[nextElem].assign(3,0);
				elems[nextElem][0] = NODE_BASE_INDEX + firstVertex;
				elems[nextElem][1] = NODE_BASE_INDEX + currentEdge->getSourceVertex();
				elems[nextElem][2] = NODE_BASE_INDEX + currentEdge->getTargetVertex();
				++nextElem;
				currentEdge = currentEdge->getNextEdgeOfFace();
			}
		}
		delete ch;
		delete[] points;
		return 0;
	}

	// create a basic (very coarse!) box mesh
	int FractureRB::getMeshFromBox(btBoxShape* box, node_map& nodes, elem_map& elems){
		btVector3 hExt = box->getHalfExtentsWithoutMargin();
		elems.clear();
		nodes.clear();
		unsigned int e=ELEM_BASE_INDEX, n=NODE_BASE_INDEX;
		elems[e].assign(3,0); elems[e][0]=n  ; elems[e][1]=n+3; elems[e][2]=n+2; ++e;
		elems[e].assign(3,0); elems[e][0]=n+2; elems[e][1]=n+3; elems[e][2]=n+6; ++e;
		elems[e].assign(3,0); elems[e][0]=n+1; elems[e][1]=n+4; elems[e][2]=n+5; ++e;
		elems[e].assign(3,0); elems[e][0]=n+4; elems[e][1]=n+7; elems[e][2]=n+5; ++e;
		elems[e].assign(3,0); elems[e][0]=n  ; elems[e][1]=n+1; elems[e][2]=n+3; ++e;
		elems[e].assign(3,0); elems[e][0]=n+1; elems[e][1]=n+5; elems[e][2]=n+3; ++e;
		elems[e].assign(3,0); elems[e][0]=n+2; elems[e][1]=n+6; elems[e][2]=n+4; ++e;
		elems[e].assign(3,0); elems[e][0]=n+4; elems[e][1]=n+6; elems[e][2]=n+7; ++e;
		elems[e].assign(3,0); elems[e][0]=n  ; elems[e][1]=n+2; elems[e][2]=n+1; ++e;
		elems[e].assign(3,0); elems[e][0]=n+1; elems[e][1]=n+2; elems[e][2]=n+4; ++e;
		elems[e].assign(3,0); elems[e][0]=n+3; elems[e][1]=n+5; elems[e][2]=n+6; ++e;
		elems[e].assign(3,0); elems[e][0]=n+5; elems[e][1]=n+7; elems[e][2]=n+6; ++e;

		nodes[n].assign(3,0.0); nodes[n][0]=-hExt[0]; nodes[n][1]=-hExt[1]; nodes[n][2]=-hExt[2]; ++n;
		nodes[n].assign(3,0.0); nodes[n][0]= hExt[0]; nodes[n][1]=-hExt[1]; nodes[n][2]=-hExt[2]; ++n;
		nodes[n].assign(3,0.0); nodes[n][0]=-hExt[0]; nodes[n][1]= hExt[1]; nodes[n][2]=-hExt[2]; ++n;
		nodes[n].assign(3,0.0); nodes[n][0]=-hExt[0]; nodes[n][1]=-hExt[1]; nodes[n][2]= hExt[2]; ++n;
		nodes[n].assign(3,0.0); nodes[n][0]= hExt[0]; nodes[n][1]= hExt[1]; nodes[n][2]=-hExt[2]; ++n;
		nodes[n].assign(3,0.0); nodes[n][0]= hExt[0]; nodes[n][1]=-hExt[1]; nodes[n][2]= hExt[2]; ++n;
		nodes[n].assign(3,0.0); nodes[n][0]=-hExt[0]; nodes[n][1]= hExt[1]; nodes[n][2]= hExt[2]; ++n;
		nodes[n].assign(3,0.0); nodes[n][0]= hExt[0]; nodes[n][1]= hExt[1]; nodes[n][2]= hExt[2]; ++n;
		printf("\n%% ... created coarse tri-mesh from box shape - use remeshing!");
		return 0;
	}

	// compute element areas and centroids as well as nodal curvatures on surface elements
	void FractureRB::precompMeshData(){
		elemArea.clear();
		elemCtr.clear();
		nodeCurv.clear();

		elem_map& el = fractSim->getElems();
		node_map& co = fractSim->getNodes();
		id_set& cracks = fractSim->getCracks();
		id_map& regions=fractSim->getRegions();
		
		vect3d_map curvNormals; // temporary store the mean curvature normals
		double cot1,cot2,cot3, e1,e2,e3;

		for(elem_map::iterator it = el.begin(); it != el.end(); ++it ){
			if(cracks.count(regions[it->first])==0){
				Eigen::Vector3d
					a (co[it->second[0]][0], co[it->second[0]][1], co[it->second[0]][2]),
					b (co[it->second[1]][0], co[it->second[1]][1], co[it->second[1]][2]),
					c (co[it->second[2]][0], co[it->second[2]][1], co[it->second[2]][2]);
				elemCtr[it->first] = (a+b+c)/3.0;
				elemArea[it->first] = 0.5*(a-c).cross(b-c).norm();

				// to compute the mean curvature we need to first sum up the cotangent-weighted edge vectors
				// -Kn = 1/(4A) sum( (cot alpha + cot beta)(x-y) ) 
				// then scale by vertex-associated area and take the magnitude
				// we are going to associate 1/3 of a triangle's area with each of it's nodes
				// which should be reasonably accurate as stated in http://mrl.nyu.edu/~dzorin/papers/grinspun2006cds.pdf (Computing discrete shape operators on general meshes)
				// more accurate variants are in http://www.multires.caltech.edu/pubs/diffGeoOps.pdf (Discrete Differential-Geometry Operators for Triangulated 2-Manifolds)
				// the original cotangent formula uses the sum of adjacent triangle areas http://www.geometry.caltech.edu/pubs/DMSB_SIG99.pdf (Implicit Fairing of Irregular Meshes using Diffusion and Curvature Flow)
				
				// first make sure everything is properly initialized, as we'll use a lot of += soon
				if(nodeCurv.count(it->second[0])==0) nodeCurv[it->second[0]]=0.0;
				if(nodeCurv.count(it->second[1])==0) nodeCurv[it->second[1]]=0.0;
				if(nodeCurv.count(it->second[2])==0) nodeCurv[it->second[2]]=0.0;
				if(curvNormals.count(it->second[0])==0) curvNormals[it->second[0]].setZero();
				if(curvNormals.count(it->second[1])==0) curvNormals[it->second[1]].setZero();
				if(curvNormals.count(it->second[2])==0) curvNormals[it->second[2]].setZero();

				// we accumulate the vertex-associated area in the nodeCurv map
				// and do the same for the (unscaled) curvature normals in the curvNormals map
				nodeCurv[it->second[0]] += elemArea[it->first]/3.0;
				nodeCurv[it->second[1]] += elemArea[it->first]/3.0;
				nodeCurv[it->second[2]] += elemArea[it->first]/3.0;

				// edg1 is b-a, edg2 is c-b, edg3 is a-c
				// phi1 is btwn edg2 and edg3 and opposite edg1, phi2 is opposite edg2, phi3 is opposite edg3
				// we need to add edg1*cot(phi1)-edg3*cot(phi3) at a, edg2*cot(phi2)-edg1*cot(phi1) at b, and edg3*cot(phi3)-edg2*cot(phi2) at c
				// since phi1 = acos( dot(edg2,edg3) / ( norm(edg2)*norm(edg3) ) ), we'll use cot(acos(x)) = x/sqrt(1-x*x) (see wolfram alpha)
				e1 = (b-a).norm(); e2 = (c-b).norm(); e3 = (a-c).norm();
				cot1 = (c-b).dot(a-c)/(e2*e3); cot1 /= sqrt(1.0-cot1*cot1);
				cot2 = (b-a).dot(a-c)/(e1*e3); cot2 /= sqrt(1.0-cot2*cot2);
				cot3 = (b-a).dot(c-b)/(e1*e2); cot3 /= sqrt(1.0-cot3*cot3);
				curvNormals[it->second[0]] += ((b-a)*cot1 - (a-c)*cot3);
				curvNormals[it->second[1]] += ((c-b)*cot2 - (b-a)*cot1);
				curvNormals[it->second[2]] += ((a-c)*cot3 - (c-b)*cot2);

				//if(it == el.begin()){ // some debugging ...
				//	printf("\n***\nedge lengths: %lf, %lf, %lf\ncotangents: %lf, %lf, %lf\n",e1,e2,e3,cot1,cot2,cot3);
				//	printf(" cn1 ( %lf %lf %lf ), cn2 ( %lf %lf %lf ), cn3 ( %lf %lf %lf )\n",
				//		curvNormals[it->second[0]][0],curvNormals[it->second[0]][1],curvNormals[it->second[0]][2],
				//		curvNormals[it->second[1]][0],curvNormals[it->second[1]][1],curvNormals[it->second[1]][2],
				//		curvNormals[it->second[2]][0],curvNormals[it->second[2]][1],curvNormals[it->second[2]][2]);
				//}
			}
		}
		// finally we need to get the curvature normal's magnitude and divide by 4*area, remember that we stored the areas in the nodeCurv map
		for(vect3d_map::iterator it=curvNormals.begin(); it!=curvNormals.end(); ++it){
			//printf("\n*** node-area %lf, curvature normal norm %lf", nodeCurv[it->first],it->second.norm());
			nodeCurv[it->first] = it->second.norm() / (4.0*nodeCurv[it->first]);
			//printf(" = curvature %lf\n", nodeCurv[it->first]);
		}

		// testing: write VTK output with these data fields
        id_map em1; id_set em2; elem_map em3; state_map em4; value_map em5; // empty maps - somewhat abusing the PostProcessor class here ...
        PostProcessor pp(co, el, em1, em2, em3, em3, em4, em5);
        vector_type areaVector(el.size()), ctrVector(3*el.size()), curvVector(co.size());
		areaVector.setZero(); ctrVector.setZero(); curvVector.setZero();
        unsigned int idx=0;
		for(elem_map::iterator it = el.begin(); it != el.end(); ++it, ++idx){
			if(elemArea.count(it->first)){
				areaVector[idx]=elemArea[it->first];
				ctrVector[3*idx  ]=elemCtr[it->first][0];
				ctrVector[3*idx+1]=elemCtr[it->first][1];
				ctrVector[3*idx+2]=elemCtr[it->first][2];
			}
		}
		idx=0;
		for(node_map::iterator it = co.begin(); it != co.end(); ++it, ++idx){
			if(nodeCurv.count(it->first)) curvVector[idx]=nodeCurv[it->first];
		}
        pp.vtkAddData("area",1,true,areaVector);
		//pp.vtkAddData("centroid",3,true,ctrVector);
		pp.vtkAddData("curvature",1,false,curvVector);
		if( rb && rb->getUserPointer() )
			pp.vtkWrite(outDir+(*(string*)rb->getUserPointer())+"_precomp.vtk",VTK_TRIANGLE);
		//printf("el=[\n"); printMap(el); printf("];\n");
		//printf("co=[\n"); printMap(co); printf("]; x=co(:,1); y=co(:,2); z=co(:,3);\n");
	}
	
	double FractureRB::getCurvatureAndElementID(const Eigen::Vector3d& p, unsigned int& elem){
		if(!haveMesh) return 0.0;
		elem = fractSim->findClosestSurfaceTri(p);
		if( elem == ELEM_BASE_INDEX-1 ) return 0.0;
		vector<unsigned int>& nd = fractSim->getElems()[elem];
		//ToDo: could use barycentric projection of p to interpolate the curvature, but we probably don't need it to be so accurate
		return (nodeCurv[nd[0]] + nodeCurv[nd[1]] + nodeCurv[nd[2]])/3.0;
	}
	
	double FractureRB::getEovNuSq(){
		return mat->getE() / (1.0 - mat->getNu()*mat->getNu());
	}
}
