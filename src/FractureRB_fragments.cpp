#include "FractureRB.h"
#include "VDBWrapper.h"
#include "FractureBEM.h"
#include "MaterialModel.h"
#include "FractureRB_config.h"
#include "distPtToTri.h"

#include "vcgHelper.h"
#include "SubsampledCrackTip.h"
using namespace myVCG;

#include <omp.h> // just for omp_get_wtime

#include <btBulletDynamicsCommon.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>
//#include <GIMPACTUtils/btGImpactConvexDecompositionShape.h> // not quite working as intended

#include <openvdb/tools/LevelSetMeasure.h>
#include <openvdb/tools/LevelSetRebuild.h>
#include <openvdb/tools/LevelSetSphere.h>

using namespace std;

namespace FractureSim{
	// file-local defs
	btRigidBody* createConvexRigidBody(vdb::FloatGrid::Ptr grid, FractureRB& parent, std::string writeMeshFile, double gridVolume=-1.0); // creates a rigid body from the given grid (used for small fragments - no longer breakable)
	btRigidBody* createFragmentRB(btRigidBody* parent, btCollisionShape* shape, const btVector3& btCOM, btScalar mass);
	double getCenteredMeshFromVBD(VcgMesh* mesh, Eigen::Vector3d& com, vdb::FloatGrid::Ptr grid); // invokes VDB to mesh conversion and moves all vertices such that the center of mass is at the coordinate origin
	Eigen::Vector3d getCollMeshFromVDB(ColliderData& coll, vdb::FloatGrid::Ptr grid, std::string writeMeshFile, unsigned int targetTris=0); // computes a collision mesh used for updating the parent or creating large fragments
	int getTrisAndDistance(
		value_map& currentTriMinDist,  
		FractureBEM* source, FractureBEM* target, id_map newCrackIDs,
		vdb::FloatGrid::Ptr mask, vdb::Coord ijk
	);
	int getFrontEdgesFromMarking(id_map& markedElems, elem_map& sourceElems, edge_imap& edgs);
	void findPinchedNodes(const edge_imap& newCrackEdges, id_set& pinchedNodes);

	// actual implementations start here ...


	int FractureRB::splitFragments(set<FractureRB*>& largeFragments, set<btRigidBody*>& smallFragments){
		largeFragments.clear();
		smallFragments.clear();
		
		// 1. segment the level-set, get grids
		// 2. check size - can be almost original, large or small
		// 3. if we need to update original, do that
		// 4. mesh the grids, make RBs
		// 5. make breakable RBs for large ones
		
		vector<vdb::FloatGrid::Ptr> segments = fractSim->getLevelSet().getSegments(
			SEG_HDL_THR, false /*use tiles*/
		); // we don't need to use tiles as we don't allow interior cracks atm.
		printf("\n%% have %d fragments ... ", segments.size()); // a bit of debug output
//		printf("%% parent margin is %.3lg\n", rb->getCollisionShape()->getMargin());

		if(segments.size()<=1) return -1; // if there's only one segment we don't need to do anything, if there are no segments, something went wrong and we can't do anything ... should probably print a warning or something
		
		bool keepParent=false;
		double parentVolume = vdb::tools::levelSetVolume(*fractSim->getLevelSet().getObjectGrid());
		double segmentVolume; unsigned int ignoreCount=0;

		if( parentVolume < 0.0 ){
			printf("\n%% !!! negative parent volume !!!");
			parentVolume*=-1.0;
		}
		
		for(unsigned int i=0; i<segments.size(); ++i){
			//ToDo: write output files for creating the render-quality meshes of all the fragments (it's not very convenient with the current implementation)
			//      in particular avoid doing the segmentation twice

			segmentVolume = vdb::tools::levelSetVolume(*segments[i]);

			if( segmentVolume < 0.0 ){
				//printf("\n%% !!! negative segment volume (%.2lg) !!!", segmentVolume);
				segmentVolume*=-1.0;
			}

			stringstream name;
			name << *(string*) getRB()->getUserPointer();
			//printf(" vol %.1le  ", segmentVolume / originalVolume );

			if( segmentVolume <= FRAGMENT_IGNORE_THR* originalVolume ){
				// do nothing
				++ignoreCount;
			}else
			if( segmentVolume <= FRAGMENT_SMALL_THR* originalVolume){ //parentVolume ){
				// make a small fragment
				name << "_f" << ++fragmentCounter;
				// use segmentVolume for mass computation
				// we want all masses to be VDB based so we have consistent measurements
				btRigidBody* fragmentRB;
				if( segmentVolume <= FRAGMENT_SMALL_CONV* originalVolume ) fragmentRB = createConvexRigidBody(
					segments[i],
					*this,
					outDir+name.str(),
					segmentVolume
				);
				else{
					//printf("\n%% creating small meshed fragment\n");
					ColliderData collMesh;
					segments[i]=vdb::tools::levelSetRebuild(*segments[i],0.0,fractSim->getLevelSet().getBandwidth()); // make sure the narrow-band is properly resolved (also removes most interior cracks)
					segments[i]->setName("objectGrid");
					segments[i]->setGridClass(vdb::GRID_LEVEL_SET);

					Eigen::Vector3d tmp, com = getCollMeshFromVDB( collMesh, segments[i], outDir+name.str(), COLL_RES*remeshTarget ); // limit the number of triangles to 20* the surface elements in the parent's BEM mesh
					tmp = com - vdbCOM;
					fragmentRB = createFragmentRB( rb, collMesh.shape, btVector3(tmp[0],tmp[1],tmp[2]), segmentVolume*mat->getDensity() );
					ColliderData::storeData(collMesh);
//					printf("%% meshed ");
				}
				if( fragmentRB ){
					fragmentRB->setUserPointer( new string( name.str() ) );
					fragmentRB->getCollisionShape()->setMargin( rb->getCollisionShape()->getMargin() );
					smallFragments.insert( fragmentRB );
					//printf("\n%% !!! small fragment has rest.coeff. %.3lf (%s)\n",fragmentRB->getRestitution(),name.str().c_str());
//					printf("%% small fragment margin is %.3lg\n", fragmentRB->getCollisionShape()->getMargin());
				}
			}else
			if( segmentVolume >= FRAGMENT_ORI_THR* parentVolume){// originalVolume){// replace the implicit surface of the parent
				name << "_" << ++updateCounter;
				keepParent=true;
				segments[i]=vdb::tools::levelSetRebuild(*segments[i],0.0,fractSim->getLevelSet().getBandwidth()); // make sure the narrow-band is properly resolved (also removes most interior cracks)
				segments[i]->setName("objectGrid");
				segments[i]->setGridClass(vdb::GRID_LEVEL_SET);
				fractSim->getLevelSet().setObjectGrid( segments[i] );
				
				printf("\n%% collision shape update on parent");
				//queue the update of the collision shape to a mesh of the remaining part ...
				if( pendingCollisionUpdate() ) collUpdate.deleteAll(); // we may have skipped an older update
				updateCOM = getCollMeshFromVDB( collUpdate, segments[i], outDir+name.str(), COLL_RES*remeshTarget );
				updateVolume = segmentVolume;
				collUpdate.shape->setMargin( rb->getCollisionShape()->getMargin() );
				//printf("\n%% collision shape update pending ... ");

			}else{ // make a large fragment
				// ... get collision mesh
				ColliderData collMesh;
				name << "_F" << ++fragmentCounter;
				printf("\n%% building breakable fragment %s",name.str().c_str());

				segments[i]=vdb::tools::levelSetRebuild(*segments[i],0.0,fractSim->getLevelSet().getBandwidth()); // make sure the narrow-band is properly resolved (also removes most interior cracks)
				segments[i]->setName("objectGrid");
				segments[i]->setGridClass(vdb::GRID_LEVEL_SET);

				Eigen::Vector3d tmp, com = getCollMeshFromVDB( collMesh, segments[i], outDir+name.str(), COLL_RES*remeshTarget );
				// ... make rigid body
				tmp = com - vdbCOM;
				btRigidBody* fragmentRB = createFragmentRB( rb, collMesh.shape, btVector3(tmp[0],tmp[1],tmp[2]), segmentVolume*mat->getDensity() );

				if( fragmentRB ){
					fragmentRB->setUserPointer( new string( name.str() ) ); // need to set this before calling initFractureSim on the FractureRB object
					fragmentRB->getCollisionShape()->setMargin( rb->getCollisionShape()->getMargin() );
					// ... create breakable rigid body from VDB
					FractureRB* breakableFragment = new FractureRB( fragmentRB );
					breakableFragment->setOutputDir( outDir );
					breakableFragment->setEstSIFsThreshold( useEstSIFs );

					breakableFragment->initFractureSim(
						collMesh,
						segments[i],
						com,
						fractSim->getTargetMeshSize(),
						fractSim->getLevelSet().getVoxelSize(),
						*mat,
						remeshTarget, // ToDo? could use a fraction of this depending on the volume
						originalVolume,
						fractSim,
						maxCracks,
						samplingDensity
					);

					largeFragments.insert( breakableFragment );
//					printf("%% large fragment margin is %.3lg\n", fragmentRB->getCollisionShape()->getMargin());
				}
			}
		}
		printf("\n%% ignored %u tiny fragments", ignoreCount);
		return keepParent?updateCounter:0;
	}
	
	void FractureRB::doCollisionUpdate(){
		if( pendingCollisionUpdate() ){
			//printf("\n%% processing collision shape update ... ");
			if( collInUse.shape ) collInUse.deleteAll();
			collInUse.set(collUpdate);
			collUpdate.setNull();

			rb->setCollisionShape( collInUse.shape );
			// update mass and inertia ...
			double mass = updateVolume*mat->getDensity();
			btVector3 inertia;
			collInUse.shape->calculateLocalInertia(mass, inertia);
			rb->setMassProps(mass, inertia);

			// correct for the change in COM between the old and new collision shapes
			Eigen::Vector3d com = updateCOM - vdbCOM;
			btVector3 btCOM(com[0],com[1],com[2]); // assuming the COM was at the origin in the original mesh, this is the COM shift
			rb->getWorldTransform().setOrigin(
				rb->getWorldTransform().getOrigin() + rb->getWorldTransform().getBasis()*btCOM
			);
			vdbCOM=updateCOM;
			//printf("\n%% collision shape update done ... ");
			//printf("\n%% !!! new rest.coeff. is %.3lf (%s)\n",rb->getRestitution(), ((string*)rb->getUserPointer())->c_str());
//			printf("%% new margin is %.3lg\n", rb->getCollisionShape()->getMargin());
		}
	}
	
	Eigen::Vector3d getCollMeshFromVDB(ColliderData& coll, vdb::FloatGrid::Ptr grid, string writeMeshFile, unsigned int targetTris){
		VcgMesh mesh;
		Eigen::Vector3d com;
		getCenteredMeshFromVBD(&mesh, com, grid);
		
		// simplify mesh if targetTris > 0
		if( targetTris > 0 ){
            TriEdgeCollapseQuadricParameter qparams;
            qparams.PreserveTopology = true;
			qparams.QualityQuadric   = true;
			qparams.QualityThr       = DBL_MAX;			
            vcg::LocalOptimization<VcgMesh> DeciSession(mesh,&qparams);
            DeciSession.Init<MyTriEdgeCollapse>();
            DeciSession.SetTargetSimplices(targetTris);
            DeciSession.DoOptimization();
		}

		// since VCG uses a complicated structure of object containers,
		// we copy vertex and triangle arrays to a Bullet-friendly format
		int nTris=mesh.FN(), nVerts=mesh.VN(),j,k;
		coll.tridata = new char[3*sizeof(int)*nTris];
		coll.vertdata= new char[3*sizeof(btScalar)*nVerts];
		int* tris = (int*)(coll.tridata);
		btScalar* verts = (btScalar*)(coll.vertdata);
		id_map vertexId;
		
		j=-1; k=-1;
		for(int i=0; i<mesh.vert.size(); ++i) if(!mesh.vert[i].IsD()){
			verts[++j]=mesh.vert[i].P()[0];
			verts[++j]=mesh.vert[i].P()[1];
			verts[++j]=mesh.vert[i].P()[2];
			vertexId[i]=++k;
		}
		j=-1;
		for(int i=0; i<mesh.face.size(); ++i) if(!mesh.face[i].IsD()){
			tris[++j] = vertexId[ vcg::tri::Index(mesh, mesh.face[i].V(0)) ];
			tris[++j] = vertexId[ vcg::tri::Index(mesh, mesh.face[i].V(1)) ];
			tris[++j] = vertexId[ vcg::tri::Index(mesh, mesh.face[i].V(2)) ];
		}

		if( !writeMeshFile.empty() ){
			vcg::tri::io::ExporterOBJ<VcgMesh> writerObj;
			writerObj.Save(mesh,writeMeshFile.append(".obj").c_str(), vcg::tri::io::Mask::IOM_NONE);
		}
		mesh.Clear();
		
		coll.iface = new btTriangleIndexVertexArray(
		  nTris, tris, 3*sizeof(int)/*stride*/, nVerts, verts, 3*sizeof(btScalar)/*stride*/
		);
		// directly using btBvhTriangleMeshShape will NOT WORK for mesh-mesh collisions!
		// could probably use a convex decomposition (this functionality is allegedly part of the Bullet lib)
		//coll.shape = new btBvhTriangleMeshShape(coll.iface,true,true);
		/**/
		btGImpactMeshShape* gimesh = new btGImpactMeshShape(coll.iface);
		gimesh->updateBound();
		coll.shape = gimesh;
		/*/
		// this is quite slow and not actually working
		printf("\n%% ... building convex decomposition of collider mesh ...");
		btVector3 scale(1.0,1.0,1.0);
		btGImpactConvexDecompositionShape* cd = new btGImpactConvexDecompositionShape(coll.iface,scale);
		coll.shape = cd;
		printf(" done");
		/**/
		return com;
	}

	btRigidBody* createConvexRigidBody(vdb::FloatGrid::Ptr grid, FractureRB& parent, string writeMeshFile, double gridVolume){
		VcgMesh mesh;
		Eigen::Vector3d com;
		double volume = getCenteredMeshFromVBD(&mesh, com, grid);

		if( !writeMeshFile.empty() ){
			vcg::tri::io::ExporterOBJ<VcgMesh> writerObj;
			writerObj.Save(mesh,writeMeshFile.append(".obj").c_str(), vcg::tri::io::Mask::IOM_NONE);
		}

		// this method is supposed to be used for small fragments
		// a convex hull collision shape should be accurate enough ...
		btConvexHullShape* shape = new btConvexHullShape();
		for(int i=0; i<mesh.vert.size(); ++i) if(!mesh.vert[i].IsD()){
			shape->addPoint( btVector3( mesh.vert[i].P()[0], mesh.vert[i].P()[1], mesh.vert[i].P()[2] ), false ); // don't update AABB of shape when inserting points
		}
		mesh.Clear();
		shape->recalcLocalAabb();

		com-=parent.getCOMshift();
		btRigidBody* rb = createFragmentRB(
			parent.getRB(),
			shape,
			btVector3(com[0],com[1],com[2]),
			(btScalar)( (gridVolume>0.0?gridVolume:volume)*parent.getMaterial()->getDensity() )
		);
		return rb;
	}

	btRigidBody* createFragmentRB(btRigidBody* parent, btCollisionShape* shape, const btVector3& btCOM, btScalar mass){
		btVector3 inertia;
		shape->calculateLocalInertia(mass, inertia);

		// compute the position of the fragment in world space
		// based on the shift in COM from the parent's material space to it's own
		// as well as the parent's rotation
		btVector3 position = parent->getWorldTransform().getOrigin();
		btQuaternion rotation = parent->getWorldTransform().getRotation();
		position+=quatRotate(rotation,btCOM);
		
		btRigidBody::btRigidBodyConstructionInfo rbci(
			mass,
			new btDefaultMotionState( btTransform(rotation, position)),
			shape,
			inertia
		);
//		printf("%% ... new fragment has volume %.3le, inertia (%.1le,%.1le,%.1le), com (%.3lf %.3lf %.3lf)\n",volume,inertia[0],inertia[1],inertia[2],com[0],com[1],com[2]);
		btRigidBody* rb = new btRigidBody(rbci);
		rb->setRestitution(    parent->getRestitution()    );
		rb->setFriction(       parent->getFriction()       );
		rb->setLinearVelocity( parent->getVelocityInLocalPoint(btCOM) );
		rb->setAngularVelocity(parent->getAngularVelocity());
		return rb;
	}
	
	double getCenteredMeshFromVBD(VcgMesh* mesh, Eigen::Vector3d& com, vdb::FloatGrid::Ptr grid){
		vdb::tools::VolumeToMesh* mesher =
            new vdb::tools::VolumeToMesh(0.0,1.0); // first param is isosurface (which level-set), second is adaptivity of mesh in [0,1]
		(*mesher)(*grid);
//		//printf("%% meshed\n");
		double volume=0.0, vol;
		Eigen::Vector3d a,b,c;
		com.setZero();
		copyVDBtoVCG(*mesher,*mesh);
		delete mesher;
		
		for(int i=0; i<mesh->face.size(); ++i) if(!mesh->face[i].IsD()){
			mesh->face[i].V(0)->P().ToEigenVector(a);
			mesh->face[i].V(1)->P().ToEigenVector(b);
			mesh->face[i].V(2)->P().ToEigenVector(c);
			vol = a.dot( b.cross( c )) / 6.0;
			volume += vol;
			com += vol*( a+b+c )/4.0; // volume weighted centroid of tet (a,b,c,0)
		}
		com/=volume;
		for(int i=0; i<mesh->vert.size(); ++i) if(!mesh->vert[i].IsD()){
			mesh->vert[i].P().ToEigenVector(a);
			a-=com;
			mesh->vert[i].P().FromEigenVector(a);
		}
		return volume;
	}

	double FractureRB::updateMass(){
		double volume = vdb::tools::levelSetVolume(*fractSim->getLevelSet().getObjectGrid());
		double mass = volume*mat->getDensity();
		btVector3 inertia;
		rb->getCollisionShape()->calculateLocalInertia(mass, inertia);
		rb->setMassProps(mass, inertia);
		return volume;
	}

	int FractureRB::initVDBfromSphere(double voxelSize, bool ignoreCrackSelfIntersect){
		btSphereShape* sphere = dynamic_cast<btSphereShape*>(collInUse.shape);
		collInUse.setNull();
		if( sphere==NULL ) return -1;
		fractSim->getNodes().clear();
		fractSim->getElems().clear();
		fractSim->initVDB(voxelSize,ignoreCrackSelfIntersect);
		fractSim->getLevelSet().setObjectGrid(
			vdb::tools::createLevelSetSphere<vdb::FloatGrid>(
				sphere->getRadius(), vdb::Vec3f::zero(), voxelSize
			)
		);
		fractSim->getLevelSet().getObjectGrid()->setName("objectGrid");
		printf("\n%% ... built VDB sphere with radius %.3lg and voxel size %.3lg",sphere->getRadius(),voxelSize);
		return 0;
	}

	//***********************************************************
	//** stuff below is for creating a breakable fragment ...  **
	//***********************************************************

	template<typename GridPtrType>
	int FractureRB::initFractureSim(
		ColliderData& collMesh,
		GridPtrType vdbImplicitSurface,
		const Eigen::Vector3d& com,
		double crackMeshSize,
		double voxelSize,
		const MaterialModel& matMdl,
		int    remesh,
		double originalVolume_,
		FractureBEM* parent,
		int maxCracks_,
		double samplingDensity_
	){
		double t0=omp_get_wtime(), t1;
		FRACTURE_METHOD mtd = FULL_BEM; // if we switch based on (number of elems > useEstSIFs) start with FULL_BEM (we'll check once we have the mesh, just before initializing the BEM)
		if( useEstSIFs < 0 ) mtd = parent->getMethod(); 
		fractSim = new FractureBEM(crackMeshSize,mtd);
		fractSim->setCPVersion(2);
		
		collInUse.set(collMesh);
		vdbCOM=com; // store center of mass shift
		
		printf("\n%% ... initializing implicit surface");
		if( fractSim->initVDB(voxelSize,parent->getLevelSet().ignoreCrackSelfIntersection())==0 ){
			fractSim->getLevelSet().setObjectGrid(vdbImplicitSurface);
			parent->getLevelSet().getCrackIDs( fractSim->getCracks() ); // fetch crack-IDs from parent (to be copied)
		}

		printf("\n%% ... building BEM mesh");
		remeshTarget=remesh;
		if( fractSim->remesh(remesh) > 3 ) haveMesh=true; // meshes with 3 or less elements don't make any sense
		printf("\n%% ... have %d elements in the BEM mesh",fractSim->getElems().size());

		mat = matMdl.clone();
		samplingDensity=samplingDensity_;
		fractSim->initFractureModel(*mat,false,samplingDensity); // requires initVDB before calling, also initializes nearTriGrid
		haveParams=true;

		// re-assign region IDs for cracks to be copied
		unsigned int nextCrackID=fractSim->getElems().size()+ELEM_BASE_INDEX;
		id_map newCrackIDs;
		for(id_set::iterator it=fractSim->getCracks().begin(); it!=fractSim->getCracks().end();++it){
			newCrackIDs[*it]=nextCrackID++; //post-inc!
		}
		// copy crack level-sets from parent and apply re-numbering, then update the crack-map
		parent->getLevelSet().copyCrackGrids( fractSim->getLevelSet(), newCrackIDs, fractSim->getLevelSet().getObjectGrid() );
		fractSim->getLevelSet().getCrackIDs( fractSim->getCracks() ); // write the IDs of copied cracks to the crack-map
		
		// prepare regions for boundary conditions on the surface
		for(elem_map::iterator it = fractSim->getElems().begin();
			it != fractSim->getElems().end(); ++it
		){  // this assumes that there are no crack regions with IDs less or equal to the highest element ID
			// copied crack have been renumbered to higher IDs
			fractSim->getRegions()[it->first]=it->first; // each surface element is treated as it's own region when placing boundary conditions
		}
		maxCracks=maxCracks_;

		if(haveMesh){
			precompMeshData();

			// add BEM mesh for copied cracks (if inside)
			copyInteriorCrackMesh( parent, fractSim, newCrackIDs );

			// build boundary conditions (types only, not values) for the surface
			bndCnds.clear();
			for(id_map::iterator it = fractSim->getRegions().begin(); it != fractSim->getRegions().end(); ++it){
				if( fractSim->getCracks().count( it->second ) ==0){ // this is not a (copied) crack region
					stringstream bc;
					bc << it->second << "(0,0,0)"; // zero-traction BC, we'll fill in the values later ...
					bndCnds.push_back( bc.str() );
				}
			}
			// add boundary conditions for copied cracks
			for(id_set::iterator it=fractSim->getCracks().begin(); it!=fractSim->getCracks().end(); ++it){
				stringstream bc;
				bc << (*it) << "(crack)";
				bndCnds.push_back( bc.str() );
			}

			// we've delayed this until the first fracture simluation is run
			//printf("\n%% ... initializing BEM solver");
			//if( useEstSIFs >=0 && fractSim->getElems().size() > useEstSIFs ) fractSim->setMethodEstimated();
			//try{
			//	fractSim->initBEM(*mat,false,bndCnds);
			//}catch(...){
			//	printf(" ERROR !!! ");
			//	haveParams=false; haveMesh=false;
			//}
			
			if( haveParams ){
				stringstream outfile;
				outfile << outDir << *(string*)rb->getUserPointer();
				fractSim->writeVDB(outfile.str(),false);

				//printf("\n%% writing BEM mesh for debugging ... ");
				//fractSim->computeBEM(); //just for debug output
				//fractSim->writeMesh(outfile.str() + "_cp");
				//fractSim->writeCrackTip(outfile.str() + "_ct.vtk");

				//more debug output
				fractSim->getHiResCrackTip().setSIFsOutputFile(outfile.str()+"_sifs");
			}
		}
		t1=omp_get_wtime();
		originalVolume = originalVolume_;
		printf("\n%% ... mass is %.3le", 1.0/rb->getInvMass());
		printf("\n%% built breakable rigid body (%.4lfs)",t1-t0);
		return 0;
	}


	// copy elements and nodes with associated data
	// then re-build the crack-front edges and re-initialize the high-res crack-front
	// finally, copy high-res crack-front data from the source if it exists
	int FractureRB::copyInteriorCrackMesh(FractureBEM* source, FractureBEM* target, id_map newCrackIDs ){
		if( target->getNodes().empty() ||
			target->getElems().empty() ){
			printf("\n!!! cannot copy crack meshes into empty target, write surface mesh first!");
			return -1;
		}
	
		// strategy:
		// - build high-res mask on a VDB grid
		// - run through all active voxels
		// -- find containing voxel in source nearTriGrid
		// -- mark all crack-elements listed there
		// -- deactivate mask-voxels in bounding box of nearTri-voxel

		// for the surface, we need to update nodes, elems, and regions
		// but we already have cracks in the target (match to newCrackIDs)
		// we'll also keep entries in fracturedNodes/fracturedElems if they are marked
		// for the crack-front we need to update crackTips, crackTipParents and crackTipStates
		// but crackTipFaceNormals and crackTipTangents are re-built by the post-processor
		
		id_map markedNodes, markedElems; // a node/element is considered "marked" if present in the key-set, we use the value to assign new IDs later
		
		vdb::FloatGrid::Ptr mask = source->getLevelSet().maskCracks(target->getLevelSet().getObjectGrid());
		vdb::FloatGrid::Accessor acc=mask->getAccessor();
		value_map currentTriMinDist;
		while( mask->activeVoxelCount() > 0 ){
			vdb::Coord ijk = mask->beginValueOn().getCoord();
			if( acc.getValue(ijk) >= /**0.0/*/ -0.1*halfDiag*mask->voxelSize()[0] /**/  ){
				acc.setValueOff(ijk);
			}else{
				getTrisAndDistance( // also deactivates voxels that have the same bounding box in the source's near-tri-grid
					currentTriMinDist,  
					source, target, newCrackIDs, mask, ijk
				);

				for(value_map::iterator el_it=currentTriMinDist.begin();
					el_it!=currentTriMinDist.end(); ++el_it
				){
					//printf("\n tri %d has dist %.3lg ", el_it->first, sqrt(el_it->second));
					if(el_it->second <= mask->voxelSize()[0]*mask->voxelSize()[0]){
						markedElems[el_it->first]=0; // just put it in the map, we're not using the mapped value yet
						//printf("**");
					}
				}
			}
		}
		
		//ToDo: move stuff below (where we don't need access to VDB methods) into the FractureBEM class
		
		//compute crack-front edges (and node-set) before node-marking and copying
		//      -- un-mark all triangles that have no DOFs (i.e. tris that have 3 crack-front nodes) --> this will NOT affect the free nodes!
		//      -- update crack-front edges, check that we have no "pinches" ie. each crack-front node can only be part of exactly 2 edges
		//      -- resolve these cases by marking an additional triangle that contains the problem-node (and a free node so we don't go backwards!)
		edge_imap newCrackEdges;
		id_set    newCrackNodes;
		getFrontEdgesFromMarking(markedElems, source->getElems(), newCrackEdges);
		newCrackNodes = nodeSet(newCrackEdges);
		id_set unmarkElems;
		// check if each marked tri contains a free (ie. non-crack-front) node, otherwise un-mark it
		for(id_map::iterator el_it = markedElems.begin(); el_it != markedElems.end(); ++el_it){
			bool haveFreeNode=false;
			for(elem_map::mapped_type::iterator nd_it = source->getElems()[el_it->first].begin(); !haveFreeNode && (nd_it != source->getElems()[el_it->first].end()); ++nd_it){
				if( newCrackNodes.count( *nd_it ) == 0 ) haveFreeNode=true;
			}
			if(!haveFreeNode) unmarkElems.insert(el_it->first); // delayed un-mark
		}
		for(id_set::iterator it=unmarkElems.begin(); it!=unmarkElems.end(); ++it){
			markedElems.erase(*it);
		}
		// first part done, no more elements without DOFs marked, now update crack-front edges
		// check for "pinched" nodes
		id_set pinchedNodes;
		getFrontEdgesFromMarking(markedElems, source->getElems(), newCrackEdges);
		findPinchedNodes(newCrackEdges, pinchedNodes);
		while(!pinchedNodes.empty()){
//			printf("\n%%!!! WARNING: HAVE %d PINCHED CRACK-FRONT NODES !!!\n", pinchedNodes.size());
			// if we found a pinched node, we mark its whole 1-ring (all elements containing the node)
			unsigned int pNd = *(pinchedNodes.begin());
			for(elem_map::iterator el_it=source->getElems().begin(); el_it!=source->getElems().end(); ++el_it){
				for(elem_map::mapped_type::iterator nd_it=el_it->second.begin(); nd_it!=el_it->second.end(); ++nd_it){
					if( pNd==(*nd_it) ){
						markedElems[el_it->first]=0; // target ID will be assigned later
//						printf("%% --> marked element %d\n", el_it->first);
					}
				}
			}
			
			getFrontEdgesFromMarking(markedElems, source->getElems(), newCrackEdges);
			findPinchedNodes(newCrackEdges, pinchedNodes);
		}
		
		// build node-marking for marked elements
		for(id_map::iterator el_it = markedElems.begin(); el_it != markedElems.end(); ++el_it){
			for(elem_map::mapped_type::iterator nd_it = source->getElems()[el_it->first].begin(); nd_it != source->getElems()[el_it->first].end(); ++nd_it){
				markedNodes[*nd_it]=0; // just put it in the map, we're not using the mapped value yet
			}
		}
		//// debug output
		//printf("\n nodes marked for copying:\n");
		//for(id_map::iterator nd_it = markedNodes.begin(); nd_it != markedNodes.end(); ++nd_it){
		//	std::vector<double>& nd = source->getNodes()[nd_it->first];
		//	printf("   %d => %d%c (inside %d) at (%.3lf %.3lf %.3lf)\n", nd_it->first, nd_it->second, extraNodes.count(nd_it->first)?'*':' ', target->getLevelSet().isInside( nd[0], nd[1], nd[2] ), nd[0], nd[1], nd[2]);
		//}
		// assign new elem IDs and node IDs to all marked elems and nodes
		unsigned int nextNode=(target->getNodes().rbegin()->first)+1;
		unsigned int nextElem=(target->getElems().rbegin()->first)+1;
		for(id_map::iterator el_it = markedElems.begin(); el_it != markedElems.end(); ++el_it){
			el_it->second = nextElem++; //post-inc!
//			printf("\n%% el %d --> %d",el_it->first,el_it->second);
			for(elem_map::mapped_type::iterator nd_it = source->getElems()[el_it->first].begin(); nd_it != source->getElems()[el_it->first].end(); ++nd_it){
				if( markedNodes.count(*nd_it) >0 && markedNodes[*nd_it]==0 )
					markedNodes[*nd_it]=nextNode++; //post-inc!
				// note that HyENA requires the nodes to be numbered
				// in order of appearance in the element list
			}
		}
		//printf("\n%% new mesh data:\n");
		// copy node coordinates
		for(id_map::iterator nd_it = markedNodes.begin(); nd_it != markedNodes.end(); ++nd_it){
			target->getNodes()[nd_it->second] = source->getNodes()[nd_it->first];
			//printf("%% nd %u: (%.3lf %.3lf %.3lf)\n",nd_it->second,
			//	target->getNodes()[nd_it->second][0],target->getNodes()[nd_it->second][1],target->getNodes()[nd_it->second][2]
			//);
		}
		// copy element data
		for(id_map::iterator el_it = markedElems.begin(); el_it != markedElems.end(); ++el_it){
			target->getElems()[el_it->second].assign(source->getElems()[el_it->first].size(),0);
			for(int k=0; k<source->getElems()[el_it->first].size(); ++k){
				target->getElems()[el_it->second][k]=markedNodes[source->getElems()[el_it->first][k]];
			}
			target->getRegions()[el_it->second]=newCrackIDs[source->getRegions()[el_it->first]];
			//printf("%% el %u: (%u %u %u) r %u\n", el_it->second,
			//	target->getElems()[el_it->second][0],target->getElems()[el_it->second][1],target->getElems()[el_it->second][2],
			//	target->getRegions()[el_it->second]
			//);
		}
		// copy seeding data for branches (usually only one of these two is in use, but we don't want to assume which one here)
		// only copy values which correspond to MARKED nodes or elements!
		for(id_set::iterator nd_it = source->getFracturedNodes().begin(); nd_it != source->getFracturedNodes().end(); ++nd_it){
			if( markedNodes.count( *nd_it ) ) target->getFracturedNodes().insert( markedNodes[*nd_it] );
		}
		for(id_set::iterator el_it = source->getFracturedElems().begin(); el_it != source->getFracturedElems().end(); ++el_it){
			if( markedElems.count( *el_it ) ) target->getFracturedElems().insert( markedElems[*el_it] );
		}

		// since we've (probably) truncated the mesh, we'll first copy the
		// new crack-front edges to the target mesh, then we'll check if 
		// some of them exist in the source mesh and copy data where we have it		
		target->getCrackTips().clear();
		target->getCrackTipParents().clear();
		target->getCrackTipStates().clear();
		unsigned int nextCrackTip = ELEM_BASE_INDEX;
//		printf("\n%% new crack-front edges:\n");
		for(edge_imap::iterator ed_it = newCrackEdges.begin(); ed_it != newCrackEdges.end(); ++ed_it){
//			printf("%% (%d, %d) in %d\n",
//				markedNodes[ed_it->first.first ],
//				markedNodes[ed_it->first.second],
//				markedElems[ed_it->second]
//			);  // debug output
			target->getCrackTips()[nextCrackTip].assign(2,0);
			target->getCrackTips()[nextCrackTip][0]=markedNodes[ed_it->first.first ];
			target->getCrackTips()[nextCrackTip][1]=markedNodes[ed_it->first.second];
			target->getCrackTipParents()[nextCrackTip].assign(2,0);
			target->getCrackTipParents()[nextCrackTip][0]=markedElems[ed_it->second];
			target->getCrackTipStates()[nextCrackTip]=ACTIVE;
			++nextCrackTip;
		}
		//we might want to remove crack-surface intersections by moving the nodes inside (if they are not already)
		//(either along LS-normals or in-plane mesh-normals)

		target->recomputeCrackAreas();

		// re-init the high-res crack-tip
		target->getHiResCrackTip().init(source->getTargetMeshSize(),source->getHiResCrackTip().getVertsPerSegment());

		// now find crack-front edges in the source and copy over crack-front marker positions, directions, etc. ...
		target->getHiResCrackTip().copyEdges(source->getHiResCrackTip(),newCrackEdges,markedNodes,true /*reset state*/);

		return 0;
	}

	int getFrontEdgesFromMarking(id_map& markedElems, elem_map& sourceElems, edge_imap& edgs){
		// we assume all elements are numbered consistently
		// and cracks are manifolds with boundary
		// to find edges on the crack-front we look for edges appearing exactly once
		// -- we check against the reverse-edge since each interior edge
		// must appear once in each direction in the mesh
		edgs.clear();
		for(id_map::iterator el_it = markedElems.begin(); el_it != markedElems.end(); ++el_it){
			edge
				e1 (sourceElems[el_it->first][0], sourceElems[el_it->first][1]),
				e2 (sourceElems[el_it->first][1], sourceElems[el_it->first][2]),
				e3 (sourceElems[el_it->first][2], sourceElems[el_it->first][0]),
				re1(sourceElems[el_it->first][1], sourceElems[el_it->first][0]),
				re2(sourceElems[el_it->first][2], sourceElems[el_it->first][1]),
				re3(sourceElems[el_it->first][0], sourceElems[el_it->first][2]);
			edge_imap::iterator
				e1i = edgs.find( re1 ),
				e2i = edgs.find( re2 ),
				e3i = edgs.find( re3 );
			if( e1i != edgs.end() ) edgs.erase(e1i);
			else                    edgs[e1] = el_it->first;
			if( e2i != edgs.end() ) edgs.erase(e2i);
			else                    edgs[e2] = el_it->first;
			if( e3i != edgs.end() ) edgs.erase(e3i);
			else                    edgs[e3] = el_it->first;
		}
		return 0;
	}

	void findPinchedNodes(const edge_imap& newCrackEdges, id_set& pinchedNodes){
		id_map nodeCount;
		pinchedNodes.clear();
		for(edge_imap::const_iterator ed_it = newCrackEdges.begin(); ed_it != newCrackEdges.end(); ++ed_it){
			if( nodeCount.count( ed_it->first.first ) ==0 ) nodeCount[ed_it->first.first ]=1;
			else{
				if( ++nodeCount[ed_it->first.first ] > 2 ) pinchedNodes.insert( ed_it->first.first );
			}
			if( nodeCount.count( ed_it->first.second) ==0 ) nodeCount[ed_it->first.second]=1;
			else{
				if( ++nodeCount[ed_it->first.second] > 2 ) pinchedNodes.insert( ed_it->first.second);
			}
		}
	}
	
	int getTrisAndDistance(
		value_map& currentTriMinDist,  
		FractureBEM* source, FractureBEM* target, id_map newCrackIDs,
		vdb::FloatGrid::Ptr mask, vdb::Coord ijk
	){
		currentTriMinDist.clear();

		vdb::FloatGrid::Accessor acc=mask->getAccessor();
		std::deque<vdb::Coord> coordList;
		vdb::Coord c_ijk, n_ijk;
		id_set currentTris;
		vdb::Vec3d p = source->getLevelSet().getTransform()->indexToWorld(ijk);
		vdb::BBoxd bbox=source->getLevelSet().getNearTris(p, currentTris);
		//printf("\n%% mask coord (%d, %d, %d) has value %.3lg,\n%% -- bbox (%.3lf, %.3lf, %.3lf):(%.3lf, %.3lf, %.3lf) marked\n",
		//	ijk[0],ijk[1],ijk[2],acc.getValue(ijk),bbox.min()[0],bbox.min()[1],bbox.min()[2],bbox.max()[0],bbox.max()[1],bbox.max()[2]
		//);
		for(id_set::iterator el_it=currentTris.begin(); el_it!=currentTris.end();++el_it){
			if( newCrackIDs.count( source->getRegions()[*el_it] ) >0 ){ // this is a crack surface element of a known crack region
				double s,t; bool f; // unused
				currentTriMinDist[*el_it]=sqDistPtTri(
					mask->transform().indexToWorld(ijk),
					vdb::Vec3d(
						source->getNodes()[source->getElems()[*el_it][0]][0],
						source->getNodes()[source->getElems()[*el_it][0]][1],
						source->getNodes()[source->getElems()[*el_it][0]][2] ),
					vdb::Vec3d(
						source->getNodes()[source->getElems()[*el_it][1]][0],
						source->getNodes()[source->getElems()[*el_it][1]][1],
						source->getNodes()[source->getElems()[*el_it][1]][2] ),
					vdb::Vec3d(
						source->getNodes()[source->getElems()[*el_it][2]][0],
						source->getNodes()[source->getElems()[*el_it][2]][1],
						source->getNodes()[source->getElems()[*el_it][2]][2] ),
					s,t,f
				);
				//printf(" tri %d has dist %.3lg to voxel (%d, %d, %d)\n", *el_it, sqrt(currentTriMinDist[*el_it]), ijk[0],ijk[1],ijk[2]);
				//markedElems[*el_it]=0; // just put it in the map, we're not using the mapped value yet
			}
		}
		// update min-dist and deactivate voxels in BBox
		acc.setValueOff(ijk);
		coordList.push_back(ijk);
		while(!coordList.empty()){
			c_ijk=coordList.back(); coordList.pop_back();
			for(int i=0; i<26; ++i){
				vdb::Coord n_ijk = c_ijk+vdb::util::COORD_OFFSETS[i];
				if( acc.isValueOn(n_ijk) && bbox.isInside( mask->transform().indexToWorld(n_ijk) ) ){
					acc.setValueOff(n_ijk);
					coordList.push_back(n_ijk);
					// update min-dist
					for(id_set::iterator el_it=currentTris.begin(); el_it!=currentTris.end();++el_it){
						if( newCrackIDs.count( source->getRegions()[*el_it] ) >0 ){ // this is a crack surface element of a known crack region
							double s,t; bool f; // unused
							currentTriMinDist[*el_it]=std::min(currentTriMinDist[*el_it], sqDistPtTri(
								mask->transform().indexToWorld(n_ijk),
								vdb::Vec3d(
									source->getNodes()[source->getElems()[*el_it][0]][0],
									source->getNodes()[source->getElems()[*el_it][0]][1],
									source->getNodes()[source->getElems()[*el_it][0]][2] ),
								vdb::Vec3d(
									source->getNodes()[source->getElems()[*el_it][1]][0],
									source->getNodes()[source->getElems()[*el_it][1]][1],
									source->getNodes()[source->getElems()[*el_it][1]][2] ),
								vdb::Vec3d(
									source->getNodes()[source->getElems()[*el_it][2]][0],
									source->getNodes()[source->getElems()[*el_it][2]][1],
									source->getNodes()[source->getElems()[*el_it][2]][2] ),
								s,t,f
							));
						}
					}
				}
			}
		}
	return currentTriMinDist.size();
	}
}




// old stuff ...

//		// OLD VERSION of marking
//		// pre-mark all nodes that are part of crack surfaces in the source mesh
//		for(elem_map::iterator el_it = source->getElems().begin(); el_it != source->getElems().end(); ++el_it){
//			//if( source->getCracks().count( source->getRegions()[el_it->first] ) >0 ){ // this is a crack surface element
//			if( newCrackIDs.count( source->getRegions()[el_it->first] ) >0 ){ // this is a crack surface element of a known crack region
//				for(elem_map::mapped_type::iterator nd_it = el_it->second.begin(); nd_it != el_it->second.end(); ++nd_it){
//					extraNodes[*nd_it]=0; // just make sure it's in the map, we're not using the mapped value yet
//				}
//			}
//		}
//		// mark all nodes that are inside the target's implicit surface
//		for(id_map::iterator nd_it = extraNodes.begin(); nd_it != extraNodes.end(); ++nd_it){
//			std::vector<double>& nd = source->getNodes()[nd_it->first]; // reference to the current node's raw coordinate data
//			if( target->getLevelSet().isInside( nd[0], nd[1], nd[2] ) ){ // test if the node is inside
//				markedNodes[nd_it->first]=0; // put it in the marked map, we're not using the mapped value yet
//			}
//		}
//		extraNodes.clear(); // has been used as intermediate storage
//
//		// mark all elements that contain a marked node, if it also contains unmarked nodes, put these into extraNodes
//		for(elem_map::iterator el_it = source->getElems().begin(); el_it != source->getElems().end(); ++el_it){
//			bool hasMarkedNode=false;
//			for(elem_map::mapped_type::iterator nd_it = el_it->second.begin(); nd_it != el_it->second.end() && !hasMarkedNode; ++nd_it){
//				if( markedNodes.count( *nd_it ) >0 ) hasMarkedNode=true;
//			}
//			if( hasMarkedNode ){
//				markedElems[el_it->first]=0; // again not using the mapped value just yet
//				for(elem_map::mapped_type::iterator nd_it = el_it->second.begin(); nd_it != el_it->second.end(); ++nd_it){
//					if( markedNodes.count( *nd_it ) ==0 ) extraNodes[*nd_it]=0;
//				}
//			}
//		}
//		// merge markedNodes and extraNodes
//		markedNodes.insert(extraNodes.begin(), extraNodes.end());