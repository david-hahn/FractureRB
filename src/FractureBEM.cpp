/* 
 * File:   FractureBEM.cpp
 * Author: David
 * 
 * Created on 15. Jan 2014, 13:49
 */

#include "FractureBEM.h"
#include "hyena/wrap/HyENAWrapper.h"
#include "PostProcessor.h"
#include "VDBWrapper.h"
#include "FractureModel.h"
#include "SubsampledCrackTip.h"
#include "MaterialModel.h"

#include <cfloat>
#include <cmath>
#include <omp.h> // just for omp_get_wtime

namespace FractureSim{
	int FractureBEM::initBEM(MaterialModel& material, bool is2D, const std::vector<std::string>& bndCnds){
		return initBEM(material.getE(), material.getNu(), is2D, bndCnds);
	}
    int FractureBEM::initBEM(
        double youngsModulus, double poissonsRatio, bool is2D,
        const std::vector<std::string>& bndCnds
    ){
		bemSolver.reset( HyENAWrapper::constructHyENAWrapper(
            nodes, elems, regions, cracks,
            youngsModulus, poissonsRatio, is2D
        ));
        bemSolver->buildBoundaryCnd(bndCnds); // fills cracks based on bnd. cnds.
	    postPro.reset( new PostProcessor(
			nodes, elems, regions, cracks, crackTips, crackTipParents, crackTipStates, crackAreas, crackMeshSize
        ));

		if(!checkMesh()) return -1;
		updateCrackAreas(elems);

		id_set ctNodes = nodeSet(crackTips);
        int ret=1;
		if(method==FULL_BEM)    ret = bemSolver->init(ctNodes);
		if(method==INITIAL_BEM) ret = bemSolver->init(ctNodes, true /*ignore cracks*/);
		if(ret==0) bemInitialized=true;
		bemFirst=true;
        return ret;
    }

    int FractureBEM::initVDB(double voxelSize, bool noSI, double nbHWidth){
        levelSet.reset( new VDBWrapper(voxelSize,nbHWidth,noSI) );
		int ret = levelSet->init(nodes,elems,regions,cracks);
		if(ret==0){
            // initialize crackTipStates based on inside/outside object
            for(elem_map::const_iterator it = crackTips.begin();
                it!= crackTips.end(); ++it
            ){
                unsigned int nd_a, nd_b;
                nd_a = it->second[0]; nd_b = it->second[1];
                vdb::Vec3d a(nodes[nd_a][0],nodes[nd_a][1],nodes[nd_a][2]);
                vdb::Vec3d b(nodes[nd_b][0],nodes[nd_b][1],nodes[nd_b][2]);
                if((levelSet->isInside(a) ||
                    levelSet->isInside(b))&&
                    levelSet->isInside((a+b)/2.0)
                )
                    crackTipStates[it->first]=ACTIVE;
                else
                    crackTipStates[it->first]=INACTIVE;
            }
			vdbInitialized=true;
        }
        return ret;
    }

	int FractureBEM::initFractureModel(double strength, double toughness, double density, double compressive){
		if(!bemInitialized) return -1;
		MaterialModel* mat = createMaterialModel( "default",
			bemSolver->getYoungsModulus(), bemSolver->getPoissonsRatio(), density,
			strength, toughness, compressive
		);
		int ret=initFractureModel(*mat);
		delete mat;
		postPro->setCompressiveFactor(compressive);
		return ret;
	}
	int FractureBEM::initFractureModel(MaterialModel& material, bool substepOutput, double samplingDensity){
		if(!vdbInitialized) return -1;
		materialModel.reset( material.clone() );
		fractureModel.reset( new DynamicVelocityFractureModel(*materialModel,  crackMeshSize));
        fineCrackTip.reset( new SubsampledCrackTip(
            nodes,elems,regions,
            crackTips,crackTipParents,crackTipStates,
            levelSet.get(),fractureModel.get()
        ));
		fineCrackTip->init(
            crackMeshSize,
            std::max<int>(1, samplingDensity*ceil( crackMeshSize / levelSet->getVoxelSize() ) ),
            substepOutput
        ); //number of subsamples per crack tip
        fractureInitialized=true;
        
        //also initialize VDB near triangle grid
        levelSet->addToNearTriGrid(nodes,elems,crackMeshSize);

		return 0;
	}

	int FractureBEM::computeBEM(){
		if(!bemInitialized) return -1;
		if(method==FULL_BEM){
			bemSolver->compute();

			postPro->computeNodeSIFs(
				bemSolver->getDisplacements(), bemSolver->getYoungsModulus(), bemSolver->getPoissonsRatio(),
				nodeSIFs, crackTipFaceNormals, crackTipTangents
			);
		}
		if(method==INITIAL_BEM){
			if(bemFirst) bemSolver->compute();
			bemFirst=false;

			state_map nodeStates;
			int activeNodes = fineCrackTip->getCoarseCrackFrontStates( nodeStates );
			//printf("\n%% ... %d of %d crack front nodes are active", activeNodes, nodeStates.size());
			postPro->estimateNodeSIFs(
				bemSolver->getDisplacements(), bemSolver->getTractions(), bemSolver->getYoungsModulus(),
				bemSolver->getPoissonsRatio(), nodeSIFs, crackTipFaceNormals, crackTipTangents, nodeStates
			);
		}
		return 0;
	}

	int FractureBEM::computeAddedCracks(){
		if(method==FULL_BEM){
			id_set ctNodes = nodeSet(crackTips);
			int ret = bemSolver->addCrack(ctNodes);
			if(ret!=-1) ret=computeBEM();
			return ret;
		}
		return 0;
	}

	int FractureBEM::seedCracksAndPropagate(int maxSeed){
		if(!vdbInitialized || !bemInitialized || !fractureInitialized) return -1;
		double t0=omp_get_wtime(), t1;
		int newCracks=0; // count how many new cracks have been seeded
		// - compute max. pricipal stresses and their directions for all triangles
		// - check which of these should start a new edge-crack
		// - seed those cracks
		// - update the BEM solution
		// - run propagateCracks()
        printf("\n%% crack initiation  ...");
		id_set ctNodes = nodeSet(crackTips);
		vector_type empty;
		node_map psn; // maps node-/element-IDs to 4 double values each, first value is the max. principal stress, last 3 values are it's direction (unit vector)
		int ret = postPro->computeSurfaceStresses(
			psn,bemSolver->getDisplacements(), (method==INITIAL_BEM)? empty : bemSolver->computeCrackBaseDisplacements(), //don't use crack-base disp in method INIT_BEM !!!
			bemSolver->getTractions(), materialModel->getE(),materialModel->getNu(),
			(method==INITIAL_BEM)?firstCrackElement:-1 //if method==INIT_BEM skip crack-elements in all of this
		);
		if(ret!=0) return ret;
        //printf("\n%% ... ... surface stresses done");

		typedef std::multimap<double, long> mm_t; mm_t mm; // map stress excess to SIGNED element ID
		vect3d_map stressNormals;
		std::set<unsigned int> negSidesMax, negSidesMin; // store which cracks should start in the reversed direction (branches from the "negative" fracture surface)
		for(node_map::iterator it=psn.begin(); it!=psn.end(); ++it){
			Eigen::Vector3d c ( // centroid of element
				(nodes[elems[it->first][0]][0]+nodes[elems[it->first][1]][0]+nodes[elems[it->first][2]][0])/3.0,
				(nodes[elems[it->first][0]][1]+nodes[elems[it->first][1]][1]+nodes[elems[it->first][2]][1])/3.0,
				(nodes[elems[it->first][0]][2]+nodes[elems[it->first][1]][2]+nodes[elems[it->first][2]][2])/3.0
			);
			if( it->second[0] > materialModel->tensileStrength(c) && //strength and toughness are quite unrelated material params - depends on mat. type, production process, etc...
				fracturedElems.count(it->first)==0 && // don't fracture an element more than once
				ctNodes.count(elems[it->first][0])==0 && // don't fracture elements node-adjacent to the crack tip
				ctNodes.count(elems[it->first][1])==0 &&
				ctNodes.count(elems[it->first][2])==0 ){ // candidate element to seed a crack
				mm.insert(mm_t::value_type(it->second[0]-materialModel->tensileStrength(c),it->first)); // insert into multimap -> sort by stress excess asc.
                Eigen::Vector3d	n(it->second[1],it->second[2],it->second[3]);
                stressNormals[it->first]=n;
				if( it->second.size()>8){ // check for negative side flag on max. principal stress (2 or 4)
					if( std::abs(it->second[8]-2.0)<0.1 || std::abs(it->second[8]-4.0)<0.1 ){
						negSidesMax.insert(it->first);
					}
				}
			}
			if(it->second.size()>7) // compressive fracture seeding works if computeSurfaceStresses returns min. principal stress along with max. ps.
                if( fracturedElems.count(it->first)==0 && // don't fracture more than once in the same location
                    ctNodes.count(elems[it->first][0])==0 && // don't fracture elements node-adjacent to the crack tip
                    ctNodes.count(elems[it->first][1])==0 &&
                    ctNodes.count(elems[it->first][2])==0 &&
                    -it->second[4] > materialModel->tensileStrength(c)*materialModel->compressiveFactor(c)
                ){ // candidate to seed a crack
                    mm.insert(mm_t::value_type((-it->second[4]/materialModel->compressiveFactor(c))-materialModel->tensileStrength(c),-(it->first))); // insert negative elem-ID into multimap -> sort by stress excess asc.
                    Eigen::Vector3d	n(it->second[5],it->second[6],it->second[7]);
                    stressNormals[it->first]=n;
					if( it->second.size()>8){ // check for negative side flag on min. principal stress (3 or 4)
						if( std::abs(it->second[8]-3.0)<0.1 || std::abs(it->second[8]-4.0)<0.1 ){
							negSidesMin.insert(it->first);
						}
					}
                }

		}
        //printf("\n%% ... ... candidates listed");

		// reverse-iterate candidates (max stress first),     create at most this --v number of cracks (maxSeed < 0 is unlimited)
        for(mm_t::reverse_iterator it=mm.rbegin(); it!=mm.rend() && ((newCracks < maxSeed) || (maxSeed<0)); ++it){
			unsigned int elem = std::abs(it->second);
			bool negSide=false;
			if(checkDistanceEl(elem,ctNodes)){ // check whether the element is far enough from all crack-tip nodes
                //printf("\n%% ... ... seeding from element %d", it->second);
				if(it->second > 0 && negSidesMax.count(elem)>0) negSide=true;
				if(it->second < 0 && negSidesMin.count(elem)>0) negSide=true;
				ret = startCrack(elem,stressNormals[elem],negSide);
				if(ret!=0) return -1;
				fracturedElems.insert(elem); // don't fracture this element again
				++newCracks;
			}
			ctNodes = nodeSet(crackTips);

			//printf("\n%% fractured element %d", it->second);
			// problems: this would potentially start (almost) coincident cracks in each step
			//           * taken care of by checking fracturedElems --> each element may only start a single crack
			//           * also do not branch off triangles that contain a crack-tip node (these touch the stress singularity)
			//           this would potentially start a lot of cracks either at once or in similar places at each step
			//           * taken care of by checking the distances from the triangle's nodes to all crack-tip nodes
			//           * and only allowing a new crack if all these distances are longer than the triangle's longest edge
			// could do: 
			//           check the level-set if we're too close to an existing crack (limited by NBHW to background value)
			//           tri-tri distance computation (expensive?)
			//           tri-point distance from candidate tris to crack-tip nodes?
			//           ^-- could do this based only on the min vertex-point distance since anything shorter than
			//               the longest edge of the triangle is too close anyway
		}
        //printf("\n%% ... ... seeding done");

        if(newCracks>0){ // added at least 1 new crack
			// update the last written VDB file
			if(!lastVDBfile.empty()) levelSet->writeGrids(lastVDBfile);
            // update BEM solution
            //printf("\n%% ... ... updating BEM");
            //ctNodes = nodeSet(crackTips);
			if(method == FULL_BEM){
				ret = bemSolver->addCrack(ctNodes); // update matrices
				computeBEM(); // update solution
			}
			// if method==INITIAL_BEM
			// we don't update here, since SIF evaluations are the expensive part (and we don't need the matrix assembly anymore)
			//ToDo: evaluate SIFs ONLY at new cracks here, without changing already done ones
        }
		t1=omp_get_wtime();
		printf("\t%.4lfs",t1-t0);
		t0=t1;

		ret = propagateCracks();
		if(ret!=0) return ret;
		return newCracks;
	}

	int FractureBEM::propagateCracks(){
		if(!vdbInitialized || !bemInitialized || !fractureInitialized) return -1;
		double t0=omp_get_wtime(), t1;
		node_map newNodes;
		elem_map newElems, newCrackTips, newCrackTipParents;
        id_map newRegions;
        state_map newCrackTipStates;
		//printf("\n%% crack propagation ...");

		//// for testing *********************
		//if(0 && method==FULL_BEM){
		//	postPro->estimateNodeSIFs(
		//		bemSolver->getDisplacements(), bemSolver->getTractions(), bemSolver->getYoungsModulus(),
		//		bemSolver->getPoissonsRatio(), nodeSIFs, crackTipFaceNormals, crackTipTangents
		//	);
		//	fineCrackTip->addToSIFsOutput(nodeSIFs,"reg_");
		//	postPro->computeNodeSIFs(
		//		bemSolver->getDisplacements(), bemSolver->getYoungsModulus(), bemSolver->getPoissonsRatio(),
		//		nodeSIFs, crackTipFaceNormals, crackTipTangents
		//	);
		//}
		//else if(0 && method==FULL_BEM){ // reverse test
		//	fineCrackTip->addToSIFsOutput(nodeSIFs,"reg_");
		//	postPro->estimateNodeSIFs(
		//		bemSolver->getDisplacements(), bemSolver->getTractions(), bemSolver->getYoungsModulus(),
		//		bemSolver->getPoissonsRatio(), nodeSIFs, crackTipFaceNormals, crackTipTangents
		//	);
		//}
		//// ********************* end of testing

		printf("\n%% sub-stepping      ...");
		int ret = fineCrackTip->propagate(
			newNodes,newElems,newRegions,
			newCrackTips,newCrackTipParents,newCrackTipStates,
			crackTipFaceNormals,crackTipTangents,nodeSIFs,
            std::max<int>(1, ceil( crackMeshSize / levelSet->getVoxelSize() ) )
		);
		t1=omp_get_wtime();
		printf("\t%.4lfs",t1-t0);
		t0=t1;

		// even if no new geometry is added, these maps can update
		crackTips=newCrackTips; // these are copy-assign operations, "old" references to maps stay valid (with updated values)
		crackTipParents=newCrackTipParents;
		crackTipStates=newCrackTipStates;
		if( ret!=0 ) printf("\n%% !!! propagate(...) returned %d !!!\n",ret);
		else if( !newElems.empty() ){
            nodes.insert(newNodes.begin(), newNodes.end());
            elems.insert(newElems.begin(), newElems.end());
            regions.insert(newRegions.begin(), newRegions.end());

            levelSet->addToNearTriGrid(nodes,newElems);
            
			updateCrackAreas(newElems);

			t1=omp_get_wtime();
			printf("\n%% BEM mesh update   ...\t%.4lfs",t1-t0);
			t0=t1;

			id_set ctNodes = nodeSet(crackTips);
			if(method == FULL_BEM){
				printf("\n%% matrix assembly   ...");
				ret = bemSolver->addCrack(ctNodes);
				t1=omp_get_wtime();
				printf("\t%.4lfs",t1-t0);
				t0=t1;
			}else if(method == INITIAL_BEM){
				vect3d_map newCOD;
				postPro->estimateCOD(newElems,nodeSIFs,ctNodes,materialModel->getE(),materialModel->getNu(), newCOD);

				elem_map support_elems;
				id_set support_crack_tip_nodes;
				for(elem_map::iterator it=elems.find(firstCrackElement); it!=elems.end(); ++it){
					if( newCOD.count( it->second[0] )>0 ||
						newCOD.count( it->second[1] )>0 ||
						newCOD.count( it->second[2] )>0 )
						support_elems[it->first]=it->second; // this is a support element if it contains a node for which we have new COD
				}
				support_crack_tip_nodes = nodeSet(support_elems);
				for(vect3d_map::iterator it=newCOD.begin(); it!=newCOD.end(); ++it){
					support_crack_tip_nodes.erase(it->first); // virtual crack tip nodes in the support elems are all nodes for which we don't have new COD
				}
				bemSolver->updateRHSwithCOD(newCOD, support_elems, support_crack_tip_nodes);

				//for testing: write new CODs into displacement result
				vector_type& u = bemSolver->getDisplacements();
				unsigned int u_tmp=u.rows();
				u.conservativeResize(3*nodes.size());
				u.bottomRows(u.rows()-u_tmp).setZero();
				// add disp update to displacements
				for(vect3d_map::iterator it=newCOD.begin(); it!=newCOD.end(); ++it){
					u_tmp = 3*(it->first-NODE_BASE_INDEX);
					u[u_tmp  ]+=it->second[0];
					u[u_tmp+1]+=it->second[1];
					u[u_tmp+2]+=it->second[2];
				}
				t1=omp_get_wtime();
				printf("\n%% RHS update        ...\t%.4lfs",t1-t0);
				t0=t1;
			}
		}
		printf("\n%%                    total ...");
		return ret;
	}

	unsigned int FractureBEM::getActiveCount(){
		unsigned int count=0;
		for(state_map::iterator it = crackTipStates.begin();
            it != crackTipStates.end(); ++it
        ){
			if( it->second == ACTIVE   ||
				it->second == ACTIVE_A ||
				it->second == ACTIVE_B ) ++count;
		}
		return count;
	}

	void FractureBEM::setCracksInactive(const id_set& offRegions){
        for(elem_map::iterator it=crackTipParents.begin();
            it!=crackTipParents.end(); ++it){
            if(offRegions.count( regions[it->second[0]] )>0)
                crackTipStates[it->first]=INACTIVE;
        }
    }

	int FractureBEM::writeMesh(std::string filename,
            double visualQuality, bool visDisplace, bool visCOD, bool visClose, bool visOBJ
    ){
		int ret=0;
		if(method==FULL_BEM){
			vector_type& u = bemSolver->getDisplacements();
			vector_type& q = bemSolver->getTractions();
			vector_type ps, pn, psm, pnm; //principal stress and normals
			vector_type u_std, u_half; //post-processed displacements (split into continous part and half-opening)
			vector_type regionVect; //region IDs per elem

			postPro->vtkClearData(); // just in case we left old stuff in there
			postPro->vtkAddData("displacement",3,false,u);
			postPro->vtkAddData("traction",    3,true ,q);

			// write region-IDs to mesh
			regionVect.resize( elems.size() );
			for(int i=0; i<elems.size(); ++i){
				regionVect[i]=regions[i+ELEM_BASE_INDEX];
			}
			postPro->vtkAddData("region",1,true,regionVect);

			// generate principal stress data for the VTK output
			node_map maxPrincipalStress;
			ret = postPro->computeSurfaceStresses(
				maxPrincipalStress,
				bemSolver->getDisplacements(),
				bemSolver->computeCrackBaseDisplacements(),
				bemSolver->getTractions(),
				materialModel->getE(),
				materialModel->getNu(), true /*addToVTK -> also adds cartesian stresses to the VTK data*/
			);
			if( ret<0 ) printf(" ** computeSurfaceStresses returned %d ", ret);
			else{
				ps.resize(elems.size()); pn.resize(3*elems.size());
				psm.resize(elems.size()); pnm.resize(3*elems.size());
				int i=0;
				for(elem_map::iterator it=elems.begin(); it!=elems.end(); ++it,++i){
					ps[i]     =maxPrincipalStress[it->first][0]; // max. principal stress
					pn[3*i  ] =maxPrincipalStress[it->first][1]; // normal ...
					pn[3*i+1] =maxPrincipalStress[it->first][2];
					pn[3*i+2] =maxPrincipalStress[it->first][3];
					psm[i]    =maxPrincipalStress[it->first][4]; // min. princiap stress
					pnm[3*i  ]=maxPrincipalStress[it->first][5]; // normal ...
					pnm[3*i+1]=maxPrincipalStress[it->first][6];
					pnm[3*i+2]=maxPrincipalStress[it->first][7];
				}
				postPro->vtkAddData("max_principal_stress"       ,1,true,ps);
				postPro->vtkAddData("max_principal_stress_normal",3,true,pn);
				postPro->vtkAddData("min_principal_stress"       ,1,true,psm);
				postPro->vtkAddData("min_principal_stress_normal",3,true,pnm);
			}

			double t=omp_get_wtime();
			/**/ // extrapolate cod of inactive crack tips from interior nodes
			id_map count;
			for(elem_map::iterator it=crackTips.begin(); it!=crackTips.end(); ++it)
				if( crackTipStates[it->first] == INACTIVE ){ // extrapolate cod of inactive crack tips from interior nodes
					unsigned int nd_a=it->second[0], nd_b=it->second[1];
					unsigned int nd_c = postPro->findInteriorNode(nd_a,nd_b,elems[crackTipParents[it->first][0]]);
					u[3*(nd_a-NODE_BASE_INDEX)  ]+=u[3*(nd_c-NODE_BASE_INDEX)  ];
					u[3*(nd_a-NODE_BASE_INDEX)+1]+=u[3*(nd_c-NODE_BASE_INDEX)+1];
					u[3*(nd_a-NODE_BASE_INDEX)+2]+=u[3*(nd_c-NODE_BASE_INDEX)+2];
					u[3*(nd_b-NODE_BASE_INDEX)  ]+=u[3*(nd_c-NODE_BASE_INDEX)  ];
					u[3*(nd_b-NODE_BASE_INDEX)+1]+=u[3*(nd_c-NODE_BASE_INDEX)+1];
					u[3*(nd_b-NODE_BASE_INDEX)+2]+=u[3*(nd_c-NODE_BASE_INDEX)+2];
					++count[nd_a]; ++count[nd_b];
				}
			for(id_map::iterator it=count.begin(); it!=count.end(); ++it){
				u[3*(it->first - NODE_BASE_INDEX)  ]/=it->second;
				u[3*(it->first - NODE_BASE_INDEX)+1]/=it->second;
				u[3*(it->first - NODE_BASE_INDEX)+2]/=it->second;
			}
			//printf("\n%% extrapolating COD ...\t%.4lfs",omp_get_wtime()-t); t=omp_get_wtime(); /**/

			vector_type& u_c = bemSolver->computeCrackBaseDisplacements(); // do this after extrapolation to get additional crack-crack interactions
			postPro->vtkAddData("crack_base_displacement",3,false,u_c);
			/**/ // add full coarse displacement data to VTK output
			u_std=u; u_half=u; //copy-assign
			id_set done; unsigned int nd3;
			for(elem_map::iterator it=elems.begin(); it!=elems.end(); ++it){
				for(int k=0; k < it->second.size(); ++k){
					if( done.count( it->second[k] ) ==0 ){ // if a node is not done yet
						nd3 = 3*(it->second[k]-NODE_BASE_INDEX);
						if( cracks.count( regions[it->first] ) ==0 ){ // ordinary surface
							u_half[nd3   ] =0.0;
							u_half[nd3 +1] =0.0;
							u_half[nd3 +2] =0.0;
						}else{ // crack
							u_std[nd3   ]  =0.0;
							u_std[nd3 +1]  =0.0;
							u_std[nd3 +2]  =0.0;
							u_half[nd3   ]*=0.5;
							u_half[nd3 +1]*=0.5;
							u_half[nd3 +2]*=0.5;
						}
						done.insert( it->second[k] );
					}
				}
			}
			postPro->vtkAddData("crack_half_opening_displacement",3,false,u_half);
			postPro->vtkAddData("surface_only_displacement",3,false,u_std);
			//printf("\n%% processing full u ...\t%.4lfs\n%%            total ...         ",omp_get_wtime()-t); /**/
			ret=postPro->vtkWrite(filename+".vtk",VTK_TRIANGLE);
			postPro->vtkClearData(); // clean up as most references will not be valid beyond the scope of this function

			if(visualQuality<=1.0){
				double vdbQual=visualQuality, vcgError=-1.0;
				if(visualQuality<0.0){ vdbQual=0.001; vcgError=-visualQuality; }
				levelSet->writeVisualMesh(
					nodes, elems, regions, cracks,
					u, u_c, filename, vdbQual, vcgError,
					visDisplace, visCOD, visClose, visOBJ
				);
			}
		}
		if(method==INITIAL_BEM){
			vector_type u = bemSolver->getDisplacements(); // need copies here!
			vector_type q = bemSolver->getTractions();
			vector_type regionVect; //region IDs per elem

			unsigned int u_vals=u.rows(), q_vals=q.rows();
			u.conservativeResize(3*nodes.size());
			q.conservativeResize(3*elems.size());
			u.bottomRows(u.rows()-u_vals).setZero();
			q.bottomRows(q.rows()-q_vals).setZero();
			//printf("\n\n%% **** u size was %d is now %d, diff %d ****",  u_vals,u.rows(),u.rows()-u_vals);
			//printf(  "\n%% **** q size was %d is now %d, diff %d ****\n",q_vals,q.rows(),q.rows()-q_vals);

			postPro->vtkClearData(); // just in case we left old stuff in there
			postPro->vtkAddData("displacement",3,false,u);
			postPro->vtkAddData("traction",    3,true ,q);

			vector_type& u_c = bemSolver->computeCrackBaseDisplacements();
			postPro->vtkAddData("crack_base_displacement",3,false,u_c);

			// write region-IDs to mesh
			regionVect.resize( elems.size() );
			for(int i=0; i<elems.size(); ++i){
				regionVect[i]=regions[i+ELEM_BASE_INDEX];
			}
			postPro->vtkAddData("region",1,true,regionVect);
			
			ret=postPro->vtkWrite(filename+".vtk",VTK_TRIANGLE);
			postPro->vtkClearData(); // clean up as most references will not be valid beyond the scope of this function
		}
		return ret;
	}

	int FractureBEM::startCrack(unsigned int elem, Eigen::Vector3d normal, bool negSide){
		if(!fractureInitialized) return -1;
		// to add a crack do the following:
		// - find an anchor point (i.e. the centre of the given element)
		// - find an unused region number
		// - create a geometry around the anchor point in the plane given by normal
		// -- such that it has at least 1 free (not-crack-tip) node
		// -- has edges roughly of the length of the mesh-resolution
		// -- add this to the mesh with a new region-ID
		// - update the BEM solver to account for the new crack regions
		// - update fine crack tip
		// -- generate vertices
		// -- initialize their state
		// -- create a new level-set grid for the new crack
		Eigen::Vector3d // a,b,c are node coordinates of the given element
			a (nodes[elems[elem][0]][0], nodes[elems[elem][0]][1], nodes[elems[elem][0]][2]),
			b (nodes[elems[elem][1]][0], nodes[elems[elem][1]][1], nodes[elems[elem][1]][2]),
			c (nodes[elems[elem][2]][0], nodes[elems[elem][2]][1], nodes[elems[elem][2]][2]),
			p, f; // p stores the anchor point of the new crack, f its "forward" direction (in-plane-normal to crack-tip)

		p=(a+b+c)/3.0;
        
        normal.normalize(); // just to be sure it's a unit vector
		f=-(a-c).cross(b-c).normalized(); // start with the (inward==negative) face-normal of the triangle
		if(negSide) f*=-1.0;
		f-=f.dot(normal)*normal; // remove out-of-plane component
		if( f.dot(f) < FLT_EPSILON ) return -1; // if f and normal are almost-parallel the whole idea of starting such a crack is nonsense
		f.normalize();
		return startCrack(p,normal,f);
	}
	
	int FractureBEM::startCrack(Eigen::Vector3d p, Eigen::Vector3d normal, Eigen::Vector3d f){
		if(!fractureInitialized) return -1;
		unsigned int s=0;
		for(id_map::iterator it=regions.begin(); it!=regions.end(); ++it){
			s=std::max(s,it->second); // find highest region-id in regions map
		}
		for(id_set::iterator it=cracks.begin(); it!=cracks.end(); ++it){
			s=std::max(s,*it); // find highest region-id in cracks
		}
		++s; // this will be the next region-id
        //printf(" new crack region %d", s);
		
		elem_map newElems;
		int ret = fineCrackTip->startCrack( newElems,
			p+f*levelSet->getVoxelSize()*0.5,normal,f,s
		); // this writes the geometry update into the BEM mesh
        //printf(" ... done, returned %d", ret);
		if(ret==0){
            // apply crack BC to new surface
            bemSolver->addCrackBoundaryCnd(s); // this also adds the new region-id to the cracks set!
			//the caller is responsible for updating the BEM solution after starting new crack(s)
			updateCrackAreas(newElems);
		}
		return ret;
	}

	int FractureBEM::remesh(int nTris, double adaptive, double offsetVoxels){
		if(vdbInitialized && nTris>3)
			return levelSet->mesh(nodes,elems,nTris,adaptive,offsetVoxels);
		return -1;
	}

	int FractureBEM::writeVDB(std::string filename, bool updateSeeds){ //NEW FOR FractureRB: option to suppress updates when seeding cracks
        fineCrackTip->setFilename(filename);
		if(vdbInitialized){
			if( updateSeeds ) lastVDBfile = filename;
			else lastVDBfile.clear();
			return levelSet->writeGrids(filename);
		}
		return -1;
	}
    int FractureBEM::writeCrackTip(std::string filename){
        if(fractureInitialized)
            return fineCrackTip->vtkWrite(filename+".vtk");
        return -1;
    }
	int FractureBEM::innerEvalTest(std::string filename){
		if(bemInitialized)
			return bemSolver->innerEvalTest(filename);
		return -1;
	}
	// when using smart-pointers (boost) with incomplete types, constructor and destructor should NOT be inlined in the header!
	FractureBEM::FractureBEM(double resolution, FRACTURE_METHOD m){
		crackMeshSize      =resolution;
		bemInitialized     =false;
		vdbInitialized     =false;
        fractureInitialized=false;
		if( m==FULL_BEM || m==INITIAL_BEM ){ // currently we only support these two methods
			method = m;
		}else{
			printf("\n!!! invalid fracture method %d, choose either full BEM (%d) or initial BEM solution (%d). Switching to FULL BEM method.",
				m,FULL_BEM,INITIAL_BEM
			);
			method=FULL_BEM;
		}
		lastVDBfile.clear();
	}
	FractureBEM::~FractureBEM(){} // smart-pointers take care of object deletions

	bool FractureBEM::checkDistanceEl(unsigned int elem, const id_set& nodeSet){
		// return true iff all vertices of the given element are far enough away from all given nodes
		Eigen::Vector3d // a,b,c are node coordinates in material space; ua,ub,uc are nodal displacements, so world coordinates are a+ua etc.
			a (nodes[elems[elem][0]][0], nodes[elems[elem][0]][1], nodes[elems[elem][0]][2]),
			b (nodes[elems[elem][1]][0], nodes[elems[elem][1]][1], nodes[elems[elem][1]][2]),
			c (nodes[elems[elem][2]][0], nodes[elems[elem][2]][1], nodes[elems[elem][2]][2]),
			tmp;
		double lengthSqr = std::max((b-a).squaredNorm(),std::max((c-a).squaredNorm(),(c-b).squaredNorm())); // use longest edge
		//double lengthSqr = std::min((b-a).squaredNorm(),std::min((c-a).squaredNorm(),(c-b).squaredNorm())); // use shortest edge
		//double lengthSqr = crackMeshSize*crackMeshSize; // use target mesh-resolution as cutoff
        for(id_set::iterator it=nodeSet.begin(); it!=nodeSet.end(); ++it){
			tmp[0]=nodes[*it][0]; tmp[1]=nodes[*it][1]; tmp[2]=nodes[*it][2];
			if( (tmp-a).squaredNorm() < lengthSqr) return false;
			if( (tmp-b).squaredNorm() < lengthSqr) return false;
			if( (tmp-c).squaredNorm() < lengthSqr) return false;
		}
		return true;
	}

	bool FractureBEM::checkDistanceNd(unsigned int node, const id_set& nodeSet){
		Eigen::Vector3d p(nodes[node][0], nodes[node][1], nodes[node][2]), tmp;
        double lengthSqr = crackMeshSize*crackMeshSize;
		for(id_set::iterator it=nodeSet.begin(); it!=nodeSet.end(); ++it){
			tmp[0]=nodes[*it][0]; tmp[1]=nodes[*it][1]; tmp[2]=nodes[*it][2];
			if( (tmp-p).squaredNorm() < lengthSqr) return false;
		}
		return true;
	}

	void FractureBEM::setCPVersion(unsigned int v){
		if(v<=2) FractureModel::cpVersion=v;
	}

	
	bool FractureBEM::checkMesh(){
		unsigned int r;
		bool surface=true;
		firstCrackElement=0;
		//printf("\n%% checking the mesh for consistency ...");
		for(elem_map::iterator it=elems.begin(); it!=elems.end(); ++it){
			if( regions.count(it->first)==1 ) r=regions[it->first];
			else{ printf("\n%% element %d has no region", it->first); return false;}

			if( cracks.count(r)==1 ){ // this is a crack element
				if( surface ){ // still within surface elements range
					firstCrackElement=it->first; // this is the first crack element we found
					surface=false;
				}
			}else{ // not a crack
				if(!surface ){ printf("\n%% surface element %d has ID higher than first crack element (%d) - reorder the mesh!",it->first,firstCrackElement); return false;}
			}

			if( it->second.size()!=3 ){ printf("\n%% element %d has %d nodes (should have 3)",it->first, it->second.size()); return false;}
			for(unsigned int k=0; k<3; ++k){
				if( nodes.count(it->second[k])!=1 ){ printf("\n%% element %d references non-existent node %d",it->first, it->second[k]); return false;}
			}
		}
		if( firstCrackElement==0 ) firstCrackElement = elems.size()+ELEM_BASE_INDEX; // we have no cracks in the mesh yet, the next element will become the first crack element
		return true; // the mesh is ok
	}

	// NEW FOR FractureRB
	void FractureBEM::updateCrackAreas(const elem_map& newElems){
		unsigned int region;
		double area;
		for(elem_map::const_iterator it=newElems.begin(); it!=newElems.end(); ++it){
			// check if region is a crack
			// if so, compute the area of the triangle
			// and add it to crackAreas[region]
			region=regions[it->first];
			if( cracks.count(region)!=0){ // this is a crack element
				Eigen::Vector3d // a,b,c are node coordinates in material space
					a (nodes[it->second[0]][0], nodes[it->second[0]][1], nodes[it->second[0]][2]),
					b (nodes[it->second[1]][0], nodes[it->second[1]][1], nodes[it->second[1]][2]),
					c (nodes[it->second[2]][0], nodes[it->second[2]][1], nodes[it->second[2]][2]);
				area = (a-c).cross(b-c).norm()*0.5;
				if( crackAreas.count(region)==0 ) crackAreas[region]=0.0; // initialize on first access
				crackAreas[region]+=area;
			}
		}
	}

	// NEW FOR FractureRB
	unsigned int FractureBEM::findClosestSurfaceTri(const Eigen::Vector3d& p){
		return levelSet->findClosestSurfaceTri(p, nodes, elems, regions, cracks);
	}
	// NEW FOR FractureRB
	vector_type& FractureBEM::getDisplacements(){return bemSolver->getDisplacements();}
	vector_type& FractureBEM::getTractions(){return bemSolver->getTractions();}
	// NEW FOR FractureRB
	double FractureBEM::getTimeStep(){
		if( fractureInitialized ) return fractureModel->getTimeStep();
		return 0.0;
	}
	// NEW FOR FractureRB
	int FractureBEM::updateBoundaryData(vect3d_map& newBcData){
		if( bemInitialized ){
			bemFirst=true;
			id_set ctNodes = nodeSet(crackTips);
			return bemSolver->updateRHS(newBcData,ctNodes);
		}
		return -1;
	}
	// NEW FOR FractureRB
	VDBWrapper& FractureBEM::getLevelSet(){
		return *levelSet;
		//return levelSet.get();
	}
	// NEW FOR FractureRB
	SubsampledCrackTip& FractureBEM::getHiResCrackTip(){
		return *fineCrackTip;
	}

	// just for testing ...
	void FractureBEM::dumpDebugData(){
		printf("\nel=[\n");
		printMap(elems);
		printf("];\nco=[\n");
		printMap(nodes);
		printf("];\n u=[\n");
		vector_type& u = bemSolver->getDisplacements();
		for(int i=0;i<u.size();++i) printf("%12.4le%c",u[i],(i%3==2)?'\n':' ');
		printf("];\n q=[\n");
		vector_type& q = bemSolver->getTractions();
		for(int i=0;i<q.size();++i) printf("%12.4le%c",q[i],(i%3==2)?'\n':' ');
		printf("];\n");
		
		// test stress eval ...
		vect3d_map sigma,tau;
		std::vector<Eigen::Vector3d> eval_points;
		Eigen::Vector3d h,cmin(nodes.begin()->second[0],nodes.begin()->second[1],nodes.begin()->second[2]),
			cmax(nodes.begin()->second[0],nodes.begin()->second[1],nodes.begin()->second[2]);
		for(node_map::iterator ni=nodes.begin(); ni!=nodes.end(); ++ni){
			cmin[0]=std::min(cmin[0], ni->second[0]); cmin[1]=std::min(cmin[1], ni->second[1]); cmin[2]=std::min(cmin[2], ni->second[2]);
			cmax[0]=std::max(cmax[0], ni->second[0]); cmax[1]=std::max(cmax[1], ni->second[1]); cmax[2]=std::max(cmax[2], ni->second[2]);
		}
		unsigned int N_pts=16, tmp=0;
		h=(cmax-cmin)/(N_pts);
		//The Visualization Toolkit supports five different dataset formats: structured points, structured grid, rectilinear 
		//grid, unstructured grid, and polygonal data. Data with implicit topology (structured data such as vtkImageData and 
		//vtkStructuredGrid) are ordered with x increasing fastest, then y, then z. http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.html
		eval_points.assign(N_pts*N_pts*N_pts, Eigen::Vector3d::Zero());
		for(        int i=0; i<N_pts;++i) {
			for(    int j=0; j<N_pts;++j) {
				for(int k=0; k<N_pts;++k) {
					eval_points[tmp][0]=cmin[0]+h[0]*(0.5+k);
					eval_points[tmp][1]=cmin[1]+h[1]*(0.5+j);
					eval_points[tmp][2]=cmin[2]+h[2]*(0.5+i);
					++tmp;
				}
			}
		}

		double t1,t0 = omp_get_wtime();
		vect3d_map unused;
		postPro->computeInteriorStresses(sigma,tau,unused,eval_points,u,q,materialModel->getE(),materialModel->getNu());
		t1 = omp_get_wtime();
		printf("\n%% stress eval for %d points done in %.4lfs",eval_points.size(), t1-t0);
		//for(vect3d_map::iterator it=eval_points.begin(); it!=eval_points.end(); ++it){
		//	printf("\n%% stress at (%.2lf,%.2lf,%.2lf) is\n%% "
		//		   "(%12.4le, %12.4le, %12.4le\n%% %14c %12.4le, %12.4le\n%% %14c %13c %12.4le )"
		//		   "\n%% evaluation done in %.4lfs",
		//		it->second[0],it->second[1],it->second[2],
		//		sigma[it->first][0],tau[it->first][0],tau[it->first][1],
		//		' ', sigma[it->first][1],tau[it->first][2],
		//		' ', ' ', sigma[it->first][2], t1-t0
		//	);
		//}
		std::ofstream out("test_int_stress.vtk");
        if( out.is_open()){
			out.precision(12);
			//out.setf(std::ios::scientific);
			//write header
			out << "# vtk DataFile Version 2.0" << std::endl;
			out << "VTK inner evaluation" << std::endl;
			out << "ASCII" << std::endl;
			out << "DATASET STRUCTURED_POINTS" << std::endl;
			out << "DIMENSIONS " << N_pts << " " << N_pts << " " << N_pts << std::endl;
			out << "ORIGIN " << cmin[0]+h[0]*0.5 << " " << cmin[1]+h[1]*0.5 << " " << cmin[2]+h[2]*0.5 << std::endl;
			out << "SPACING " << h[0] << " " << h[1] << " " << h[2] << std::endl;

			out << "POINT_DATA " << N_pts*N_pts*N_pts << std::endl;
			//out << "TENSORS stress double" << std::endl;
			//for(vect3d_map::iterator it=eval_points.begin();it!=eval_points.end(); ++it){
			//	out << sigma[it->first][0] << " " <<   tau[it->first][0] << " " <<   tau[it->first][1] << std::endl;
			//	out <<   tau[it->first][0] << " " << sigma[it->first][1] << " " <<   tau[it->first][2] << std::endl;
			//	out <<   tau[it->first][1] << " " <<   tau[it->first][2] << " " << sigma[it->first][2] << std::endl;
			//}
			out << "VECTORS tension double" << std::endl;
			for(int i=0;i<eval_points.size();++i){
				out << sigma[i][0] << " " << sigma[i][1] << " " << sigma[i][2] << std::endl;
			}
			out << "VECTORS shear double" << std::endl;
			for(int i=0;i<eval_points.size();++i){
				out <<   tau[i][0] << " " <<   tau[i][1] << " " <<   tau[i][2] << std::endl;
			}
			out << "VECTORS position double" << std::endl;
			for(int i=0;i<eval_points.size();++i){
				out << eval_points[i][0] << " " << eval_points[i][1] << " " << eval_points[i][2] << std::endl;
			}
			out.close();
		}
	}
}
