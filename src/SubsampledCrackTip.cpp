/* 
 * File:   SubsampledCrackTip.cpp
 * Author: David
 * 
 * Created on 24. Feb. 2014, 13:49
 */

#include "SubsampledCrackTip.h"
#include "VDBWrapper.h"
#include "FractureModel.h"
#include "PostProcessor.h"

#include <openvdb/util/Util.h> // needed only for vdb::util::INVALID_IDX

#include <cmath>
#include <sstream>
#include <fstream> // for debug output
#include <float.h>

namespace FractureSim{
	const double SubsampledCrackTip::collapseThr=0.2; // collapse if edge is shorter than collapseThr*meshSize
    const double SubsampledCrackTip::subdivThr  =1.4; // subdiv   if edge is longer  than subdivThr*meshSize

    // propagate all vertices for n substeps
    // store the crack-path in the level-set
    // and build the BEM mesh update
	int SubsampledCrackTip::propagate(
		node_map& newNodes, elem_map& newElems, id_map& newRegions,
		elem_map& newCrackTips, elem_map& newParents, state_map& newCrackTipStates,
		vect3d_map& crackTipFaceNormals, vect3d_map& crackTipTangents,
		vect3d_map& nodeSIFs, unsigned int steps
	){
        //updateVertexStates();
        //init(vertsPerSegment);
        movedVerts.clear();
		// store small scale mesh for level-set updates (built from propagating vertices)
        //store vertex positions, and which entry in the list is the first of the current ones
		std::vector<vdb::Vec3d> pointList;
		std::map<unsigned int, std::vector<vdb::Vec4I> > crackLists; // elements using quads!
		unsigned int pointListFirstCurrent=0;
        unsigned int pointListFirstPrevious;
		//pointList.clear(); crackLists.clear();
		//store initial vertex positions for level-set update
		for(vect3d_map::iterator it=verts.begin(); it!=verts.end(); ++it){
			pointList.push_back(vdb::Vec3d(it->second[0],it->second[1],it->second[2]));
		}

		// debug-output SIFs on node-verts
		addToSIFsOutput(nodeSIFs,"");

		// propagate vertices
        vect3d_map oldVerts; // store positions from the previous substep
        for(int step=0; step<steps; ++step){ // take substeps
			double stepsInv = 1.0/(double)steps , weight=0.0, old_th=0.0, old_th_min=0.0, old_th_max=0.0;
			double stepPart = (step+1.0)*stepsInv;
            //consitencyCheck();
            oldVerts = verts;
			//propagate each vertex
            for(vect3d_map::iterator it=verts.begin(); it!=verts.end(); ++it){
//#pragma omp parallel for schedule(static,128) // probably unsafe, also not really worth it in terms of runtime ...
//			for(int vi=verts.begin()->first; vi<=verts.rbegin()->first; ++vi){
//				vect3d_map::iterator it=verts.find(vi); if( it==verts.end() ) continue;
				//interpolate SIFs and coordinate system
                Eigen::Vector3d vertSIFs(0.0,0.0,0.0),n1(0.0,0.0,0.0),n2,n3(0.0,0.0,0.0),dn;
                for(int k=0; k<vertNodes[it->first].size();++k){
                    vertSIFs += nodeSIFs[vertNodes[it->first][k]]*nodeWeights[it->first][k];
                    n1+=crackTipFaceNormals[vertNodes[it->first][k]]*nodeWeights[it->first][k];
                    n3+=crackTipTangents[vertNodes[it->first][k]]*nodeWeights[it->first][k];
                }
                n1.normalize(); n3.normalize();
				n2=n3.cross(n1); n2.normalize(); // normalize should be obsolete

				if( FractureModel::cpVersion==2 ){
					// bool haveDir=false, haveMin=false, haveMax=false; // for debug only
					if( direction.count(it->first) ){
						weight = 1.0-stepPart;
						// weight*=weight; // quadratic blending is not so great for homogeneous materials
						old_th = atan2(n1.dot(direction[it->first]) , n2.dot(direction[it->first]));
						old_th_min = DBL_MAX; // set to invalid value so it won't be used if we don't have it
						old_th_max = DBL_MAX; // set to invalid value so it won't be used if we don't have it
						// haveDir=true;
					}
					if( dir_min.count(it->first) ){
						old_th_min = atan2(n1.dot(dir_min[it->first]) , n2.dot(dir_min[it->first]));
						// haveMin=true;
					}
					if( dir_max.count(it->first) ){
						old_th_max = atan2(n1.dot(dir_max[it->first]) , n2.dot(dir_max[it->first]));
						// haveMax=true;
					}
					if( old_th_min > (old_th_max+DBL_EPSILON) ){ // since the coordinate system changed we might need to correct by 2*pi
						if( old_th_min >= 0.0 ) old_th_min -= 2*M_PI;
						else if(old_th_max<0.0) old_th_max += 2*M_PI;
					}
					// if(haveDir && !(haveMin && haveMax) )
						// printf("\n !!! found old direction but not %s %s", haveMin?"":"min", haveMax?"":"max");
					// if( std::abs( old_th_min - old_th_max ) < 10*DBL_EPSILON )
						// printf("\nmin-max: ( %.6lf , %.6lf ) vert=%d w=%.4lf",old_th_min, old_th_max,it->first, weight);
				} else  weight = 0.0;

                //advance the current vertex
                if( states[it->first]!=INACTIVE &&
					fractureModel->fractureCriterion(dn,vertSIFs, it->second, n1,n2,n3, weight,old_th,old_th_min,old_th_max)
                ){
                    //printf("\n * %d --> (%.3lf, %.3lf, %.3lf)",it->first,dn[0],dn[1],dn[2]);
                    Eigen::Vector3d dx1, dx2; // transform dn to world coordinate displacements of endpoints
                    // dn[0] is radius of displacement, dn[1]=theta is angle of rotation around the segment's axis
                    // allow the crack segment to "twist" around the forward edge normal (n2) by angle phi=dn[2]

					
					if( FractureModel::cpVersion == 2 || direction.count(it->first)==0 ){
						dx1=(n2*cos(dn[1])+n1*sin(dn[1]));
					}else{ //if we have a previous direction, smoothly turn towards the new direction
						weight = (1.0-stepPart);
						if(FractureModel::cpVersion==1) weight *= weight;
						dx1=
							(1.0-weight)*(n2*cos(dn[1])+n1*sin(dn[1]))
							+weight*direction[it->first];
						dx1.normalize();
					}

					// store propagation direction on the last step
					if(step==(steps-1)){
						if( FractureModel::cpVersion == 2){ // new and experimental version
							double thH=0.0,thMin=0.0,thMax=0.0;
							thH=fractureModel->getHoopAngleAndInterval(vertSIFs,it->second,thMin,thMax);
							direction[it->first]=(n2*cos(thH)  +n1*sin(thH));
							dir_min  [it->first]=(n2*cos(thMin)+n1*sin(thMin));
							dir_max  [it->first]=(n2*cos(thMax)+n1*sin(thMax));
							//printf("\nmin-max: ( %.6lf , %.6lf )\n",thMin, thMax);
						}
						else direction[it->first]=dx1; // for versions 0 and 1 we simply store the last propagation direction
					}

                    dx1*= stepsInv*dn[0];
					// update position of vertex
                    it->second += dx1;
                    states[it->first]=UPDATE;
                    movedVerts.insert(it->first);

					// if the sample is close to the surface, also move toward the surface
					double lsval = levelSet->getValueAndGradient(it->second, dx2, ctRegions[it->first]);
					double gradnorm = dx2.norm(), bg=levelSet->getBandwidth()*levelSet->getVoxelSize();
					if( std::abs(lsval) < bg - FLT_EPSILON ){//gradnorm > FLT_EPSILON ){
						it->second += dx2/gradnorm*bg*stepsInv;/*levelSet->getVoxelSize()*(bg+lsval)/bg;*/
					}
                }
                // state-check based on level-set
                if( states[it->first]==UPDATE){
                    vdb::Vec3d
						v(it->second[0],it->second[1],it->second[2]),
						old_v(
						    oldVerts[it->first][0],
                            oldVerts[it->first][1],
                            oldVerts[it->first][2]);
                    states[it->first]=
                        levelSet->checkCrackState(v, old_v, ctRegions[it->first]);
                }
            } // end vertex loop
			crackTipSmoothing();

			// prepare level-set update:
			// add a quad connecting the previous position of this vertex and the current position of this vertex
			// to the previous and current positions of the left neighbor (if it exists)
			// since the map is key-sorted, and vertices are continously numbered starting at NODE_BASE_INDEX
			// each vertex will be at position i=vertex-ID - NODE_BASE_INDEX relative to the first vertex in the pointList
			// current vertex is at i+pointListFirstCurrent, its previous position at i+pointListFirstPrevious
			//store current vertex positions for level-set update
			pointListFirstPrevious=pointListFirstCurrent;
			pointListFirstCurrent = pointList.size();
			id_map consecutiveIndex; unsigned int k=0;
			for(vect3d_map::iterator it=verts.begin(); it!=verts.end(); ++it,++k){
				consecutiveIndex[it->first]=k; // VDB needs points to be indexed consecutively and 0-based
				pointList.push_back(vdb::Vec3d(it->second[0],it->second[1],it->second[2]));
			}
			vdb::Vec3d n,n2; // n is used to avoid updates to the level-set with 0-area elements
			unsigned int a,b,c,d; // tmp storage for point IDs
			for(vect3d_map::iterator it=verts.begin(); it!=verts.end(); ++it){
				a=pointListFirstPrevious + consecutiveIndex[it->first];
				b=pointListFirstPrevious + consecutiveIndex[neighbors[it->first][0]];
				c=pointListFirstCurrent  + consecutiveIndex[neighbors[it->first][0]];
				d=pointListFirstCurrent  + consecutiveIndex[it->first];

				n = ( pointList[b] - pointList[a] ).cross( pointList[c] - pointList[a] );
				n2= ( pointList[c] - pointList[a] ).cross( pointList[d] - pointList[a] );
				if((n.lengthSqr() !=0.0  ||
					n2.lengthSqr()!=0.0) && (
					movedVerts.count(it->first) !=0 ||
					movedVerts.count(neighbors[it->first][0]) !=0 ))
					crackLists[ctRegions[it->first]].push_back(vdb::Vec4I(a,b,c,d));

			}
            if(substepOutput){ // this makes the substepping slow, and self-intersection tests less reliable but can be used to generate intermediate output
                if((step<(steps-1)) && !outFileName.empty()){ // do not write an update file for the last substep as it will be followed by a full step output anyway
                    std::stringstream sstr;
                    sstr << outFileName << "_";
                    sstr.fill('0'); sstr.width(4);
                    sstr << step;
                    //levelSet->writeGrids(sstr.str());
                    levelSet->addCrack(pointList,crackLists,false,true,sstr.str());
                    vtkWrite(sstr.str() + "_ct.vtk");
                }else{
                    levelSet->addCrack(pointList,crackLists); 
                }
                crackLists.clear();
            }
			if( movedVerts.empty() ) step=steps; // terminate the loop since obviously no vertex has moved (since SIFs are constant, the next iteration won't move anything either)
       } // end substepping loop

        //update the level-set
		if(!substepOutput) levelSet->addCrack(pointList,crackLists);

		//printf("\nlsup_nodes=[\n"); //debug output
		//for(std::vector<vdb::Vec3d>::iterator n_it=pointList.begin(); n_it!=pointList.end(); ++n_it)
		//	printf("%.3le %.3le %.3le\n",(*n_it)[0],(*n_it)[1],(*n_it)[2]);
		//for(std::map<unsigned int, std::vector<vdb::Vec4I> >::iterator s_it=crackLists.begin(); s_it!=crackLists.end(); ++s_it){
		//	printf("];\nlsup_%d_elems=[\n",s_it->first);
		//	for(std::vector<vdb::Vec4I>::iterator e_it=s_it->second.begin(); e_it!=s_it->second.end(); ++e_it)
		//		printf("%d %d %d %d\n", (*e_it)[0],(*e_it)[1],(*e_it)[2],(*e_it)[3]);
		//}   printf("];\n%%");
        int ret = buildMeshUpdate(
            newNodes, newElems, newRegions, newCrackTips, newParents, newCrackTipStates
        );
        
//        printf("\n%% vertices by neighbor[0]");
//        unsigned int v=verts.begin()->first;
//        do{
//            printf("\n%%..v%3d(%+.3lf,%+.3lf,%+.3lf)l:%3d,r:%3d,s:%d,nd(%3d,%3d),w(%.3lf,%.3lf)",
//                    v, verts[v][0], verts[v][1], verts[v][2],
//                    neighbors[v][0],neighbors[v][1],states[v],
//                    vertNodes[v][0],vertNodes[v].size()>1?vertNodes[v][1]:-1,
//                    nodeWeights[v][0],nodeWeights[v].size()>1?nodeWeights[v][1]:0);
//            v=neighbors[v][0];
//        }while( v!=verts.begin()->first);
//        printf("\n%% vertices by neighbor[1]");
//        v=verts.begin()->first;
//        do{
//            printf("\n%%..v%3d(%+.3lf,%+.3lf,%+.3lf)l:%3d,r:%3d,s:%d,nd(%3d,%3d),w(%.3lf,%.3lf)",
//                    v, verts[v][0], verts[v][1], verts[v][2],
//                    neighbors[v][0],neighbors[v][1],states[v],
//                    vertNodes[v][0],vertNodes[v].size()>1?vertNodes[v][1]:-1,
//                    nodeWeights[v][0],nodeWeights[v].size()>1?nodeWeights[v][1]:0);
//            v=neighbors[v][1];
//        }while( v!=verts.begin()->first);
        
        return ret;
    }
    
    int SubsampledCrackTip::buildMeshUpdate(
		node_map& newNodes, elem_map& newElems, id_map& newRegions,
		elem_map& newCrackTips, elem_map& newParents, state_map& newCrackTipStates
    ){
        //for debugging ...
        //consitencyCheck();

        newNodes.clear();
        newElems.clear();
        newRegions.clear();
        newCrackTips.clear();
        newParents.clear();
        newCrackTipStates.clear();
        
		//first check for edge collapses
        //(before generating any new nodes - keeps the node-numbering simpler)
        collapseShortEdges();

        id_map addNodeVerts; // subdivisions will add new node verts (this is used as a reverse-map!)
        id_map newNodeIDs; // map old to new node IDs
		id_map newVertNodes; // store which verts have already spawned new nodes in order to avoid creating duplicate nodes
        unsigned int nextNode=(nodes.rbegin()->first)+1; // take the last (highest) node ID and add 1
        unsigned int nextElem=(elems.rbegin()->first)+1; // take the last (highest) element ID and add 1
        unsigned int nextCtElem=ELEM_BASE_INDEX;
        bool subdivided, collapsed;
        // loop over the old crack tip, check whether it has propagated and triangulate
        for(elem_map::const_iterator it = crackTips.begin(); it!= crackTips.end(); ++it){
            unsigned int  nd_a = it->second[0] , nd_b = it->second[1];
            unsigned int vertA= nodeVerts[nd_a],vertB = nodeVerts[nd_b];

            int newA=-1, newB=-1;
            subdivided = false;
            collapsed = false;

			if(vertA==vertB){ // detect collapsed edges
				collapsed=true;
				if(newNodeIDs.count(nd_a)!=0 && newNodeIDs.count(nd_b)!=0)
					if(newNodeIDs[nd_a]!=newNodeIDs[nd_b]){ //should be unreachable
                        printf("\nERROR: duplicate nodes found!!!\n");
                        newNodes.clear();
                        newElems.clear();
                        newRegions.clear();
                        newCrackTips.clear();
                        newParents.clear();
                        newCrackTipStates.clear();
                        return -1;
                    }
				if(newNodeIDs.count(nd_a)!=0) newNodeIDs[nd_b]=newNodeIDs[nd_a];
				if(newNodeIDs.count(nd_b)!=0) newNodeIDs[nd_a]=newNodeIDs[nd_b];
			}
			
			if(newVertNodes.count(vertA)==0){
	            if(movedVerts.count(vertA)!=0){
                    newNodes[nextNode].assign(3,0.0);
                    setNodeCoords(newNodes[nextNode],vertA);
                    newNodeIDs[nd_a]=nextNode;
					newVertNodes[vertA]=nextNode;
                    newA=nextNode; ++nextNode;
				}
            }else{
                newA=newVertNodes[vertA];
            }
			if(newVertNodes.count(vertB)==0){
				if(movedVerts.count(vertB)!=0){
					newNodes[nextNode].assign(3,0.0);
                    setNodeCoords(newNodes[nextNode],vertB);
					newNodeIDs[nd_b]=nextNode;
					newVertNodes[vertB]=nextNode;
					newB=nextNode; ++nextNode;
				}
			}else{
				newB=newVertNodes[vertB];
			}

			if(collapsed){ //just add the triangle connecting the old crack-tip edge to the point it has been collapsed into
				nextElem=addTri(
					newElems,newRegions,regions[crackTipParents[it->first][0]],
					newParents,nextElem,-1,nd_a,newA,nd_b
				);
			}else{ // not collapsed
				// build new crack tip element
				newCrackTips[nextCtElem].assign(2,0);
				newCrackTips[nextCtElem][0]= (newA!=-1)?newA:(nd_a);
				newCrackTips[nextCtElem][1]= (newB!=-1)?newB:(nd_b);

				if(newA==-1 && newB==-1){ //neither node propagated
					// simply copy state and parent info
					newParents[nextCtElem]=crackTipParents[it->first];
					newCrackTipStates[nextCtElem]=crackTipStates[it->first];
				}else{ // at least one node propagated
					//set state to update
					newCrackTipStates[nextCtElem]=updateState(crackTipStates[it->first]);
                    // subdivide edge if too long
                    subdivided=checkSubdivideEdge(nd_a, nd_b, newA, newB,
                        nextCtElem, nextNode, nextElem, newNodes, newElems,
                        newRegions, newCrackTips, newParents,
                        newCrackTipStates, addNodeVerts
                    );
				}
				// build triangulation behind crack tip
				if(newA==-1 && newB!=-1){ // nd_a not propagated, only "right" triangle
					nextElem=addTri(
						newElems,newRegions,regions[crackTipParents[it->first][0]],
						newParents,nextElem,nextCtElem,nd_a,newB,nd_b
					);
				}else if(newA!=-1 && newB==-1){ // b not propagated, only "left" triangle
					nextElem=addTri(
						newElems,newRegions,regions[crackTipParents[it->first][0]],
						newParents,nextElem,nextCtElem,newA,nd_b,nd_a
					);
				}else if(newA!=-1 && newB!=-1){ // both triangles
					nextElem=addTri(
						newElems,newRegions,regions[crackTipParents[it->first][0]],
						newParents,nextElem,nextCtElem,newA,newB,nd_a
					);
					nextElem=addTri(
						newElems,newRegions,regions[crackTipParents[it->first][0]],
						newParents,nextElem,-1,nd_a,newB,nd_b
					);
				}
				++nextCtElem;
				if(subdivided) subdivCleanup(
                    nextCtElem, nextElem, newNodes, newElems,
                    newCrackTips, newNodeIDs, newVertNodes, addNodeVerts
                );
			}
		}
        
        // update nodeVerts and vertNodes to access the new nodes
		// update mapping from new nodes to generating vertices
        id_map newNodeVerts;
        for(id_map::iterator it=nodeVerts.begin(); it!=nodeVerts.end(); ++it){
            if(newNodeIDs.count(it->first)!=0)
                newNodeVerts[newNodeIDs[it->first]]=it->second;
            else
                newNodeVerts[it->first]=it->second;
        }
        nodeVerts=newNodeVerts;
        for(id_map::iterator it=addNodeVerts.begin(); it!=addNodeVerts.end(); ++it){
            nodeVerts[it->second]=it->first;
//            printf("\n added nodeVert node %d vertex %d",it->second,it->first);
        }
//        printf("\nnew nodeVerts:");
//        for(id_map::iterator it=nodeVerts.begin(); it!=nodeVerts.end(); ++it){
//            printf("\n %d --> %d", it->first, it->second);
//        }        
		// update interpolation nodes for all vertices based on updated nodeVerts
		for(elem_map::const_iterator it = newCrackTips.begin(); it!= newCrackTips.end(); ++it){
			unsigned int nd_a = it->second[0]  , nd_b = it->second[1];
			if(nodeVerts.count(nd_a)==0 || nodeVerts.count(nd_b)==0)
				printf("\nERROR: no vertex assigned to node %s%s", (nodeVerts.count(nd_a)==0)?"a":"", (nodeVerts.count(nd_b)==0)?"b":"");
			unsigned int vertA= nodeVerts[nd_a], vertB= nodeVerts[nd_b];
			updateWeights(vertA,vertB,nd_a,nd_b);
		}
		// finally, update crack-tip states based on vertex states
		updateCrackTipStates(newCrackTips, newCrackTipStates);

		return 0;
	}

	int SubsampledCrackTip::init(double targetMeshSize, unsigned int samples, bool substepOutput_){
        meshSize = targetMeshSize;
        substepOutput=substepOutput_;
        if(samples < 1) return -1; // samples has to be at least 1 to make sense
        verts.clear();
        states.clear();
        neighbors.clear();
        vertNodes.clear();
        nodeWeights.clear();
        ctRegions.clear();
        nodeVerts.clear();

        unsigned int nextVert=NODE_BASE_INDEX;
        vertsPerSegment = samples;
        generateVertices(crackTips, nextVert); // generate verts for all crack tip elements
        
		return verts.size();
	}

	int SubsampledCrackTip::startCrack(
		elem_map& addElems, Eigen::Vector3d p, Eigen::Vector3d n1,
		Eigen::Vector3d n2, unsigned int region
	){
        double scale=collapseThr*meshSize*startScale;
		node_map addNodes;
        elem_map addCrackTips, addParents;
        id_map addRegions;
		state_map addStates;
		unsigned int n = nodes.rbegin()->first, // highest node
                     m = elems.rbegin()->first, // element
                     l;                         // and crack-tip numbers in mesh
        if(crackTips.empty()) l=ELEM_BASE_INDEX-1;
        else l = crackTips.rbegin()->first; 

		addElems.clear();
        /**
		// the new crack geometry will resemble a half-disc with 5 tris and 1 interior node
		// in MatLab notation:
		// t = 1/sqrt(2);
		// el=[ 1 2 3
		//		1 3 4
		//		1 4 5
		//		1 5 6
		//		1 6 2 ];
		// co=[  0 0.5 0
		//		-1   0 0
		//		 1   0 0
		//	 	 t   t 0
		//		 0   1 0
		//		-t   t 0 ];
		// M=[cross(n1,n2) n1 n2].*s;
		// co=(repmat(p,1,size(co,1))+M*co')';
		addElems[m+1].assign(3,0); addElems[m+1][0]=n+1; addElems[m+1][1]=n+2; addElems[m+1][2]=n+3;
		addElems[m+2].assign(3,0); addElems[m+2][0]=n+1; addElems[m+2][1]=n+3; addElems[m+2][2]=n+4;
		addElems[m+3].assign(3,0); addElems[m+3][0]=n+1; addElems[m+3][1]=n+4; addElems[m+3][2]=n+5;
		addElems[m+4].assign(3,0); addElems[m+4][0]=n+1; addElems[m+4][1]=n+5; addElems[m+4][2]=n+6;
		addElems[m+5].assign(3,0); addElems[m+5][0]=n+1; addElems[m+5][1]=n+6; addElems[m+5][2]=n+2;
		addRegions[m+1]=region;
		addRegions[m+2]=region;
		addRegions[m+3]=region;
		addRegions[m+4]=region;
		addRegions[m+5]=region;

		double t = 1.0/sqrt(2.0);  addNodes.clear();
		addNodes[n+1].assign(3,0.0);                        addNodes[n+1][1]=0.5;
		addNodes[n+2].assign(3,0.0); addNodes[n+2][0]=-1.0;
		addNodes[n+3].assign(3,0.0); addNodes[n+3][0]= 1.0;
		addNodes[n+4].assign(3,0.0); addNodes[n+4][0]=  t;  addNodes[n+4][1]= t;
		addNodes[n+5].assign(3,0.0);                        addNodes[n+5][1]=1.0;
		addNodes[n+6].assign(3,0.0); addNodes[n+6][0]= -t;  addNodes[n+6][1]= t;
		Eigen::Matrix3d rs; // rotation and scale matrix
		rs.col(0) = n1.cross(n2); // assuming n1,n2 are already normalized by caller
		rs.col(1) = n2;
		rs.col(2) = n1;
		rs*=scale;
		// apply rs and translation to p
		Eigen::Vector3d tmp;
		for(node_map::iterator it=addNodes.begin(); it!=addNodes.end(); ++it){
			tmp[0]=it->second[0]; tmp[1]=it->second[1]; tmp[2]=it->second[2];
			tmp = p + rs*tmp;
			it->second[0]=tmp[0]; it->second[1]=tmp[1]; it->second[2]=tmp[2];
		}
		// elems, nodes, and regions done, next up: crack-tips, parents and states
        addCrackTips.clear();
        addCrackTips[l+1].assign(2,0); addCrackTips[l+1][0]=n+2; addCrackTips[l+1][1]=n+3;
        addCrackTips[l+2].assign(2,0); addCrackTips[l+2][0]=n+3; addCrackTips[l+2][1]=n+4;
        addCrackTips[l+3].assign(2,0); addCrackTips[l+3][0]=n+4; addCrackTips[l+3][1]=n+5;
        addCrackTips[l+4].assign(2,0); addCrackTips[l+4][0]=n+5; addCrackTips[l+4][1]=n+6;
        addCrackTips[l+5].assign(2,0); addCrackTips[l+5][0]=n+6; addCrackTips[l+5][1]=n+2;
        addParents.clear();
        addParents[l+1].assign(2,0); addParents[l+1][0]=m+1;
        addParents[l+2].assign(2,0); addParents[l+2][0]=m+2;
        addParents[l+3].assign(2,0); addParents[l+3][0]=m+3;
        addParents[l+4].assign(2,0); addParents[l+4][0]=m+4;
        addParents[l+5].assign(2,0); addParents[l+5][0]=m+5;
        addStates.clear();
        addStates[l+1]=ACTIVE; //=INACTIVE; // this is the "back" edge of the crack - decided to init all ACTIVE since the BEM mesh may not coincide with the level-set surface
        addStates[l+2]=ACTIVE; //=ACTIVE_B; // this is "left-adjacent" to the back edge
        addStates[l+3]=ACTIVE; //=ACTIVE;
        addStates[l+4]=ACTIVE; //=ACTIVE;
        addStates[l+5]=ACTIVE; //=ACTIVE_A; // this is "right-adjacent" to the back edge
        /*/
		// different geometry: a hexagon with 6 tris and 1 interior node
		// in MatLab notation:
        //el=[ 1 2 3
        //     1 3 4
        //     1 4 5
        //     1 5 6
        //     1 6 7
        //     1 7 2 ];
        //co=[ 0 1 0
        //     sin(  pi/3) 1+cos(  pi/3) 0
        //     sin(2*pi/3) 1+cos(2*pi/3) 0
        //     sin(  pi  ) 1+cos(  pi  ) 0
        //     sin(4*pi/3) 1+cos(4*pi/3) 0
        //     sin(5*pi/3) 1+cos(5*pi/3) 0
        //     sin(2*pi  ) 1+cos(2*pi  ) 0 ];
		// M=[cross(n1,n2) n1 n2].*s;
		// co=(repmat(p,1,size(co,1))+M*co')';
		// ToDo: ? add different initial geometries, and an ENUM parameter to select one
		addElems[m+1].assign(3,0); addElems[m+1][0]=n+1; addElems[m+1][1]=n+2; addElems[m+1][2]=n+3;
		addElems[m+2].assign(3,0); addElems[m+2][0]=n+1; addElems[m+2][1]=n+3; addElems[m+2][2]=n+4;
		addElems[m+3].assign(3,0); addElems[m+3][0]=n+1; addElems[m+3][1]=n+4; addElems[m+3][2]=n+5;
		addElems[m+4].assign(3,0); addElems[m+4][0]=n+1; addElems[m+4][1]=n+5; addElems[m+4][2]=n+6;
		addElems[m+5].assign(3,0); addElems[m+5][0]=n+1; addElems[m+5][1]=n+6; addElems[m+5][2]=n+7;
		addElems[m+6].assign(3,0); addElems[m+6][0]=n+1; addElems[m+6][1]=n+7; addElems[m+6][2]=n+2;
		addRegions[m+1]=region; addRegions[m+2]=region;
		addRegions[m+3]=region; addRegions[m+4]=region;
		addRegions[m+5]=region; addRegions[m+6]=region;

		addNodes.clear();
		addNodes[n+1].assign(3,0.0);                                     addNodes[n+1][1]=1.0;
		addNodes[n+2].assign(3,0.0); addNodes[n+2][0]=sin(    M_PI/3.0); addNodes[n+2][1]=1.0+cos(    M_PI/3.0);
		addNodes[n+3].assign(3,0.0); addNodes[n+3][0]=sin(2.0*M_PI/3.0); addNodes[n+3][1]=1.0+cos(2.0*M_PI/3.0);
		addNodes[n+4].assign(3,0.0); addNodes[n+4][0]=sin(    M_PI    ); addNodes[n+4][1]=1.0+cos(    M_PI    );
		addNodes[n+5].assign(3,0.0); addNodes[n+5][0]=sin(4.0*M_PI/3.0); addNodes[n+5][1]=1.0+cos(4.0*M_PI/3.0);
		addNodes[n+6].assign(3,0.0); addNodes[n+6][0]=sin(5.0*M_PI/3.0); addNodes[n+6][1]=1.0+cos(5.0*M_PI/3.0);
		addNodes[n+7].assign(3,0.0); addNodes[n+7][0]=sin(2.0*M_PI    ); addNodes[n+7][1]=1.0+cos(2.0*M_PI    );
		Eigen::Matrix3d rs; // rotation and scale matrix
		rs.col(0) = n1.cross(n2); // assuming n1,n2 are already normalized by caller
		rs.col(1) = n2;
		rs.col(2) = n1;
		rs*=scale;
		// apply rs and translation to p
		Eigen::Vector3d tmp;
		for(node_map::iterator it=addNodes.begin(); it!=addNodes.end(); ++it){
			tmp[0]=it->second[0]; tmp[1]=it->second[1]; tmp[2]=it->second[2];
			tmp = p + rs*tmp;
			it->second[0]=tmp[0]; it->second[1]=tmp[1]; it->second[2]=tmp[2];
		}
		// elems, nodes, and regions done, next up: crack-tips, parents and states
        addCrackTips.clear();
        addCrackTips[l+1].assign(2,0); addCrackTips[l+1][0]=n+2; addCrackTips[l+1][1]=n+3;
        addCrackTips[l+2].assign(2,0); addCrackTips[l+2][0]=n+3; addCrackTips[l+2][1]=n+4;
        addCrackTips[l+3].assign(2,0); addCrackTips[l+3][0]=n+4; addCrackTips[l+3][1]=n+5;
        addCrackTips[l+4].assign(2,0); addCrackTips[l+4][0]=n+5; addCrackTips[l+4][1]=n+6;
        addCrackTips[l+5].assign(2,0); addCrackTips[l+5][0]=n+6; addCrackTips[l+5][1]=n+7;
        addCrackTips[l+6].assign(2,0); addCrackTips[l+6][0]=n+7; addCrackTips[l+6][1]=n+2;
        addParents.clear();
        addParents[l+1].assign(2,0); addParents[l+1][0]=m+1;
        addParents[l+2].assign(2,0); addParents[l+2][0]=m+2;
        addParents[l+3].assign(2,0); addParents[l+3][0]=m+3;
        addParents[l+4].assign(2,0); addParents[l+4][0]=m+4;
        addParents[l+5].assign(2,0); addParents[l+5][0]=m+5;
        addParents[l+6].assign(2,0); addParents[l+6][0]=m+6;
        addStates.clear();
        addStates[l+1]=ACTIVE; addStates[l+2]=ACTIVE;
        addStates[l+3]=ACTIVE; addStates[l+4]=ACTIVE;
        addStates[l+5]=ACTIVE; addStates[l+6]=ACTIVE;
        /**/
        // update level-set (ToDo: could start the level-set with a high-res half-disc?)
		levelSet->addCrack(addNodes,addElems,addRegions,n+1);
		
        // update BEM mesh
        nodes.insert(addNodes.begin(), addNodes.end());
        elems.insert(addElems.begin(), addElems.end());
        regions.insert(addRegions.begin(), addRegions.end());
        crackTips.insert(addCrackTips.begin(), addCrackTips.end());
        crackTipParents.insert(addParents.begin(), addParents.end());
        crackTipStates.insert(addStates.begin(), addStates.end());
        levelSet->addToNearTriGrid(nodes,addElems);

        //finally, subsample
        unsigned int vn;
        if(verts.empty()) vn=NODE_BASE_INDEX;
        else vn = verts.rbegin()->first+1;
        generateVertices(addCrackTips,vn,n+1);
		return 0;
	}

    int SubsampledCrackTip::vtkWrite(std::string fileName){
        id_map e1; id_set e2; elem_map e3; state_map e4; value_map e5; // empty maps to be given to the temp-postpro object
        node_map ctNodes; elem_map ctElems; unsigned int next=ELEM_BASE_INDEX;
		id_map consecutiveIndex; unsigned int k=NODE_BASE_INDEX;
		for(vect3d_map::iterator it=verts.begin(); it!=verts.end(); ++it,++k){
			consecutiveIndex[it->first]=k;
		} // VTK needs points to be indexed consecutively
        for(vect3d_map::iterator it=verts.begin(); it!=verts.end(); ++it){
            ctNodes[it->first].assign(3,0.0);
            ctNodes[it->first][0]=it->second[0];
            ctNodes[it->first][1]=it->second[1];
            ctNodes[it->first][2]=it->second[2];
            ctElems[next].assign(3,0);
            ctElems[next][0]=consecutiveIndex[it->first];
            ctElems[next][1]=consecutiveIndex[neighbors[it->first][0]];
            ctElems[next][2]=consecutiveIndex[neighbors[it->first][1]];
            ++next;
        }
        PostProcessor ctOut(ctNodes, ctElems, e1, e2, e3, e3, e4, e5);
        vector_type stateVector(ctNodes.size());
        int idx=0;
        for(state_map::iterator it=states.begin(); it!=states.end(); ++it, ++idx)
            stateVector[idx]=it->second;
        ctOut.vtkAddData("state",1,false,stateVector);
        return ctOut.vtkWrite(fileName,VTK_TRIANGLE);
    }
    
    int SubsampledCrackTip::generateVertices(
        elem_map& theseCrackTips, unsigned int nextVert, int sourceNode
    ){
        unsigned int firstNewVert=nextVert;
		for(elem_map::const_iterator it = theseCrackTips.begin();
			it!= theseCrackTips.end(); ++it
		){
            // create subsamples at the crack-tip nodes of the current segment
            for(int k=0; k<2; ++k) if(nodeVerts.count(it->second[k])==0 ){
                //k==0 the left node, k==1 is the right node
                verts[nextVert].setZero();
                verts[nextVert][0]=nodes[it->second[k]][0];
                verts[nextVert][1]=nodes[it->second[k]][1];
                verts[nextVert][2]=nodes[it->second[k]][2];
                states[nextVert]=INACTIVE; //initialize inactive
                nodeVerts[it->second[k]]=nextVert;
                neighbors[nextVert].assign(2,NODE_BASE_INDEX-1); //initialize with invalid index (might not have neighbors of node)
                vertNodes[nextVert].assign(1,it->second[k]);
                nodeWeights[nextVert].assign(1,1.0);
				ctRegions[nextVert]=regions[crackTipParents[it->first][0]];
				//printf("\n node  vert %3d has node (%3d) with weight (%.3lf)",
				//	nextVert, ctVertNodes[nextVert][0], ctVertNodeWeights[nextVert][0]
				//);
				++nextVert;
            }
            // initialize neighbors to nodes, ctVertNeighbors[...][0] is the left
            // and ...[1] is the right neighbor
            neighbors[nodeVerts[it->second[1]]][0] = nodeVerts[it->second[0]];
            neighbors[nodeVerts[it->second[0]]][1] = nodeVerts[it->second[1]];
            // now the neighbor-map points to vertices at the endpoints of the segment
            // later we only update when inserting intermediate vertices
            // activate vertices at nodes, if segment is active
            if(crackTipStates[it->first]!=INACTIVE){
                if(crackTipStates[it->first]!=ACTIVE_B)
                    states[nodeVerts[it->second[0]]]=ACTIVE;
                if(crackTipStates[it->first]!=ACTIVE_A)
                    states[nodeVerts[it->second[1]]]=ACTIVE;
            }

            // next: create subsamples along the segment
            if(vertsPerSegment>1 && (crackTipStates[it->first]!=INACTIVE) ){
                Eigen::Vector3d a,b;
                //postPro->copyNode(it->second[0], a);
				a[0]=nodes[it->second[0]][0];
				a[1]=nodes[it->second[0]][1];
				a[2]=nodes[it->second[0]][2];
                //postPro->copyNode(it->second[1], b);
				b[0]=nodes[it->second[1]][0];
				b[1]=nodes[it->second[1]][1];
				b[2]=nodes[it->second[1]][2];
                for(int k=0; k<=vertsPerSegment; ++k){
                    // k=0 corresponds to the vertex nodesToVerts[it->second[0]]
                    // k=samples takes the vertex at nodesToVerts[it->second[1]]
                    if(k==0){ // re-assign right neighbor of left node
                        neighbors[nodeVerts[it->second[0]]][1]=nextVert;
                    }else if(k==vertsPerSegment){ // re-assign left neighbor of right node
                        neighbors[nodeVerts[it->second[1]]][0]=nextVert-1;
                    }else{ // create subsample in between the nodes
						double k_s = (double)k/(double)vertsPerSegment;
						verts[nextVert]=(1.0-k_s)*a+k_s*b;
						//ctVertStates[nextVert]= (crackTipStates[it->first]!=INACTIVE) ? ACTIVE:INACTIVE;
						states[nextVert]=ACTIVE; // active check is done before creating subsamples now
						neighbors[nextVert].assign(2,nextVert-1); // usually the left neighbor is the previously created vertex
						neighbors[nextVert][1]=nextVert+1; // and the right neigh. is the one that will be created after the current one
						// special cases for neighbs.: first and last intermediate vert link to node-verts
						if(k==1)
							neighbors[nextVert][0] = nodeVerts[it->second[0]];
						if(k==vertsPerSegment-1)
							neighbors[nextVert][1] = nodeVerts[it->second[1]];

						vertNodes[nextVert].assign(2,it->second[0]); // link the vertex to the endpoints (nodes) of the current segment
						vertNodes[nextVert][1]=it->second[1];

						nodeWeights[nextVert].assign(2,1.0-k_s); // determine the interpolation weights of the nodes
						nodeWeights[nextVert][1]=k_s;
						ctRegions[nextVert]=regions[crackTipParents[it->first][0]];
						//printf("\n inter vert %3d has neighs (%3d,%3d) and nodes (%3d,%3d) with weights (%.3lf,%.3lf)",
						//	nextVert, ctVertNeighbors[nextVert][0],ctVertNeighbors[nextVert][1], ctVertNodes[nextVert][0],ctVertNodes[nextVert][1],
						//	ctVertNodeWeights[nextVert][0],ctVertNodeWeights[nextVert][1]
						//);
						++nextVert;
                    }
                }
            }
		}
        if(sourceNode>=NODE_BASE_INDEX){
            Eigen::Vector3d s(
                nodes[sourceNode][0],
                nodes[sourceNode][1],
                nodes[sourceNode][2]);
            for(unsigned int v=firstNewVert; v<nextVert; ++v){
                // if source node is given for these crack-tip vertices
                // initialze their direction to the vector vertex-sourceNode
                direction[v]=verts[v]-s;
                direction[v].normalize();
				if( FractureModel::cpVersion==2 ){ // also initialize their min and max directions ot the same vector
					dir_min[v]=direction[v];
					dir_max[v]=direction[v];
				}
            }
        }
        return 0;
    }

    
    void SubsampledCrackTip::crackTipSmoothing(){
		//updated Feb.2015: rewritten with delayed assignments
		double w=0.25; // smoothing weight for neighbor's values
		if( FractureModel::cpVersion==2 ) w=0.1;
		vect3d_map tmpPos, tmpDir, tmpMin, tmpMax;
        for(vect3d_map::iterator it=verts.begin(); it!=verts.end(); ++it){
			if( states[it->first]==INACTIVE &&
				movedVerts.count(neighbors[it->first][0])!=0 &&
				movedVerts.count(neighbors[it->first][1])!=0 &&
				states[neighbors[it->first][0]]==ACTIVE &&
				states[neighbors[it->first][1]]==ACTIVE
			){ states[it->first]==ACTIVE; } // sometimes there are false positives on the self-intersection test, try to fix them here
            if( states[it->first]!=INACTIVE && (
                movedVerts.count(neighbors[it->first][0])!=0 ||
                movedVerts.count(neighbors[it->first][1])!=0 )
            ){ // average with neighbors
                // it->second*=0.5;
                // it->second+=0.25*(verts[neighbors[it->first][0]]+verts[neighbors[it->first][1]]);
				tmpPos[it->first]=(1.0-2.0*w)*it->second+w*(verts[neighbors[it->first][0]]+verts[neighbors[it->first][1]]);

				// average directions if stored
				if( direction.count(it->first)>0 &&
					direction.count(neighbors[it->first][0])>0 &&
					direction.count(neighbors[it->first][1])>0 ){
					// direction[it->first]*=0.5;
					// direction[it->first]+=0.25*direction[neighbors[it->first][0]]+0.25*direction[neighbors[it->first][1]];
					// direction[it->first].normalize();
					tmpDir[it->first]=(1.0-2.0*w)*direction[it->first]+w*(direction[neighbors[it->first][0]]+direction[neighbors[it->first][1]]);
					tmpDir[it->first].normalize();
				}
				//also do this for min and max directions if present
				if( dir_min.count(it->first)>0 &&
					dir_min.count(neighbors[it->first][0])>0 &&
					dir_min.count(neighbors[it->first][1])>0 ){
					// dir_min[it->first]*=0.5;
					// dir_min[it->first]+=0.25*dir_min[neighbors[it->first][0]]+0.25*dir_min[neighbors[it->first][1]];
					// dir_min[it->first].normalize();
					tmpMin[it->first]=(1.0-2.0*w)*dir_min[it->first]+w*(dir_min[neighbors[it->first][0]]+dir_min[neighbors[it->first][1]]);
					tmpMin[it->first].normalize();
				}
				if( dir_max.count(it->first)>0 &&
					dir_max.count(neighbors[it->first][0])>0 &&
					dir_max.count(neighbors[it->first][1])>0 ){
					// dir_max[it->first]*=0.5;
					// dir_max[it->first]+=0.25*dir_max[neighbors[it->first][0]]+0.25*dir_max[neighbors[it->first][1]];
					// dir_max[it->first].normalize();
					tmpMax[it->first]=(1.0-2.0*w)*dir_max[it->first]+w*(dir_max[neighbors[it->first][0]]+dir_max[neighbors[it->first][1]]);
					tmpMax[it->first].normalize();
				}
            }
        }
		for(vect3d_map::iterator it=tmpPos.begin(); it!=tmpPos.end(); ++it){
			verts[it->first]=it->second;
		}
		for(vect3d_map::iterator it=tmpDir.begin(); it!=tmpDir.end(); ++it){
			direction[it->first]=it->second;
		}
		for(vect3d_map::iterator it=tmpMin.begin(); it!=tmpMin.end(); ++it){
			dir_min[it->first]=it->second;
		}
		for(vect3d_map::iterator it=tmpMax.begin(); it!=tmpMax.end(); ++it){
			dir_max[it->first]=it->second;
		}
        //make sure the crack-tip is nicely sampled
        for(elem_map::const_iterator it = crackTips.begin(); it!= crackTips.end(); ++it){
            unsigned int  nd_a = it->second[0] , nd_b = it->second[1];
            unsigned int vertA= nodeVerts[nd_a],vertB = nodeVerts[nd_b];
            redistributeVertices(vertA, vertB); //equalize sampling along the crack-tip
        }
    }
    void SubsampledCrackTip::redistributeVertices(
        unsigned int vertA, unsigned int vertB
    ){
		//return; // for testing
        // run along crack-tip by neigbor[1] starting at vertA
        // until we reach vertB
        // count how many vertices are in between (n)
        // compute the length of each line-segment
        // sum up to get the total arc-length from vertA to vertB
        // update the positions of all in-between vertices such that
        // they sit along the crack-tip in 1/n*length intervals
        unsigned int currentVertex= vertA, tmp;
		unsigned int n;
        double length=0.0;
		//printf("\nreD:vA=%d( ",vertA);
        for(n=0; n<vertsPerSegment && currentVertex!=vertB; ++n){
            tmp=currentVertex;
            currentVertex=neighbors[currentVertex][1];
            length+=(verts[currentVertex]-verts[tmp]).norm();
			//printf("%d ",currentVertex);
        }
		//printf("), vB=%d, n=%d",vertB,n);
        if(currentVertex!=vertB){
            printf("\nERROR: vertex b is not within vertsPerSegment hops of vertex a");
            return;
        }
        // now n holds the number of hops from vertA to vertB (there is 1 less
        // in-between vertex than hops), length holds the length of the path
        // next: update positions
        vect3d_map newVerts, newDirs,newMins,newMaxs;
		unsigned int k=1; // the next position will be at k/n*length
		unsigned int nextVertex=neighbors[vertA][1]; // which vertex is next to update position
        double s,stepLength, lastLength=0.0; // cumulate arc-length up the currentVertex
        currentVertex=vertA;
		//printf("\nreA: ");
        for(unsigned int i=0; i<n; ++i){
            tmp=currentVertex;
            currentVertex=neighbors[currentVertex][1];
            stepLength=(verts[currentVertex]-verts[tmp]).norm();
            
            // find which verts should get new positions along the segment (tmp, currentVertex)
			while(k<n && (double)k/n*length <= (lastLength+stepLength)){
				s=((double)k/n*length-lastLength)/stepLength;
				if(states[nextVertex]!=INACTIVE)
					newVerts[nextVertex]=(s*verts[currentVertex]+(1-s)*verts[tmp]);
                if( direction.count(currentVertex)!=0 &&
                    direction.count(tmp)          !=0)
                    newDirs[nextVertex]=
                        s*direction[currentVertex]+(1-s)*direction[tmp];
				//also do this for min and max directions if present
                if( dir_min.count(currentVertex)!=0 &&
                    dir_min.count(tmp)          !=0)
                    newMins[nextVertex]=
                        s*dir_min[currentVertex]+(1-s)*dir_min[tmp];
                if( dir_max.count(currentVertex)!=0 &&
                    dir_max.count(tmp)          !=0)
                    newMaxs[nextVertex]=
                        s*dir_max[currentVertex]+(1-s)*dir_max[tmp];
				//printf("%d(%d/%d,%ds%.3lf) ",nextVertex,k,n,i,s);
				nextVertex=neighbors[nextVertex][1];
				++k;
			}
            
            lastLength+=stepLength;
        }
		//update vertex positions (note can't use .insert(...) here, as that would not update existing elements)
		for(vect3d_map::iterator it=newVerts.begin(); it!=newVerts.end(); ++it){
			verts[it->first]=it->second;
		}
		for(vect3d_map::iterator it=newDirs.begin(); it!=newDirs.end(); ++it){
			direction[it->first]=it->second;
		}
		//also do this for min and max directions if present
		for(vect3d_map::iterator it=newMins.begin(); it!=newMins.end(); ++it){
			dir_min[it->first]=it->second;
		}
		for(vect3d_map::iterator it=newMaxs.begin(); it!=newMaxs.end(); ++it){
			dir_max[it->first]=it->second;
		}
    }
    
    int SubsampledCrackTip::collapseShortEdges(){
		id_map collapsedVertices;
		std::vector<vdb::Vec3d> pointList; //storage for level-set update
		std::map<unsigned int, std::vector<vdb::Vec4I> > crackLists;
		unsigned int pt=0; //number of entries in pointList

		for(elem_map::const_iterator it = crackTips.begin(); it!= crackTips.end(); ++it){
			unsigned int nd_a = it->second[0]  , nd_b = it->second[1];
			unsigned int vertA= nodeVerts[nd_a], vertB= nodeVerts[nd_b];
            if(vertA==vertB)
				continue; //already collapsed this ct elem
			if(verts.count(vertA)==0 || verts.count(vertB)==0)
				continue; //inconsistent node-vert mapping (will be fixed after the loop) - this can happen for adjacent collapses

			double length=(verts[vertB]-verts[vertA]).norm();
            //refuse collapse if resulting segments would be too long
            bool reject = false; bool foundNdVert=false;
            unsigned int v=vertB;
            for(int i=1; i<vertsPerSegment && !foundNdVert; ++i){
                v=neighbors[v][1];
                foundNdVert = ( vertNodes[v].size()!=2 );
            }
            if( (verts[v]-verts[vertA]).norm() > subdivThr*meshSize )
                reject=true; // reject if first neighboring segment would become too long
            v=vertA; foundNdVert=false;
            for(int i=1; i<vertsPerSegment && !foundNdVert; ++i){
                v=neighbors[v][0];
                foundNdVert = ( vertNodes[v].size()!=2 );
            }
            if( (verts[v]-verts[vertB]).norm() > subdivThr*meshSize )
                reject=true; // reject if second neighboring segment would become too long
			//if(reject) printf("\n*************************** collapse reject *************\n");


			if( !reject &&
				length < collapseThr*meshSize &&
				movedVerts.count(vertA)>0 && movedVerts.count(vertB)>0 // never collapse to a non-propagated node (don't change old tris!)
            ){  //printf("\n*(%d,%d)",nd_a,nd_b);
				// push old positions into pointList
				unsigned int v=vertA;
				for(int i=0; i<vertsPerSegment && v!=vertB; ++i){
					pointList.push_back(vdb::Vec3d(verts[v][0],verts[v][1],verts[v][2]));
					//printf("\n added pt i=%d, size is now %d, vB%d", i,pointList.size(), v==vertB);
					v=neighbors[v][1];
				} // when the loop finishes v==vertB should hold
				pointList.push_back(vdb::Vec3d(verts[vertB][0],verts[vertB][1],verts[vertB][2]));
				//printf("\n added pt vB%d, size is now %d", v==vertB,pointList.size());
				// remove all vertices on this element
				removeEdgeVerts(nd_a,nd_b);
				if( (states[vertA]!=INACTIVE && states[vertB]!=INACTIVE) ||
					(states[vertA]==INACTIVE && states[vertB]==INACTIVE) ){
					verts[vertA]*=0.5;
					verts[vertA]+=0.5*verts[vertB];
				}
				if( states[vertA]!=INACTIVE && states[vertB]==INACTIVE){
					//keep state and position of vertA
				}
				if( states[vertA]==INACTIVE && states[vertB]!=INACTIVE){
					//take state and position of vertB
					verts[vertA]=verts[vertB];
					states[vertA]=states[vertB];
				}
				//printf(" v%d=(%.1le,%.1le,%.1le)",vertA,verts[vertA][0],verts[vertA][1],verts[vertA][2]);
				// delete vertex
				deleteVertex(vertB);
				nodeVerts[nd_b]=vertA;
				//printf(" deleted %d: nodeVerts[%d]=%d ",vertB,nd_b,nodeVerts[nd_b]);
				collapsedVertices[vertB]=vertA; // map invalid vertex-ID to valid one
				//printf(" deleted: collapsedVertices[%d]=%d ",vertB,collapsedVertices[vertB]);
				//re-distributing vertexes on adjacent edges to vertA is done in buildMeshUpdate!
				// push new position into pointList
				pointList.push_back(vdb::Vec3d(verts[vertA][0],verts[vertA][1],verts[vertA][2]));
				//printf("\n added pt vNew, size is now %d",pointList.size());
				//add triangle (old_a,old_b,new_ab) to crackLists
				for(int i=0; i<vertsPerSegment; ++i){
					crackLists[regions[crackTipParents[it->first][0]]].push_back(vdb::Vec4I(
						pt+i,pt+vertsPerSegment+1,pt+i+1,vdb::util::INVALID_IDX //done: check that normals are fine here
					));
					//printf("\nadded tri (%d,%d,%d)",pt+i,pt+i+1,pt+vertsPerSegment+1);
				}
				pt=pointList.size();
			}
		}
		//for(id_map::iterator it=collapsedVertices.begin(); it!=collapsedVertices.end(); ++it)
		//	printf("\ncv: %d --> %d",it->first, it->second);

		//update nodeVerts in case we collapsed many adjacent edges into one vertex
		for(id_map::iterator it=nodeVerts.begin(); it!=nodeVerts.end(); ++it){
			while(collapsedVertices.count(it->second)!=0){
				it->second=collapsedVertices[it->second];
				//printf("\nnodeVerts[%d]=%d",it->first,it->second);
			}
		}

		//update the level-set such that edge-collapsing can't produce holes in crack surfaces
		levelSet->addCrack(pointList,crackLists,false,true); // don't update signed crack grid with these - normals are unreliable

		//ToDo? renumber vertices such that numbering is consecutive
		//done: the only occasions where we need consecutive numbering are
		//      the level-set update (because VDB uses a pointList rather than a map)
		//      and the VTK output (because VTK files assign node-IDs based on the line number in the file)
		//      these use a local remapping as workaround for non-consecutive numbering
		//ToDo: renumber if the highest vertex-ID is problematically high (like too close to int-max or some k*n_verts ?)
        return collapsedVertices.size();
    }
    
    bool SubsampledCrackTip::checkSubdivideEdge(
        unsigned int nd_a, unsigned int nd_b,
        int& newA, int& newB, unsigned int nextCtElem,
        unsigned int& nextNode, unsigned int& nextElem,
        node_map& newNodes, elem_map& newElems, id_map& newRegions,
        elem_map& newCrackTips, elem_map& newParents,
        state_map& newCrackTipStates, id_map& addNodeVerts
    ){
        // check if this ct element should be subdivided
        //double length=(verts[nodeVerts[nd_b]]-verts[nodeVerts[nd_a]]).norm();
		// use arc-length to check for subdivision
		unsigned int currentVertex= nodeVerts[nd_a], tmp;
		unsigned int n; bool haveCentre=false, reject=false;
        double length=0.0;
        for(n=0; n<vertsPerSegment && currentVertex!=nodeVerts[nd_b]; ++n){
            tmp=currentVertex;
            currentVertex=neighbors[currentVertex][1];
            length+=(verts[currentVertex]-verts[tmp]).norm();
            if( ((2*(n+1)) >= vertsPerSegment) && !haveCentre ){
                // just reached or crossed the half-way point
                haveCentre=true;
                Eigen::Vector3d oldA(nodes[nd_a][0],nodes[nd_a][1],nodes[nd_a][2]);
                Eigen::Vector3d oldB(nodes[nd_b][0],nodes[nd_b][1],nodes[nd_b][2]);
                // check that tri(nd_a,nd_b, tmp) has non-zero area
                double area= ( oldA-oldB ).cross( (verts[tmp]-oldB) ).norm();
                //printf("\nsubdiv area %.3le",area);
                reject = (area < FLT_EPSILON);
            }
        }
        
        //if(reject) printf("\nreject subdiv");
        //reject=true;
        if( !reject && //(newA!=-1 && newB!=-1) &&
            length > subdivThr*meshSize ){
            // increase sampling on this edge
            int centreVertex = increaseSampling(nd_a,nd_b);
            // add new node in the middle
            newNodes[nextNode].assign(3,0.0);
            setNodeCoords(newNodes[nextNode],centreVertex);
            if(newB!=-1){ // add tri (b,c,newB)
                // for a subdivision, we add an additional triangle
                // if nd_b has propagated, which will be the parent of the
                // second half of the crack-tip element we are subdividing
                // this tri should have a higher number than the ones added by
                // buildMeshUpdate to keep the element/node numbering consistent
                // with HyENA. There may be either 1 or 2 tris added later,
                // depending on whether nd_a has also propagated (newA!=-1)
                int addElems = (newA!=-1)?2:1;
                addTri( // nextElem must be updated in buildMeshUpdate later
                    newElems,newRegions, ctRegions[centreVertex],
                    newParents,nextElem+addElems,nextCtElem+1, newB,nd_b,nextNode
                );
            }else{ //newB == -1
                newParents[nextCtElem+1].assign(2,0);
                newParents[nextCtElem+1][0]=(newA==-1)? nextElem : (nextElem+1);
            }
            newB=nextNode; // restrict triangulation done in buildMeshUpdate
            newCrackTips[nextCtElem+1]=newCrackTips[nextCtElem];
            newCrackTips[nextCtElem  ][1]=nextNode;
            newCrackTips[nextCtElem+1][0]=nextNode;
            newCrackTipStates[nextCtElem+1]=newCrackTipStates[nextCtElem];

            // reverse-map the addition to nodeVerts so we can renumber nodes later
            addNodeVerts[centreVertex]=nextNode;
//            printf("\n addNodeVerts nd %d vert %d", nextNode, centreVertex);
			// rest to be done later in buildMeshUpdate
            ++nextNode;
            return true;
        }else
            return false; // no subdivision
    }

    void SubsampledCrackTip::subdivCleanup(
        unsigned int& nextCtElem, unsigned int& nextElem,
        node_map& newNodes, elem_map& newElems, elem_map& newCrackTips,
        id_map& newNodeIDs, id_map& newVertNodes, id_map& addNodeVerts
    ){
        // need to correct some numberings after subdivision
        ++nextCtElem; // subdivision added a second ct elem
        nextElem = newElems.rbegin()->first+1; // also may have added an extra triangle
        // it might still happen that the node-numbering is not HyENA compatible
        // so we'll check and fix it here ... (should only affect some subdivisions)
        id_map renumber;
        const unsigned int oldNds=(nodes.rbegin()->first);
        unsigned int next=oldNds+1;
        bool needUpdate=false;
        for(elem_map::iterator it=newElems.begin(); it!=newElems.end(); ++it){
            for(int k=0;k<3;++k){
                if(it->second[k]>oldNds){ // it's a new node
                    if(renumber.count(it->second[k])==0){
                        renumber[it->second[k]]=next;
                        if(next!=it->second[k]) needUpdate=true;
                        //printf("\n ndnr %d --> %d up%d",it->second[k],next,needUpdate);
                        ++next;
                    }
                }
            }
        }
        if(needUpdate){ //printf("\n UDPATE ND NR SUBDIV");
            // update ct elems
            for(elem_map::iterator it=newCrackTips.begin(); it!=newCrackTips.end(); ++it){
                for(int k=0;k<2;++k){
                    if(it->second[k]>oldNds){ // it's a new node
//                        printf("\n newCrackTips update %d --> %d", it->second[k],renumber[it->second[k]]);
                        it->second[k]=renumber[it->second[k]];
//                        printf(" is now (%d,%d):%d", it->second[0], it->second[1], it->second.size());
                    }
                }
            }
            // update elems
            for(elem_map::iterator it=newElems.begin(); it!=newElems.end(); ++it){
                for(int k=0;k<3;++k){
                    if(it->second[k]>oldNds){ // it's a new node
                        it->second[k]=renumber[it->second[k]];
                    }
                }
            }
            // update addNodeVerts
            for(id_map::iterator it=addNodeVerts.begin(); it!=addNodeVerts.end();++it){
//                printf("\n addNodeVerts update vert %d node %d --> %d",it->first, it->second,renumber[it->second]);
                it->second=renumber[it->second];
            }
            // update newNodeIDs
            for(id_map::iterator it=newNodeIDs.begin(); it!=newNodeIDs.end();++it){
//                printf("\n newNodeIDs update old-node %d node %d --> %d",it->first, it->second,renumber[it->second]);
                it->second=renumber[it->second];
            }
            // update newVertNodes
            for(id_map::iterator it=newVertNodes.begin(); it!=newVertNodes.end();++it){
//                printf("\n newVertNodes update vert %d node %d --> %d",it->first, it->second,renumber[it->second]);
                it->second=renumber[it->second];
            }
            // update nodes
            node_map tmpNewNodes(newNodes); //copy-construct
            newNodes.clear();
            for(node_map::iterator it=tmpNewNodes.begin();it!=tmpNewNodes.end();++it){
                newNodes[renumber[it->first]]=it->second;
            }
        }
    }

    void SubsampledCrackTip::removeEdgeVerts(
        unsigned int nd_a, unsigned int nd_b
    ){
        unsigned int vertA=nodeVerts[nd_a], vertB=nodeVerts[nd_b];
        unsigned int currentVertex= neighbors[vertA][1];
		unsigned int tmp;
        for(int i=1; i<vertsPerSegment && currentVertex!=vertB; ++i){
			tmp=neighbors[currentVertex][1];
			// delete vertex
			deleteVertex(currentVertex);
            // next iteration
            currentVertex= tmp;
        }

    }

	void SubsampledCrackTip::deleteVertex(unsigned int v){
		unsigned int n0=neighbors[v][0];
		unsigned int n1=neighbors[v][1];
		neighbors[n0][1]=n1;
		neighbors[n1][0]=n0;
		verts      .erase(v); //printf("\ndel vert %d",v);
		states     .erase(v);
		neighbors  .erase(v);
		vertNodes  .erase(v);
		nodeWeights.erase(v);
		ctRegions .erase(v);	
	}

    unsigned int SubsampledCrackTip::increaseSampling(
        unsigned int nd_a, unsigned int nd_b
    ){
        unsigned int vertA=nodeVerts[nd_a], vertB=nodeVerts[nd_b];
		//printf("\n|(%d,%d) v(%d,%d)",nd_a,nd_b,vertA,vertB);
        //state of new vertices
        //CRACK_STATE newState=INACTIVE;
        //if( states[vertA]==ACTIVE &&
        //    states[vertB]==ACTIVE) newState=ACTIVE;
        unsigned int newRegion = ctRegions[vertA];
        unsigned int currentVertex= vertA;
        unsigned int nextVertex   = neighbors[currentVertex][1]; //right neighbor
        unsigned int nextID = verts.rbegin()->first+1;
        unsigned int centreVertex=0; bool haveCentre=false;
        for(int i=0; i<vertsPerSegment; ++i){
            // add a new vertex between current and next
            verts[nextID]=0.5*(verts[currentVertex]+verts[nextVertex]);
            if( direction.count(currentVertex)!=0 &&
                direction.count(nextVertex)   !=0)
                direction[nextID]=
                    0.5*(direction[currentVertex]+direction[nextVertex]);
			//also do this for min and max directions if present
            if( dir_min.count(currentVertex)!=0 &&
                dir_min.count(nextVertex)   !=0)
                dir_min[nextID]=
                    0.5*(dir_min[currentVertex]+dir_min[nextVertex]);
			if( dir_max.count(currentVertex)!=0 &&
                dir_max.count(nextVertex)   !=0)
                dir_max[nextID]=
                    0.5*(dir_max[currentVertex]+dir_max[nextVertex]);
			//states[nextID]=newState;
			if( states[currentVertex]==ACTIVE ||
				states[   nextVertex]==ACTIVE) states[nextID]=  ACTIVE;
			else                               states[nextID]=INACTIVE;
            ctRegions[nextID]=newRegion;
            vertNodes[nextID].assign(2,nd_a);
            vertNodes[nextID][1]=nd_b;
            nodeWeights[nextID].assign(2,0.0);
            // update neighbors
            neighbors[currentVertex][1]=nextID;
            neighbors[nextVertex][0]=nextID;
            neighbors[nextID].assign(2,currentVertex);
            neighbors[nextID][1]=nextVertex;
            
            // find centre vertex
            if( ((2*(i+1)) >= vertsPerSegment) && !haveCentre ){
                // just reached or crossed the half-way point
                if( (vertsPerSegment%2)==0 ){
                    // vertsPerSegment is even --> pick the next vertex
                    centreVertex=nextVertex;
                    //printf("\n%%ctr is next at i=%d, vs=%d",i,vertsPerSegment);
                }else{
                    // vertsPerSegment is odd  --> pick the new vertex
                    centreVertex=nextID;
                    //printf("\n%%ctr is inter at i=%d, vs=%d",i,vertsPerSegment);
                }
                haveCentre=true;
            }
            
            // next iteration
            currentVertex= nextVertex;
            nextVertex   = neighbors[nextVertex][1];
            ++nextID;
        }
        return centreVertex;
    }

    void SubsampledCrackTip::updateWeights(
        unsigned int vertA, unsigned int vertB, 
        unsigned int nodeA, unsigned int nodeB 
    ){
		vertNodes[vertA].assign(1,nodeA);
		vertNodes[vertB].assign(1,nodeB);
		nodeWeights[vertA].assign(1,1.0);
		nodeWeights[vertB].assign(1,1.0);
        unsigned int v=neighbors[vertA][1];
        if(v==vertB) return; // segment may have no interior nodes (inactive)
        for(int k=1; k<vertsPerSegment; ++k){
            double k_s = (double)k/(double)vertsPerSegment;
            vertNodes[v].assign(2,nodeA); // link the vertex to the endpoints (nodes) of the current segment
            vertNodes[v][1]=      nodeB ;
            nodeWeights[v].assign(2,1.0 - k_s); // determine the interpolation weights of the nodes
            nodeWeights[v][1]=            k_s ;

            v=neighbors[v][1];
        }
		if(v!=vertB) printf("\npossibly inconsistent interpolation between nodes %d and %d, verts %d, %d, stopped at %d",nodeA,nodeB, vertA,vertB,v);
    }
    
    int SubsampledCrackTip::addTri(
        elem_map& newElems, id_map& newRegions, unsigned int regionID,
		elem_map& newParents, int nextElem, int ctElem,
        unsigned int nd_a, unsigned int nd_b, unsigned int nd_c
    ){
        newElems[nextElem].assign(3,0);
        newElems[nextElem][0]=nd_a;
        newElems[nextElem][1]=nd_b;
        newElems[nextElem][2]=nd_c;
        newRegions[nextElem]=regionID;
        if(ctElem>=0){ // this tri is a parent of a crack-tip element
            newParents[ctElem].assign(2,0); // keep Elmer format intact (each edge can have two parents)
            newParents[ctElem][0]=nextElem;
        }
        return nextElem+1;
    }

	void SubsampledCrackTip::updateCrackTipStates(elem_map& newCrackTips, state_map& newCrackTipStates){
		// run through newCrackTips in update state, fetch node-verts, check all interior verts
		// if all verts are active, the segment is active
		// if all verts are inactive, the segment is inactive
		// if any but not all verts are active, the segment is active_a
		for(elem_map::iterator it=newCrackTips.begin(); it!=newCrackTips.end(); ++it){
			if( newCrackTipStates[it->first]==UPDATE   ||
				newCrackTipStates[it->first]==UPDATE_A ||
				newCrackTipStates[it->first]==UPDATE_B ||
				newCrackTipStates[it->first]==ACTIVE_A ){
				unsigned int vertA=nodeVerts[it->second[0]];
				unsigned int vertB=nodeVerts[it->second[1]];

				CRACK_STATE newState=states[vertB];
				unsigned int v=vertA;
				for(int i=0; i<vertsPerSegment && v!=vertB && newState!=ACTIVE_A; ++i){
					if(newState==ACTIVE && states[v]==INACTIVE)
						newState=ACTIVE_A;
					if(newState==INACTIVE && states[v]==ACTIVE)
						newState=ACTIVE_A;

					v=neighbors[v][1];
				}
//                if(newCrackTipStates[it->first]==ACTIVE_A)
//                    printf("\nswitching from active_a to %s",
//                        (newState==INACTIVE)?"INACTIVE":"stuff");
				newCrackTipStates[it->first]=newState;
			}
		}
	}
        
    void SubsampledCrackTip::updateVertexStates(){
        // old version: copy states from coarse mesh
		printf("\nOLD VERSION - DON'T USE (SubsampledCrackTip::updateVertexStates)\n");
		for(elem_map::const_iterator it = crackTips.begin();
			it!= crackTips.end(); ++it
		){
            CRACK_STATE newState=INACTIVE;
            // update state of node-vertices
            states[nodeVerts[it->second[0]]]=INACTIVE;
            states[nodeVerts[it->second[1]]]=INACTIVE;
            if(crackTipStates[it->first]!=INACTIVE){
                if(crackTipStates[it->first]!=ACTIVE_B){
                    states[nodeVerts[it->second[0]]]=ACTIVE;
                    newState=ACTIVE;
                }
                if(crackTipStates[it->first]!=ACTIVE_A){
                    states[nodeVerts[it->second[1]]]=ACTIVE;
                    newState=ACTIVE;
                }
            }
            // update state of intermediate vertices
            // start with next neighbor of the first node
            int vertID=neighbors[nodeVerts[it->second[0]]][1];
//            for(int i=0; i<vertsPerSegment; ++i){
            //printf("\nsegment %d new state %d at ",it->first,newState);
            while(vertID!=nodeVerts[it->second[1]]){
                // loop trough all vertices on the current segment
                states[vertID]=newState;
                //printf("%d ",vertID);
                vertID=neighbors[vertID][1];
            }
        }
    }
    
    void SubsampledCrackTip::setNodeCoords(std::vector<double>& c, unsigned int vertex){
        // simple version: just copy from given vertex
        c[0]=verts[vertex][0];
        c[1]=verts[vertex][1];
        c[2]=verts[vertex][2];
        //return;
        // smoother version: average a few of the vertex' neighbors, unless they are inactive
        unsigned int m=1, n = 1+(vertsPerSegment-1)/2; // int div!
        unsigned int v0,v1;
        bool v0active=true, v1active=true;
        v0 = neighbors[vertex][0];
        v1 = neighbors[vertex][1];
        for(int k=0; k<n && (v0active || v1active); ++k){
            if(states[v0]==INACTIVE) v0active=false;
            if(states[v1]==INACTIVE) v1active=false;
            if(v0active){
                ++m;
                c[0]+= verts[v0][0];
                c[1]+= verts[v0][1];
                c[2]+= verts[v0][2];
                v0 = neighbors[v0][0];
            }
            if(v1active){
                ++m;
                c[0]+= verts[v1][0];
                c[1]+= verts[v1][1];
                c[2]+= verts[v1][2];
                v1 = neighbors[v1][1];
            }
        }
        c[0]/=m; c[1]/=m; c[2]/=m;
    }
    
    void SubsampledCrackTip::consitencyCheck(){
 		//debug: consistency check

		id_set ctNodes = nodeSet(crackTips);
		for(vect3d_map::iterator it = verts.begin(); it != verts.end(); ++it){
			// MSVC has no std::isnan (?)
            //if(std::isnan(it->second[0]) || std::isnan(it->second[1]) || std::isnan(it->second[2]))
            //    printf("\nvertex %d has a NaN coordinate (%.2le,%.2le,%.2le)",it->first,
            //        it->second[0],it->second[1],it->second[2]);
            
			if(states.count(it->first) ==0) printf("\nvertex %d has no state",it->first);
			if(neighbors.count(it->first) ==0) printf("\nvertex %d has no neighbors",it->first);
			else if(neighbors[it->first].size() !=2) printf("\nvertex %d has %d neighbors",it->first,neighbors[it->first].size());
			if(ctRegions.count(it->first)==0) printf("\nvertex %d has no region",it->first);

			int k = nodeWeights[it->first].size();
			if( k!= vertNodes[it->first].size())
				printf("\nvertex %d has %d interpolation nodes but %d weights",it->first,vertNodes[it->first].size(),nodeWeights[it->first].size());
			double s=0.0;
			for(int i=0; i<k; ++i){
				s+=nodeWeights[it->first][i];
				if( ctNodes.count( vertNodes[it->first][i] )==0)
					printf("\nvertex %d has interpolation node %d which is not on the crack-tip!",it->first,vertNodes[it->first][i]);
			}
			if(std::abs(1.0-s)>DBL_EPSILON) printf("\nsum of weights at vertex %d is %.3lf", it->first, s);
		}

		for(elem_map::const_iterator it = crackTips.begin(); it!= crackTips.end(); ++it){
            unsigned int nd_a, nd_b;
            nd_a = it->second[0]; nd_b = it->second[1];
			if(nodeVerts.count(nd_a)>0){
				int v=nodeVerts[nd_a];
				if(vertNodes.count(v)>0){
					int n=vertNodes[v][0];
					if(vertNodes[v].size()>1) printf("node vert has too many interpolation nodes");
					if(n!=nd_a) printf("vertNodes<->nodeVerts inconsistent\nnd_a=%d, v[nd_a]=%d, n[v]=%d\n",nd_a,v,n);
				}
			}
			if(nodeVerts.count(nd_b)>0){
				int v=nodeVerts[nd_b];
				if(vertNodes.count(v)>0){
					int n=vertNodes[v][0];
					if(vertNodes[v].size()>1) printf("node vert has too many interpolation nodes");
					if(n!=nd_b) printf("vertNodes<->nodeVerts inconsistent\nnd_b=%d, v[nd_b]=%d, n[v]=%d\n",nd_b,v,n);
				}
			}
			bool aValid=true, bValid=true;
			if(aValid && nodeVerts.count(nd_a)>0); else{ aValid=false; printf("nd %d has no nodeVert\n",nd_a);}
			if(aValid && verts.count(nodeVerts[nd_a])>0); else{ aValid=false; printf("vertex %d is invalid\n",nodeVerts[nd_a]);}
			if(aValid && neighbors.count(nodeVerts[nd_a])>0); else{ aValid=false; printf("vertex %d (nd %d) has no neighbors\n",nodeVerts[nd_a],nd_a);}
			if(bValid && nodeVerts.count(nd_b)>0); else{ bValid=false; printf("nd %d has no nodeVert\n",nd_b);}
			if(bValid && verts.count(nodeVerts[nd_b])>0); else{ bValid=false; printf("vertex %d is invalid\n",nodeVerts[nd_b]);}
			if(bValid && neighbors.count(nodeVerts[nd_b])>0); else{ bValid=false; printf("vertex %d (nd %d) has no neighbors\n",nodeVerts[nd_b],nd_b);}
			if(!aValid || !bValid){
				printf("\n error on ct-elem (%d,%d) with nodeVerts (%d,%d)",
					nd_a,nd_b,nodeVerts[nd_a],nodeVerts[nd_b]);
			}
		}
    }

	int SubsampledCrackTip::getCoarseCrackFrontStates(state_map& coarseStates){
		int activeCount=0;
		coarseStates.clear();
		CRACK_STATE currentEdgeState;

		for(elem_map::const_iterator it = crackTips.begin(); it!= crackTips.end(); ++it){
			currentEdgeState = INACTIVE;
            unsigned int nd_a, nd_b;
            nd_a = it->second[0]; nd_b = it->second[1];
			unsigned int v=nodeVerts[nd_a];

			for(int k=0; k<vertsPerSegment && v!=nodeVerts[nd_b] && currentEdgeState==INACTIVE; ++k){
				if( states[v]!=INACTIVE ) currentEdgeState=ACTIVE;
				v=neighbors[v][1];
			}
			if( currentEdgeState==ACTIVE ){
				if( coarseStates.count(nd_a)==0 || coarseStates[nd_a]==INACTIVE ) ++activeCount; // only count the first occurence
				if( coarseStates.count(nd_b)==0 || coarseStates[nd_b]==INACTIVE ) ++activeCount;
				coarseStates[nd_a]=ACTIVE;
				coarseStates[nd_b]=ACTIVE;
			}else{
				if( coarseStates.count(nd_a)==0 ) coarseStates[nd_a]=INACTIVE; // do not overwrite here
				if( coarseStates.count(nd_b)==0 ) coarseStates[nd_b]=INACTIVE;
			}
		}
		return activeCount;
	}
	
	int SubsampledCrackTip::copyEdges(SubsampledCrackTip& source, edge_imap sourceEdges, id_map targetNodes, bool resetState){
		int copyCount=0;
 
		// for each edge in sourceEdges
		// -- check if it exists in source
		// -- check if corresponding edge exists in this object (fwd. mapped through targetNodes)
		// -- if both found
		// -- check that both source and target edges are properly sampled
		// -- copy marker positions, states, directions, etc.
		
		// for now we require that both source and target (this) have same vertsPerSegment
		// -- otherwise (optional/for future work - not implemented yet)
		// -- we could interpolate marker data along source curve
		if( source.vertsPerSegment != vertsPerSegment ){
			printf("\nSubsampledCrackTip::copyEdges: source and target objects must use the same number of markers per edge\n");
			return -1;
		}

		// edge-sampling check:
		//unsigned int vertA=nodeVerts[nd_a], vertB=nodeVerts[nd_b];
		//unsigned int v=neighbors[vertA][1];
		//if(v==vertB) ... ; // segment may have no interior nodes (inactive)
		//for(int k=1; k<vertsPerSegment; ++k) v=neighbors[v][1];
		//if(v!=vertB) printf("\npossibly inconsistent interpolation between nodes %d and %d, verts %d, %d, stopped at %d",nodeA,nodeB, vertA,vertB,v);

		for(edge_imap::iterator ed_it=sourceEdges.begin(); ed_it!=sourceEdges.end(); ++ed_it){
			unsigned int
				srcA=ed_it->first.first,
				srcB=ed_it->first.second;
			unsigned int
				tgtA=targetNodes[srcA],
				tgtB=targetNodes[srcB];
			unsigned int vsa,vsb,vta,vtb; // vertex ids for source and targes nodes
			bool haveSrcEdg=false, haveTgtEdg=false;
			// check source edge
			if( source.nodeVerts.count(srcA)>0 &&
				source.nodeVerts.count(srcB)>0
			){ // have verts for source nodes
				vsa=source.nodeVerts[srcA],
				vsb=source.nodeVerts[srcB];
				unsigned int v=source.neighbors[vsa][1];
				//if(v==vertB) ... ; // segment may have no interior nodes (inactive)
				for(int k=1; k<source.vertsPerSegment; ++k) v=source.neighbors[v][1];
				if(v==vsb) haveSrcEdg=true; // found properly sampled source edge
			}
			// check target edge
			if( haveSrcEdg &&
				nodeVerts.count(tgtA)>0 &&
				nodeVerts.count(tgtB)>0
			){ // have verts for target nodes
				vta=nodeVerts[tgtA],
				vtb=nodeVerts[tgtB];
				unsigned int v=neighbors[vta][1];
				//if(v==vertB) ... ; // segment may have no interior nodes (inactive)
				for(int k=1; k<vertsPerSegment; ++k) v=neighbors[v][1];
				if(v==vtb) haveTgtEdg=true; // found properly sampled source edge
			}
			if( haveSrcEdg && haveTgtEdg ){
				// now we can copy marker data ...
				// we assume vertsPerSegment==source.vertsPerSegment !!!
				unsigned int vs=vsa,vt=vta;
				for(int k=0; k<=vertsPerSegment; ++k){
					// copy data from source vertex (vs) to target (vt) ...
					if( source.verts.count(vs)>0 )
						verts[vt] = source.verts[vs];
					if( resetState ) states[vt] = ACTIVE;
					else if( source.states.count(vs)>0 )
						states[vt] = source.states[vs];
					if( source.direction.count(vs)>0 )
						direction[vt]=source.direction[vs];
					if( source.dir_min.count(vs)>0 )
						dir_min[vt]=source.dir_min[vs];
					if( source.dir_max.count(vs)>0 )
						dir_max[vt]=source.dir_max[vs];

					vs=source.neighbors[vs][1];
					vt=neighbors[vt][1];
				}
				if( resetState ){
					redistributeVertices(vta,vtb);
					// reset the state after redistributing the vertices
					vt=vta;
					for(int k=0; k<=vertsPerSegment; ++k){
						states[vt] = levelSet->isInside(verts[vt])?ACTIVE:INACTIVE;
						vt=neighbors[vt][1];
					}
				}
				++copyCount;
			}
		}
		return copyCount;
	}

	void SubsampledCrackTip::setSIFsOutputFile(std::string filename){
		if( sifOut!=NULL){
			sifOut->close();
			delete sifOut;
			sifOut=NULL;
		}
		if(!filename.empty()){
			sifOut = new std::ofstream((filename+".m").c_str());
			sifOutSteps.clear(); sifOutValues.clear();
			if(!sifOut->is_open()){ delete sifOut; sifOut=NULL; }
		}
	}
	void SubsampledCrackTip::addToSIFsOutput(vect3d_map& nodeSIFs, std::string prefix){
		unsigned int i;
		if( sifOut!=NULL && sifOut->good()){
			if( sifOutSteps.count(prefix)==0 )  sifOutSteps[prefix]=0;
			if( sifOutValues.count(prefix)==0 ) sifOutValues[prefix]=0;
			++sifOutSteps[prefix];
			//printf("\n%% addToSIFsOutput %s%d ... ", prefix.c_str(), sifOutSteps[prefix]);
			//(*sifOut) << "step " << sifOutCount << std::endl;
			vect3d_map tmp; // map vert-id to sifs so we can sort output by vert-id
			for(id_map::iterator it=nodeVerts.begin(); it!=nodeVerts.end(); ++it){
				if(states[it->second]==ACTIVE && nodeSIFs.count(it->first)>0)
					tmp[it->second]=nodeSIFs[it->first];
			}
			for(vect3d_map::iterator it=tmp.begin(); it!=tmp.end(); ++it){ i=++sifOutValues[prefix];
				// write sifs and vert id to file
				//(*sifOut) << it->first << ": (" << it->second[0] << "," << it->second[1] << "," << it->second[2] << ")" << std::endl;
				(*sifOut) << prefix << "i(" << i << ")=" << sifOutSteps[prefix] << "; " << prefix << "j(" << i << ")=" << it->first << "; ";
				(*sifOut) << prefix << "k1("<< i << ")=" << it->second[0] << "; "
						  << prefix << "k2("<< i << ")=" << it->second[1] << "; "
						  << prefix << "k3("<< i << ")=" << it->second[2] << ";" << std::endl;
			}
		}
	}

	SubsampledCrackTip::~SubsampledCrackTip(){
		if(sifOut!=NULL){
			sifOut->close();
			delete sifOut;
			sifOut=NULL;
		}
	}
}
