#include "PostProcessor.h"
#include "hyena/wrap/QuadratureWrapper.h"
#include "distPtToTri.h"

#include <fstream>
#include <cmath>
#include <cfloat>

#include <iostream>

using namespace std;

namespace FractureSim{

	int PostProcessor::computeNodeSIFs(
            const vector_type& displacements, double E, double nu,
			vect3d_map& sifs, vect3d_map& faceNormals, vect3d_map& tangents
    ){
		// compute local coordinate systems and SIFs at all nodes on a crack-tip
		// uses information from the crack-tip parent elements (not necessarily all adjacent nodes)
		// - run along the crack-tip, for each segment, take the parent triangle
		// - compute the averaged face normal of parent-tris at each endpoint
		// - compute the averaged SIFs by displacement correlation from the interior nodes
		sifs.clear();
		getCrackTipAxis(faceNormals,tangents);
		// now we have surface normals and crack-tip tangents at all crack-tip nodes
		// next: use displacement correlation to compute SIFs
		id_map count; //printf("\n%% computing SIFs... ");
		for(elem_map::iterator it=crackTips.begin(); it!=crackTips.end(); ++it){
            unsigned int nd_a, nd_b, nd_c; // nodes of the current triangle
            nd_a = it->second[0]; nd_b = it->second[1];
            nd_c = findInteriorNode(nd_a, nd_b, elems[parents[it->first][0]]);
			Eigen::Vector3d a,b,c, n1,n2,n3, cod; // crack opening displacement at correlation point
			copyNode(nd_a, a);
            copyNode(nd_b, b);
            copyNode(nd_c, c);
            getLocalCoordFrame(a,b,c, n1,n2,n3);
            // correlation point is the interior node (nd_c)
            cod[0] = displacements[3*(nd_c-NODE_BASE_INDEX)  ];
            cod[1] = displacements[3*(nd_c-NODE_BASE_INDEX)+1];
            cod[2] = displacements[3*(nd_c-NODE_BASE_INDEX)+2];
			double r; // distance from the crack tip to the correlation point
            r=((0.5*(a+b))-c).dot(n2); //project midpoint-interior to edge normal (shortest distance from correlation pt to crack tip)
            double mu = E/(2*(1+nu));
            double u1,u2,u3; // projected cod in n1,n2,n3 directions
			
            u1=cod.dot(faceNormals[nd_a]);
			u2=cod.dot(tangents[nd_a].cross(faceNormals[nd_a]));
            u3=cod.dot(tangents[nd_a]);
			if(sifs.count(nd_a)==0) sifs[nd_a].setZero();
            sifs[nd_a][0] += -mu*sqrt(2.0*M_PI)*u1/(sqrt(r)*2.0*(1.0-nu)); //flip signs as opening displacement is opposite normals
            sifs[nd_a][1] += -mu*sqrt(2.0*M_PI)*u2/(sqrt(r)*2.0*(1.0-nu));
            sifs[nd_a][2] += -mu*sqrt(M_PI)*u3/sqrt(2.0*r);
			//printf("\n .. intermediate K1 at node %d is %.1le", nd_a, -mu*sqrt(2*M_PI)*u1/(sqrt(r)*2*(1-nu)) );
			++count[nd_a];
			
            u1=cod.dot(faceNormals[nd_b]);
			u2=cod.dot(tangents[nd_b].cross(faceNormals[nd_b]));
            u3=cod.dot(tangents[nd_b]);
			if(sifs.count(nd_b)==0) sifs[nd_b].setZero();
            sifs[nd_b][0] += -mu*sqrt(2.0*M_PI)*u1/(sqrt(r)*2.0*(1.0-nu)); //flip signs as opening displacement is opposite normals
            sifs[nd_b][1] += -mu*sqrt(2.0*M_PI)*u2/(sqrt(r)*2.0*(1.0-nu));
            sifs[nd_b][2] += -mu*sqrt(M_PI)*u3/sqrt(2.0*r);
			//printf("\n .. intermediate K1 at node %d is %.1le", nd_b, -mu*sqrt(2*M_PI)*u1/(sqrt(r)*2*(1-nu)) );
			++count[nd_b];
		}
		for(vect3d_map::iterator it=sifs.begin(); it!=sifs.end(); ++it)
			it->second /= count[it->first]; // divide SIFs by number of contributions

		return 0;
	}
	
	int PostProcessor::estimateNodeSIFs(
		const vector_type& u, const vector_type& q, double E, double nu,
		vect3d_map& sifs, vect3d_map& faceNormals, vect3d_map& tangents, state_map& nodeStates
	){
		sifs.clear();
		getCrackTipAxis(faceNormals,tangents);
		// idea: evaluate regular stress field (based on surface displacements and tractions, but ignoring fractures)
		// then project onto the local surface normal
		// and use the crack area to estimate SIFs based on this traction vector

//		for(id_set::iterator it=cracks.begin(); it!=cracks.end(); ++it){
//			printf("\n%% crack %u has surface area %.3lg", *it, crackAreas[*it]);
//		}
		int ret=-1;
		unsigned int i; id_set::iterator it;
		value_map perNodeAreas;
		for(elem_map::iterator ci=crackTips.begin(); ci!=crackTips.end(); ++ci){
			if( perNodeAreas.count(ci->second[0])==0 )
				perNodeAreas[ci->second[0]] = crackAreas[regions[parents[ci->first][0]]];
			if( perNodeAreas.count(ci->second[1])==0 )
				perNodeAreas[ci->second[1]] = crackAreas[regions[parents[ci->first][0]]];
		}
		
		id_set ctNodes = nodeSet(crackTips);
		id_map ctNodeOrder;
		vect3d_map tensions, shears, t_tmp, s_tmp, closestPoint;
		if(!ctNodes.empty()){ //ToDo: check if there are any active markers in the neighborhood of each node and skip nodes where there are none (try not to mess up the result in the smoothing below though!)
			std::vector<Eigen::Vector3d> evalPoints;
			evalPoints.assign(ctNodes.size(),Eigen::Vector3d::Zero());
			for(it=ctNodes.begin(),i=0;it!=ctNodes.end();++it) if( nodeStates[*it]==ACTIVE ) {
				evalPoints[i][0]=nodes[*it][0];
				evalPoints[i][1]=nodes[*it][1];
				evalPoints[i][2]=nodes[*it][2];
				ctNodeOrder[*it]=i;
				tensions[i]=Eigen::Vector3d::Zero();
				shears[i]=Eigen::Vector3d::Zero();
				++i;
			}
//			printf("\n%% SIF estimation: have %d eval points", ctNodes.size());
			//printf("\n u=[\n");
			//for(int i=0;i<u.size();++i) printf("%12.4le%c",u[i],(i%3==2)?'\n':' ');
			//printf("];\n q=[\n");
			//for(int i=0;i<q.size();++i) printf("%12.4le%c",q[i],(i%3==2)?'\n':' ');
			//printf("];\n");
			ret=computeInteriorStresses(t_tmp, s_tmp, closestPoint, evalPoints, u,q,E,nu);
//			printf("\n%% ret=%d\n tensions=[\n"); printMap(tensions); printf(" ]; ");
			
			// do some smoothing of the stresses along the crack-front here ...
			// (interior evals can be noisy close to boundaries
			// an we discontinously snap to a surface eval when too close)
			int passes=4; // could make this a parameter, if we need to ...
			if( passes==0){
				tensions=t_tmp; shears=s_tmp;
			}else for(int pass=1;pass<=passes;++pass){
				for(elem_map::iterator ci=crackTips.begin(); ci!=crackTips.end(); ++ci){
					unsigned int i1=ctNodeOrder[ci->second[0]];
					unsigned int i2=ctNodeOrder[ci->second[1]];
					if( nodeStates[ci->second[0]]==ACTIVE && nodeStates[ci->second[1]]==ACTIVE ){ 
						tensions[i1]+=1.0/6.0*t_tmp[i1];
						tensions[i2]+=1.0/6.0*t_tmp[i2];
						tensions[i1]+=1.0/3.0*t_tmp[i2];
						tensions[i2]+=1.0/3.0*t_tmp[i1];
						shears[i1]+=  1.0/6.0*s_tmp[i1];
						shears[i2]+=  1.0/6.0*s_tmp[i2];
						shears[i1]+=  1.0/3.0*s_tmp[i2];
						shears[i2]+=  1.0/3.0*s_tmp[i1];
					} else
					if( nodeStates[ci->second[0]]==ACTIVE && nodeStates[ci->second[1]]==INACTIVE ){ 
						tensions[i1]+=0.5*t_tmp[i1];
						shears[i1]+=  0.5*s_tmp[i1];
					} else
					if( nodeStates[ci->second[0]]==INACTIVE && nodeStates[ci->second[1]]==ACTIVE ){ 
						tensions[i2]+=0.5*t_tmp[i2];
						shears[i2]+=  0.5*s_tmp[i2];
					}
				}
				if(pass<passes){
					t_tmp=tensions; s_tmp=shears;
					for(it=ctNodes.begin(),i=0;it!=ctNodes.end();++it) if( nodeStates[*it]==ACTIVE ) {
						tensions[i]=Eigen::Vector3d::Zero();
						shears[i]  =Eigen::Vector3d::Zero();
						++i;
					}
				}
			}
		}
		
		if(ret==0) for(it=ctNodes.begin(),i=0;it!=ctNodes.end();++it) if( nodeStates[*it]==ACTIVE ) {
			// get the largest well aligned to faceNormal eigenvector of the local stress
			// compute the ratio of K2/K1, then set K1=eigenvalue, K2=ratio*K1, K3=0
			Eigen::Vector3d n2=tangents[*it].cross(faceNormals[*it]); n2.normalize(); // normalize should be obsolete
			Eigen::Matrix3d stress = Eigen::Matrix3d::Zero();
			stress(0,0)=tensions[i][0]; stress(1,1)=tensions[i][1]; stress(2,2)=tensions[i][2];
			stress(0,1)=shears[i][0]; stress(0,2)=shears[i][1]; stress(1,2)=shears[i][2];
			stress(1,0)=shears[i][0]; stress(2,0)=shears[i][1]; stress(2,1)=shears[i][2];
			Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(stress); // The eigenvalues are sorted in increasing order. (http://eigen.tuxfamily.org/dox/classEigen_1_1SelfAdjointEigenSolver.html)

			unsigned int k=2; // by default choose the largest eigenvalue
			//if( eig.eigenvalues()[0] < -eig.eigenvalues()[2]*compressive) k=0; // somehow this does NOT work
			if( std::abs(eig.eigenvectors().col(k).dot(faceNormals[*it])) < (1.0/3.0)){ // largest ev is misaligned too much (1.0/3.0)
				if( (( eig.eigenvalues()[2-k] < 0.0 && eig.eigenvalues()[2-k] < -eig.eigenvalues()[1]*compressive ) || (eig.eigenvalues()[2-k] > 0.0 && eig.eigenvalues()[2-k] > eig.eigenvalues()[1] )
					&& std::abs(eig.eigenvectors().col(2-k).dot(faceNormals[*it])) > (1.0/3.0)) 
					|| std::abs(eig.eigenvectors().col(1).dot(faceNormals[*it])) < (1.0/3.0) ){
					// have a well aligned compressive ev OR second largest is also misaligned
					k=2-k; // choose compressive direction (most negative ev) OR smallest positive (but only reasonably aligned one)
				}else k=1; // choose second largest ev
			}

			Eigen::Vector3d v = eig.eigenvectors().col(k);
			double ev = eig.eigenvalues()[k];
			if( std::abs( ev ) < (stress*faceNormals[*it]).norm() ) ev = ( ev<0.0?-1.0:1.0 ) * (stress*faceNormals[*it]).norm();
			v -= v.dot(tangents[*it])*tangents[*it]; v.normalize();
			if( v.dot(faceNormals[*it]) <0.0 ) v=-v;
			double p  = v.dot(faceNormals[*it]) / v.norm();
			double th=0.0;
			if( std::abs(p) < 0.9999 ) th = ((v.dot(n2)<0.0)?1.0:-1.0)*acos( p ); // ignore angles less than about 1 deg
			if( std::abs(th) < 1.23095941734077 ){ // 2*atan(sqrt(8)/4)) or acos(1/3.)
				// TRY NOT TO USE THIS FOR THE RESULTS - WE WANT TO AVOID HAVING THIS IN THE PAPER
				// force cracks towards boundary if the closest point on the surface is within 35deg of the forward direction
				//Eigen::Vector3d pt,cp; copyNode(*it, pt);
				//cp = closestPoint[i] - pt;
				//if( cp.norm() < (2.0*crackMeshSize) ){
				//	findClosestPointInPlane(cp,*it,tangents[*it]);
				//	cp-=pt;
				//}
				//if( cp.dot(n2)/cp.norm() >= 0.7 && cp.norm() < (2.0*crackMeshSize) ){ //closer than about 35deg to the forward direction and within 2 propagation steps
				//	cp-=cp.dot(tangents[*it])*tangents[*it];
				//	double cp_th=acos(cp.dot(n2)/cp.norm())* ((cp.dot(faceNormals[*it]))>0.0?1.0:-1.0);
				//	double w = (2.0*crackMeshSize - cp.norm())/(2.0*crackMeshSize);
				//	th = (1.0-w)*th+w*cp_th;

				//	//// only for testing
				//	//double diffTh = th - acos(cp.dot(n2)/cp.norm())* ((cp.dot(faceNormals[*it]))>0.0?1.0:-1.0);
				//	//if( diffTh > 0.45 ) diffTh = 0.45;
				//	//if( diffTh <-0.45 ) diffTh =-0.45; // change by at most about 25deg
				//	//th -= diffTh; //printf("!*!*!");
				//	//th=acos(cp.dot(n2)/cp.norm())* ((cp.dot(faceNormals[*it]))>0.0?1.0:-1.0);
				//	//printf("\n%% *%d th= %.3lf  \t dist= %.3lf",*it, th, cp.norm());
				//}

				double r  = tan(0.5*th) / (2.0*tan(0.5*th)*tan(0.5*th) -1.0);
				double r_ = sqrt(1.0/(r*r+1.0)); // scale s.t. K1^2 + K2^2 = eig.val.^2
				sifs[*it][0] =   r_*ev; 
				sifs[*it][1] = r*r_*ev;
				sifs[*it][2] = 0.0;
				
				v=eig.eigenvectors().col(k);
				v-=v.dot(n2)*n2; v.normalize();
				if( v.dot(faceNormals[*it]) <0.0 ) v=-v;
				r_=v.dot(tangents[*it]); // r_ = cos(alpha) with alpha the angle between n3 and v in the n1 x n3 plane
				sifs[*it][2] = (r_<0.0?1.0:-1.0)*sqrt(std::abs(r_)*(1.0-nu))*sifs[*it][0];
				sifs[*it][0]*= sqrt(1.0-std::abs(r_));
				// K3 = sqrt(r * (1-nu))*K, K1 = sqrt(1-r)*K such that
				// K13 combined ==> K^2 = K1^2 + K3^2/(1-nu)
			}else{
//				printf("\n%% WARNING: MISALIGNED EIGENVECTORS OF STRESS FOUND !!!\n");
//				printf("\n stress = [\n"); std::cout << stress;
//				printf("];\n eigval = [\n"); std::cout << eig.eigenvalues();
//				printf("];\n eigvec = [\n"); std::cout << eig.eigenvectors();
//				printf("];\n%% chose %d, angle %.4lf, ratio %.4lf", k+1, th*180.0/M_PI, r);
//				printf("\n n1 = [\n"); std::cout << faceNormals[*it]; printf("];");
//				printf("\n n3 = [\n"); std::cout << tangents[*it]; printf("];");
//				printf("\n%% projections: %.3lf, %.3lf, %.3lf",
//					std::abs(evecs[0].dot(faceNormals[*it])),
//					std::abs(evecs[1].dot(faceNormals[*it])),
//					std::abs(evecs[2].dot(faceNormals[*it])) );

				sifs[*it][0] = 0.0;
				sifs[*it][1] = ((th<0.0)?1.0:-1.0)*ev;
				sifs[*it][2] = 0.0;
			}
			sifs[*it] *= sqrt(sqrt(M_PI*perNodeAreas[*it])) * 1.1215/** 2.0/M_PI */;
			// constant 1.1215 from Gross&Selig Table 4.1 (5) works best for inhom Neumann problems (2/pi works well for inhom. Dirichlet problems)

			//// basic estimator (crack path oscillates!)
			//Eigen::Vector3d localTraction = stress * faceNormals[*it];
			//sifs[*it][0] =  localTraction.dot(faceNormals[*it]);
			//sifs[*it][1] =  localTraction.dot(n2);
			//sifs[*it][2] =  localTraction.dot(tangents[*it]);
			//sifs[*it]   *=  sqrt(sqrt(M_PI*perNodeAreas[*it])) * 1.1215;
			++i;
		}
		return ret;
	}
    
	// NEW FOR FractureRB
	int PostProcessor::estimateCOD(
		const elem_map& newElems, vect3d_map& nodeSIFs, const id_set& ctNodes,
		double E, double nu, vect3d_map& estCOD
	){
		// We assume that the crack has just propagated, creating the elements in newElems
		// BUT the SIFs stored in nodeSIFs are still the values used in the propagation step
		// as opposed to the ones we'll eventually prepare for the next propagation step
		// we assume the global mesh data to be up to date including newElems and crack-tip elements
		// here we aim to update regular surface data (could be displacements, tractions or BEM RHS vectors)
		// to account for the effect of the newly opened crack onto the surface.
		// ToDo: find a better name for this function!

		double mu = E/(2*(1+nu));
		unsigned int haveSIFs[3];
		unsigned int edg_x,edg_y,nd_z; bool have_x;
		unsigned int nd_z_is_old;
		//vect3d_map estCOD; // estimate crack-opening displacements
		for(elem_map::const_iterator it=newElems.begin(); it!=newElems.end(); ++it){
            // each element in newElems contains either 2 new nodes and 1 old node, or 1 new node and 2 old nodes
			// re-order node-IDs such that the correlation distance in this element is the shortest distance from edge (x,y) to node z
			for(unsigned int k=0; k<3; ++k) haveSIFs[k]=nodeSIFs.count(it->second[k]);
			if(     (haveSIFs[0]+haveSIFs[1]+haveSIFs[2])==1) nd_z_is_old=1; // only 1 old node
			else if((haveSIFs[0]+haveSIFs[1]+haveSIFs[2])==2) nd_z_is_old=0; // 2 old nodes
			else return -1; // this should never happen
			have_x=false;
			for(unsigned int k=0; k<3; ++k){
				if( haveSIFs[k]==nd_z_is_old ) nd_z=it->second[k]; // this is the opposite node in this element
				else{ // this is one of the nodes in edge (x,y)
					if(!have_x){ edg_x=it->second[k]; have_x=true; }
					else edg_y=it->second[k];
				}
			}
			//printf("\n%% haveSIFs in %d is (%d %d %d) at (%d, %d, %d)",it->first, haveSIFs[0],haveSIFs[1],haveSIFs[2], it->second[0],it->second[1],it->second[2]);
			//printf("\n%% reordered element %d to (%d--%d  :: %d)", it->first, edg_x,edg_y, nd_z);
			// now we can compute the correlation distance in this element
			Eigen::Vector3d x,y,z,tmp,edg;
			copyNode(edg_x, x);
            copyNode(edg_y, y);
            copyNode(nd_z , z);
			edg=y-x; edg.normalize();
			tmp = 0.5*(x+y) - z; // vector from edge-midpoint to opposite node
			tmp-= edg*tmp.dot(edg); // remove parallel component
			double r = tmp.norm(); // correlation distance

			if( nd_z_is_old==1 && estCOD.count(nd_z)==0 && ctNodes.count(nd_z)==0 ){ // only nd_z is old, and it is not yet done --> apply inverse displacement correlation
				estCOD[nd_z].setZero();
				estCOD[nd_z][0] = -sqrt(r)*2.0*(1.0-nu)/(mu*sqrt(2.0*M_PI))* nodeSIFs[nd_z][0]; //flip signs as opening displacement is opposite normals
				estCOD[nd_z][1] = -sqrt(r)*2.0*(1.0-nu)/(mu*sqrt(2.0*M_PI))* nodeSIFs[nd_z][1];
				estCOD[nd_z][2] = -sqrt(2.0*r)/(mu*sqrt(M_PI))             * nodeSIFs[nd_z][2];
				//printf("\n%% COD by iDCT ndz in %d: (%.1le, %.1le, %.1le)", it->first, estCOD[nd_z][0], estCOD[nd_z][1], estCOD[nd_z][2]);
			}else if( nd_z_is_old==0 ){ // have 2 old nodes, check both edg_x and edg_y
				if( estCOD.count(edg_x)==0 && ctNodes.count(edg_x)==0 ){ // apply inverse displacement correlation for edg_x
					estCOD[edg_x].setZero();
					estCOD[edg_x][0] = -sqrt(r)*2.0*(1.0-nu)/(mu*sqrt(2.0*M_PI))* nodeSIFs[edg_x][0]; //flip signs as opening displacement is opposite normals
					estCOD[edg_x][1] = -sqrt(r)*2.0*(1.0-nu)/(mu*sqrt(2.0*M_PI))* nodeSIFs[edg_x][1];
					estCOD[edg_x][2] = -sqrt(2.0*r)/(mu*sqrt(M_PI))             * nodeSIFs[edg_x][2];
					//printf("\n%% COD by iDCT edx in %d: (%.1le, %.1le, %.1le)", it->first, estCOD[edg_x][0], estCOD[edg_x][1], estCOD[edg_x][2]);
				}
				if( estCOD.count(edg_y)==0 && ctNodes.count(edg_y)==0 ){ // apply inverse displacement correlation for edg_y
					estCOD[edg_y].setZero();
					estCOD[edg_y][0] = -sqrt(r)*2.0*(1.0-nu)/(mu*sqrt(2.0*M_PI))* nodeSIFs[edg_y][0]; //flip signs as opening displacement is opposite normals
					estCOD[edg_y][1] = -sqrt(r)*2.0*(1.0-nu)/(mu*sqrt(2.0*M_PI))* nodeSIFs[edg_y][1];
					estCOD[edg_y][2] = -sqrt(2.0*r)/(mu*sqrt(M_PI))             * nodeSIFs[edg_y][2];
					//printf("\n%% COD by iDCT edy in %d: (%.1le, %.1le, %.1le)", it->first, estCOD[edg_y][0], estCOD[edg_y][1], estCOD[edg_y][2]);
				}
			}
		}
		// convert from local crack-tip coordinates to cartesian coords
		for(vect3d_map::iterator it=estCOD.begin(); it!=estCOD.end(); ++it){
			it->second = lastCrackTipFaceNormals[it->first]                                        * (it->second[0])
						+lastCrackTipTangents[it->first].cross(lastCrackTipFaceNormals[it->first]) * (it->second[1])
						+lastCrackTipTangents[it->first]                                           * (it->second[2]);
		}
		// now we have estimated crack opening displacements on all nodes that had SIFs and are not part of the crack-tip anymore
		// we want to update the RHS of the BEM solution with these new CODs, so we need to find the support elements of these nodes
		// integrate between these and the regular surface elements, then update the BEM solution

		return 0;
	}

    int PostProcessor::vtkWrite(std::string fileName, VTK_CELL_TYPE cellType){
        int nNodes = nodes.size();
        int nCoordsPerNode = nodes.begin()->second.size();
        int nElems = elems.size();
        int nNodesPerElem = elems.begin()->second.size();
		ofstream out(fileName.c_str());
        
        if(!out.is_open()) return -1;
        out.precision(12);
        //out.setf(std::ios::scientific);
        //write header
        out << "# vtk DataFile Version 2.0" << endl;
        out << "VTK exported mesh" << endl;
        out << "ASCII" << endl;
        out << "DATASET UNSTRUCTURED_GRID" << endl;
        //node coordinates
        out << "POINTS " << nNodes << " double" << endl;
//        for(int i=0;i<nNodes; ++i){
        for(node_map::iterator i=nodes.begin(); i!=nodes.end(); ++i){
//            out << nodes[3*i] << " " << nodes[3*i+1] << " " << nodes[3*i+2] << endl;
            for(int j=0; j<nCoordsPerNode; ++j) out << i->second[j] << " ";
            for(int j=nCoordsPerNode; j<3; ++j) out << "0.0 "; // fill with zeros if less than 3 coords per node given
            out << endl;
        }
        //cells
        out << "CELLS " << nElems << " " <<  (1+nNodesPerElem)*nElems << endl;
//        for(int i=0;i<nElems; ++i){
        for(elem_map::iterator i=elems.begin(); i!=elems.end(); ++i){
            out << nNodesPerElem << " ";
            // VTK Unstructured Grid identifies nodes by 0-based index into POINTS list
            for(int j=0; j<nNodesPerElem; ++j) out << i->second[j]-NODE_BASE_INDEX << " ";
            out << endl;
        }
        //cell types
        out << "CELL_TYPES " << nElems << endl;
        for(int i=0;i<nElems; ++i) out << cellType << endl;
        // mesh information is done now, next up point and cell data ...
        
        out << "POINT_DATA " << nNodes << endl;
        for(int i=0;i<vtkDataNames.size(); ++i){
            if(vtkDataIsCellData[i]==false) // this is point data
                vtkWriteData(out, i, nNodes);
        }

        out << "CELL_DATA " << nElems << endl;
        for(int i=0;i<vtkDataNames.size(); ++i){
            if(vtkDataIsCellData[i]==true) // this is cell data
                vtkWriteData(out, i, nElems);
        }

        return 0;
    }
    
    void PostProcessor::vtkWriteData(ofstream& out, int i, int n){
        if(vtkDataDimension[i]==1){ //scalars
            out << "SCALARS " << vtkDataNames[i] << " double" << endl;
            out << "LOOKUP_TABLE default" << endl; //could specify color lookup
            for(int j=0;j<n; ++j) out << (*vtkData[i])[j] << endl;
        }else if(vtkDataDimension[i]==3){ //vectors
            out << "VECTORS " << vtkDataNames[i] << " double" << endl;
            for(int j=0;j<n; ++j)
                out << (*vtkData[i])[3*j  ] << " "
                    << (*vtkData[i])[3*j+1] << " "
                    << (*vtkData[i])[3*j+2] << endl;
        }
    }
    
    int PostProcessor::vtkAddData(
            std::string name, int dimension, bool isCellData, vector_type& data
    ){
        if( //check that number of data entries matches expected value
            (isCellData==true  && data.rows()==dimension*elems.size()) ||
            (isCellData==false && data.rows()==dimension*nodes.size())
        ){
            vtkDataNames.push_back(name);
            vtkDataDimension.push_back(dimension);
            vtkDataIsCellData.push_back(isCellData);
            vtkData.push_back(&data);
            return 0;
        }else{
//            fprintf(stderr, "can't add %s data because data size is %d but expected %d (dim=%d)\n",
//                isCellData?"cell":"point",data.rows(),
//                isCellData?(dimension*elems.size()):(dimension*nodes.size()), dimension);
            return -1;
        }
    }

	int PostProcessor::computeSurfaceStresses(
        node_map& maxPrincipalStress, const vector_type& u, const vector_type& u_c,
        const vector_type& q, double E, double nu, bool addToVTK, int ignoreElems
    ){
		// compute principal and cartesian stresses based on the linear FEM formulation found in
		// Bargteil et al. 2007 "A Finite Element Method for Animating Large Viscoplastic Flow"
		// by extending every triangle to a tet with no out-of-plane deformation
		// ...
		unsigned int n = elems.size(), i=0;
		s_xx.resize(n);
		s_yy.resize(n);
		s_zz.resize(n);
		s_xy.resize(n);
		s_xz.resize(n);
		s_yz.resize(n);
        maxPrincipalStress.clear();

		// U*P*V' SVD of stress S
		Eigen::Matrix3d U,S,Vt,P;

		for(elem_map::iterator it=elems.begin(); it!=elems.end(); ++it, ++i) // loop over elements
		if(((ignoreElems < 0) || ( it->first < ignoreElems ) ) &&  q.size() > (3*(it->first - ELEM_BASE_INDEX)+2) ){
			Eigen::Vector3d // a,b,c are node coordinates in material space; ua,ub,uc are nodal displacements, so world coordinates are a+ua etc.
				a (nodes[it->second[0]][0], nodes[it->second[0]][1], nodes[it->second[0]][2]),
				b (nodes[it->second[1]][0], nodes[it->second[1]][1], nodes[it->second[1]][2]),
				c (nodes[it->second[2]][0], nodes[it->second[2]][1], nodes[it->second[2]][2]),
				ua(u[3*(it->second[0]-NODE_BASE_INDEX)], u[3*(it->second[0]-NODE_BASE_INDEX)+1], u[3*(it->second[0]-NODE_BASE_INDEX)+2]),
				ub(u[3*(it->second[1]-NODE_BASE_INDEX)], u[3*(it->second[1]-NODE_BASE_INDEX)+1], u[3*(it->second[1]-NODE_BASE_INDEX)+2]),
				uc(u[3*(it->second[2]-NODE_BASE_INDEX)], u[3*(it->second[2]-NODE_BASE_INDEX)+1], u[3*(it->second[2]-NODE_BASE_INDEX)+2]),
				uan,ubn,ucn; // for negative COD of fracture elements
            
			if(cracks.count(regions[it->first])!=0){ // this is a crack element
				if(u_c.rows() > 0){
					Eigen::Vector3d // uba,ubb,ubc are nodal crack base displacements
						uba(u_c[3*(it->second[0]-NODE_BASE_INDEX)], u_c[3*(it->second[0]-NODE_BASE_INDEX)+1], u_c[3*(it->second[0]-NODE_BASE_INDEX)+2]),
						ubb(u_c[3*(it->second[1]-NODE_BASE_INDEX)], u_c[3*(it->second[1]-NODE_BASE_INDEX)+1], u_c[3*(it->second[1]-NODE_BASE_INDEX)+2]),
						ubc(u_c[3*(it->second[2]-NODE_BASE_INDEX)], u_c[3*(it->second[2]-NODE_BASE_INDEX)+1], u_c[3*(it->second[2]-NODE_BASE_INDEX)+2]);
					ua*=0.5; ub*=0.5; uc*=0.5; // use half of crack opening displacement (because COD is distributed over 2 surfaces)
					uan=uba-ua; ubn=ubb-ub; ucn=ubc-uc; // negative COD + base
					ua+=uba; ub+=ubb; uc+=ubc;          // positive COD + base displacement
				}else{ // we're not them crack elements
					ua.setZero();ub.setZero();uc.setZero();
				}
            }
			Eigen::Vector3d qe(
				q[3*(it->first - ELEM_BASE_INDEX)], q[3*(it->first - ELEM_BASE_INDEX)+1], q[3*(it->first - ELEM_BASE_INDEX)+2]
			);
			if( computeElementPrincipalStresses(U,S,Vt,P, a,b,c, ua,ub,uc, qe, E,nu) !=0) return -1;

            maxPrincipalStress[it->first].assign(9,0.0);
			// store max. principal stress and its plane-normal vector for each element (also min. for compressive fracture?)
            // in the SVD, V is the pre-rotation and U is the post-rotation
            // so the direction of max. principal stress is the first column of V
            // which becomes the first row of V.transpose()
            maxPrincipalStress[it->first][0]=P(0,0);  // first entry is the stress value
            maxPrincipalStress[it->first][1]=Vt(0,0); // next 3 entries
            maxPrincipalStress[it->first][2]=Vt(0,1); // are the plane-
            maxPrincipalStress[it->first][3]=Vt(0,2); // normal vector
            // same for min. principal stress used for compressive fractures
			maxPrincipalStress[it->first][4]=P(2,2);
			maxPrincipalStress[it->first][5]=Vt(2,0); // next 3 entries
            maxPrincipalStress[it->first][6]=Vt(2,1); // are the plane-
            maxPrincipalStress[it->first][7]=Vt(2,2); // normal vector
            //s_xx[i]=P(0,0); s_yy[i]=P(1,1); s_zz[i]=P(2,2); // for testing

			// flag fracture elements based on sign of COD
			// [8] --> 1: both positive; 2: max. negative, min. positive; 3: max. pos., min. neg.; 4: both neg.
			if((u_c.rows() > 0) && cracks.count(regions[it->first])!=0){
				maxPrincipalStress[it->first][8]=1.0;
				bool tmp=false;
				// repeat stress computation with crack-base displacement MINUS half crack-opening displacement and see which sign of the COD produces more principal stress
				if( computeElementPrincipalStresses(U,S,Vt,P, a,b,c, uan,ubn,ucn, -qe, E,nu) !=0) return -2;
				if( std::abs(maxPrincipalStress[it->first][0]) < std::abs(P(0,0)) ){
					maxPrincipalStress[it->first][0]=P(0,0);  // first entry is the stress value
					maxPrincipalStress[it->first][1]=Vt(0,0); // next 3 entries
					maxPrincipalStress[it->first][2]=Vt(0,1); // are the plane-
					maxPrincipalStress[it->first][3]=Vt(0,2); // normal vector
					maxPrincipalStress[it->first][8]=2.0; tmp=true;
					//printf("\n%% *** have neg. ten.");
				}
				if( std::abs(maxPrincipalStress[it->first][4]) < std::abs(P(2,2)) ){
					maxPrincipalStress[it->first][4]=P(2,2);
					maxPrincipalStress[it->first][5]=Vt(2,0); // next 3 entries
					maxPrincipalStress[it->first][6]=Vt(2,1); // are the plane-
					maxPrincipalStress[it->first][7]=Vt(2,2); // normal vector
					maxPrincipalStress[it->first][8]=tmp?(4.0):(3.0);
					//printf("\n%% *** have neg. comp.");
				}
			}

			if(addToVTK){
				P = U*P*Vt; // P is now a matrix of cartesian stresses
				//printf("\n%% -- cartesian: %.3le, %.3le, %.3le", P(0,0), P(1,1), P(2,2));
				//printf("\n%% --  w/ shear: %.3le, %.3le, %.3le", P(0,1), P(0,2), P(1,2));
                s_xx[i]=P(0,0); s_yy[i]=P(1,1); s_zz[i]=P(2,2);
                s_xy[i]=0.5*( P(0,1)+P(1,0) ); s_xz[i]=0.5*( P(0,2)+P(2,0) ); s_yz[i]=0.5*( P(2,1)+P(1,2) );
            }
		}
        if(addToVTK){
            vtkAddData("stress_xx",1,true,s_xx);
            vtkAddData("stress_yy",1,true,s_yy);
            vtkAddData("stress_zz",1,true,s_zz);
            vtkAddData("stress_xy",1,true,s_xy);
            vtkAddData("stress_xz",1,true,s_xz);
            vtkAddData("stress_yz",1,true,s_yz);
        }
		return 0;
	}
	/**/

	int PostProcessor::getCrackTipAxis(vect3d_map& faceNormals, vect3d_map& tangents){
		faceNormals.clear();
		tangents.clear();
		for(elem_map::iterator it=crackTips.begin(); it!=crackTips.end(); ++it){
            unsigned int nd_a, nd_b, nd_c; // nodes of the current triangle
            nd_a = it->second[0]; nd_b = it->second[1];
            nd_c = findInteriorNode(nd_a, nd_b, elems[parents[it->first][0]]);
            Eigen::Vector3d a,b,c, n1,n2,n3;
            copyNode(nd_a, a);
            copyNode(nd_b, b);
            copyNode(nd_c, c);
            getLocalCoordFrame(a,b,c, n1,n2,n3);

            if(tangents.count(nd_a)==0) tangents[nd_a]=Eigen::Vector3d(0.0,0.0,0.0);
            if(tangents.count(nd_b)==0) tangents[nd_b]=Eigen::Vector3d(0.0,0.0,0.0);
            if(faceNormals.count(nd_a)==0) faceNormals[nd_a]=Eigen::Vector3d(0.0,0.0,0.0);
            if(faceNormals.count(nd_b)==0) faceNormals[nd_b]=Eigen::Vector3d(0.0,0.0,0.0);
            faceNormals[nd_a]+=n1;
            faceNormals[nd_b]+=n1;
            tangents[nd_a]+=n3;
            tangents[nd_b]+=n3;
		}
		for(vect3d_map::iterator it=faceNormals.begin(); it!=faceNormals.end(); ++it){ // normalize face normals
			it->second.normalize();
			//printf("\n .. normal  at node %d:\t(%.3lf, %.3lf, %.3lf)", it->first, it->second[0], it->second[1], it->second[2]);
		}
		for(vect3d_map::iterator it=tangents.begin(); it!=tangents.end(); ++it){ // normalize tangents
			it->second -= it->second.dot(faceNormals[it->first])*faceNormals[it->first];
			it->second.normalize();
			//printf("\n .. tangent at node %d:\t(%.3lf, %.3lf, %.3lf)", it->first, it->second[0], it->second[1], it->second[2]);
		}

		// new for method INIT_BEM: keep a local copy (will be used by COD estimator in the getSurfaceUpdate)
		lastCrackTipFaceNormals = faceNormals;
		lastCrackTipTangents    = tangents;
		return 0;
	}

	int PostProcessor::computeElementDefGrad(
		Eigen::Matrix3d& F,
		const Eigen::Vector3d&  a, const Eigen::Vector3d&  b, const Eigen::Vector3d&  c,
		const Eigen::Vector3d& ua, const Eigen::Vector3d& ub, const Eigen::Vector3d& uc
	){
		printf("\n!!! DON'T USE PostProcessor::computeElementDefGrad -- old version !!!\n");
		Eigen::Matrix3d B,X;
		bool flag=true;
		X.col(0) = a - c; // using world space matrix as temporary storage before inversion
		X.col(1) = b - c;
		X.col(2) = X.col(0).cross(X.col(1)).normalized();
		X.computeInverseWithCheck(B,flag);
		if(!flag) return -1; // matrix not invertible ~> possibly degenerate element
		// now build the actual world space matrix
		X.col(0) = a+ua -c-uc;
		X.col(1) = b+ub -c-uc;
		X.col(2) = X.col(0).cross(X.col(1)).normalized();
		F=X*B;
		return 0;
	}
    int PostProcessor::computeElementStresses(
		Eigen::Matrix3d& S,
		const Eigen::Vector3d&  a, const Eigen::Vector3d&  b, const Eigen::Vector3d&  c,
		const Eigen::Vector3d& ua, const Eigen::Vector3d& ub, const Eigen::Vector3d& uc,
		const Eigen::Vector3d&  q, double E, double nu
	){
		/** old version -- not accurate
		Eigen::Matrix3d F, I=Eigen::Matrix3d::Identity();
		if( computeElementDefGrad(F, a,b,c,ua,ub,uc) !=0) return -1;
		F.noalias() = 0.5*(F+F.transpose());
		S = 2*mu*(F-I) + lambda*(F-I).trace()*I;
		/*/
		// new version:
		// flatten the undeformed and deformed triangle into a plane
		// then compute in-plane stress using hat-function gradients and Hooke's law
		// finally, rotate the stress back to 3d material space
		// and add the out-of-plane traction
		Eigen::Vector3d n1,n2,n3;
		getLocalCoordFrame(a,b,c, n3,n2,n1); // we use n1 as tangent to edge ab and n3 as face-normal here
		n2*=-1.0; // usually we use n2 as "forward" direction in-plane OUTWARD normal to the edge ab, here we want an inward version
		double e1= (b-a).norm(); // length of edge ab
		double e2= (c-a).dot(b-a)/e1; // projection of edge ac onto edge ab
		double e3= ((c-a)-(b-a)/e1*e2).norm(); // normal distance of c from edge ab
		// similarly for the deformed version
		double f1= ((b+ub)-(a+ua)).norm();
		double f2= ((c+uc)-(a+ua)).dot((b+ub)-(a+ua))/f1;
		double f3= ((c+uc)-(a+ua)-((b+ub)-(a+ua))/f1*f2).norm();
		// the stress computation goes like this:
		// E1=E/(1-nu*nu);
		// G =E/2/(1+nu);
		// D =[ E1    nu*E1 0
		//	    nu*E1    E1 0
		//	     0     0    G ];
		// B =[ e3     0    -e3   0   0   0
		//	     0    e1-e2   0  e2   0 -e1
		//	    e1-e2 e3     e2 -e3 -e1   0 ]/norm(cross(c-a,b-a));
		// u =[0;0; e1-f1;0; e2-f2;e3-f3]; % in-plane displacement of nodes
		// Sv=D*B*u;
		// S =[ Sv(1) Sv(3); Sv(3) Sv(2) ];
		// which simplifies to:
		double A2  = (c-a).cross(b-a).norm(); // twice the element's area
		double sv1 = (E*(e1*e3 - e3*f1) + E*nu*(e1*e3 - e1*f3))/(A2*(nu*nu - 1.0));
		double sv2 = (E*(e1*e3 - e1*f3) + E*nu*(e1*e3 - e3*f1))/(A2*(nu*nu - 1.0));
		double sv3 = (E*(e1*f2 - e2*f1))/(2.0*A2*(nu + 1.0));
		// the 2d in-plane stress is then S = [sv1, sv3 ; sv3 , sv2];
		// finally, we transform this stress back to 3d like this:
		// A=[ n1 n2 ]; % rotates the 2d standard basis vectors onto n1 and n2 respectively
		// % then we apply the rotation A to the 2d stress to get the traction vectors
		// % wrt. the in-plane normals (n1,n2) of the 3d triangle
		// qi=A*S;
		// % now build the stress as a dyadic of the
		// % traction vectors and their normals
		// % we use the in-plane tractions from the displacement
		// % and the out-of-plane traction from the BEM result
		// qi=qi+n3*qn'*A; % correct n1-n3 and n2-n3 shearing tractions
		// qi=qi+nu/(1-nu)*qn'*n3*A; % experimental correction for Poisson contraction -- explain!
		// S=[qi qn] *[A n3]'; % build stress matrix from tractions
		// all these matrix-vector multiplications simplify to: (ask MuPAD or Mathematica or something)
		// without Poisson effect correction:
//		S(0,0) = n3[0]*q[0] + n1[0]*(n1[0]*sv1 + n2[0]*sv3 + n3[0]*(n1[0]*q[0] + n1[1]*q[1] + n1[2]*q[2])) + n2[0]*(n1[0]*sv3 + n2[0]*sv2 + n3[0]*(n2[0]*q[0] + n2[1]*q[1] + n2[2]*q[2])) ;
//		S(1,1) = n3[1]*q[1] + n1[1]*(n1[1]*sv1 + n2[1]*sv3 + n3[1]*(n1[0]*q[0] + n1[1]*q[1] + n1[2]*q[2])) + n2[1]*(n1[1]*sv3 + n2[1]*sv2 + n3[1]*(n2[0]*q[0] + n2[1]*q[1] + n2[2]*q[2])) ;
//		S(2,2) = n3[2]*q[2] + n1[2]*(n1[2]*sv1 + n2[2]*sv3 + n3[2]*(n1[0]*q[0] + n1[1]*q[1] + n1[2]*q[2])) + n2[2]*(n1[2]*sv3 + n2[2]*sv2 + n3[2]*(n2[0]*q[0] + n2[1]*q[1] + n2[2]*q[2])) ;
//		S(0,1) = n3[1]*q[0] + n1[1]*(n1[0]*sv1 + n2[0]*sv3 + n3[0]*(n1[0]*q[0] + n1[1]*q[1] + n1[2]*q[2])) + n2[1]*(n1[0]*sv3 + n2[0]*sv2 + n3[0]*(n2[0]*q[0] + n2[1]*q[1] + n2[2]*q[2])) ;
//		S(0,2) = n3[2]*q[0] + n1[2]*(n1[0]*sv1 + n2[0]*sv3 + n3[0]*(n1[0]*q[0] + n1[1]*q[1] + n1[2]*q[2])) + n2[2]*(n1[0]*sv3 + n2[0]*sv2 + n3[0]*(n2[0]*q[0] + n2[1]*q[1] + n2[2]*q[2])) ;
//		S(1,2) = n3[2]*q[1] + n1[2]*(n1[1]*sv1 + n2[1]*sv3 + n3[1]*(n1[0]*q[0] + n1[1]*q[1] + n1[2]*q[2])) + n2[2]*(n1[1]*sv3 + n2[1]*sv2 + n3[1]*(n2[0]*q[0] + n2[1]*q[1] + n2[2]*q[2])) ;
		// with Poisson effect correction
		S(0,0) = n1[0]*(n1[0]*sv1 + n2[0]*sv3 + n3[0]*(n1[0]*q[0] + n1[1]*q[1] + n1[2]*q[2]) - (n1[0]*nu*(n3[0]*q[0] + n3[1]*q[1] + n3[2]*q[2]))/(nu - 1)) + n2[0]*(n1[0]*sv3 + n2[0]*sv2 + n3[0]*(n2[0]*q[0] + n2[1]*q[1] + n2[2]*q[2]) - (n2[0]*nu*(n3[0]*q[0] + n3[1]*q[1] + n3[2]*q[2]))/(nu - 1)) + n3[0]*q[0] ;
		S(1,1) = n1[1]*(n1[1]*sv1 + n2[1]*sv3 + n3[1]*(n1[0]*q[0] + n1[1]*q[1] + n1[2]*q[2]) - (n1[1]*nu*(n3[0]*q[0] + n3[1]*q[1] + n3[2]*q[2]))/(nu - 1)) + n2[1]*(n1[1]*sv3 + n2[1]*sv2 + n3[1]*(n2[0]*q[0] + n2[1]*q[1] + n2[2]*q[2]) - (n2[1]*nu*(n3[0]*q[0] + n3[1]*q[1] + n3[2]*q[2]))/(nu - 1)) + n3[1]*q[1] ;
		S(2,2) = n1[2]*(n1[2]*sv1 + n2[2]*sv3 + n3[2]*(n1[0]*q[0] + n1[1]*q[1] + n1[2]*q[2]) - (n1[2]*nu*(n3[0]*q[0] + n3[1]*q[1] + n3[2]*q[2]))/(nu - 1)) + n2[2]*(n1[2]*sv3 + n2[2]*sv2 + n3[2]*(n2[0]*q[0] + n2[1]*q[1] + n2[2]*q[2]) - (n2[2]*nu*(n3[0]*q[0] + n3[1]*q[1] + n3[2]*q[2]))/(nu - 1)) + n3[2]*q[2] ;
		S(0,1) = n1[1]*(n1[0]*sv1 + n2[0]*sv3 + n3[0]*(n1[0]*q[0] + n1[1]*q[1] + n1[2]*q[2]) - (n1[0]*nu*(n3[0]*q[0] + n3[1]*q[1] + n3[2]*q[2]))/(nu - 1)) + n2[1]*(n1[0]*sv3 + n2[0]*sv2 + n3[0]*(n2[0]*q[0] + n2[1]*q[1] + n2[2]*q[2]) - (n2[0]*nu*(n3[0]*q[0] + n3[1]*q[1] + n3[2]*q[2]))/(nu - 1)) + n3[1]*q[0] ;
		S(0,2) = n1[2]*(n1[0]*sv1 + n2[0]*sv3 + n3[0]*(n1[0]*q[0] + n1[1]*q[1] + n1[2]*q[2]) - (n1[0]*nu*(n3[0]*q[0] + n3[1]*q[1] + n3[2]*q[2]))/(nu - 1)) + n2[2]*(n1[0]*sv3 + n2[0]*sv2 + n3[0]*(n2[0]*q[0] + n2[1]*q[1] + n2[2]*q[2]) - (n2[0]*nu*(n3[0]*q[0] + n3[1]*q[1] + n3[2]*q[2]))/(nu - 1)) + n3[2]*q[0] ;
		S(1,2) = n1[2]*(n1[1]*sv1 + n2[1]*sv3 + n3[1]*(n1[0]*q[0] + n1[1]*q[1] + n1[2]*q[2]) - (n1[1]*nu*(n3[0]*q[0] + n3[1]*q[1] + n3[2]*q[2]))/(nu - 1)) + n2[2]*(n1[1]*sv3 + n2[1]*sv2 + n3[1]*(n2[0]*q[0] + n2[1]*q[1] + n2[2]*q[2]) - (n2[1]*nu*(n3[0]*q[0] + n3[1]*q[1] + n3[2]*q[2]))/(nu - 1)) + n3[2]*q[1] ;

		S(1,0) = S(0,1); S(2,0) = S(0,2); S(2,1) = S(1,2); // symmetric by construction
		/**/
		return 0;
	}
	int PostProcessor::computeElementPrincipalStresses(
		Eigen::Matrix3d& U, Eigen::Matrix3d& S, Eigen::Matrix3d& Vt, Eigen::Matrix3d& P,
		const Eigen::Vector3d&  a, const Eigen::Vector3d&  b, const Eigen::Vector3d&  c,
		const Eigen::Vector3d& ua, const Eigen::Vector3d& ub, const Eigen::Vector3d& uc,
		const Eigen::Vector3d&  q, double E, double nu
	){
		/** // old version
		Eigen::Matrix3d F, I=Eigen::Matrix3d::Identity();
		if( computeElementDefGrad(F, a,b,c,ua,ub,uc) !=0) return -1;
		F.noalias() = 0.5*(F+F.transpose());
		// compute singular value decomposition of F --> U*S*V'
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(F,Eigen::ComputeFullU | Eigen::ComputeFullV);
		S = svd.singularValues().asDiagonal(); // these will be returned in order (highest to lowest)
		U = svd.matrixU();
		Vt= svd.matrixV().transpose();
		P = 2*mu*(S-I) + lambda*(S-I).trace()*I; // P is now a diagonal matrix of principal stresses
		//printf("\n%% diag. deform: %.3le, %.3le, %.3le", S(0,0), S(1,1), S(2,2));
        //printf("\n%% stresses for element %d:",it->first);
        //printf("\n%% -- principal: %.3le, %.3le, %.3le", P(0,0), P(1,1), P(2,2));
		/*/ // new version
		if( computeElementStresses(S, a,b,c, ua,ub,uc, q, E,nu) !=0) return -1;
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(S,Eigen::ComputeFullU | Eigen::ComputeFullV);
		P = svd.singularValues().asDiagonal(); // these will be returned in order (highest to lowest)
		U = svd.matrixU();
		Vt= svd.matrixV().transpose();
		/**/
		return 0;
	}

	int PostProcessor::computeInteriorStresses(
		vect3d_map& t, vect3d_map& s, vect3d_map& clPt, const std::vector<Eigen::Vector3d> p,
		const vector_type& u, const vector_type& q, double E, double nu
	){
		// pre-fill result maps with zeros
		t.clear(); s.clear(); clPt.clear();
		for(int pi=0;pi<p.size();++pi){
			t[pi]=Eigen::Vector3d::Zero();
			s[pi]=Eigen::Vector3d::Zero();
			clPt[pi]=Eigen::Vector3d::Zero();
		}
		// loop over evaluation points
#pragma omp parallel for schedule(static,32)
		for(int pi=0;pi<p.size();++pi){
			Eigen::Matrix3d S; S.setZero();
			bool snap=false; double clDistSq=DBL_MAX;
			// loop over elements
			for(elem_map::iterator ei=elems.begin(); ei!=elems.end(); ++ei)
			if(cracks.count(regions[ei->first])==0){ //only loop over non-crack elements
				Eigen::Vector3d a,b,c,n;
				double ar,dsq,esq; // element area, squared distance to centroid, squared length of longest edge
				unsigned int q_order;
				// fetch node coordinates
				copyNode(ei->second[0], a);
				copyNode(ei->second[1], b);
				copyNode(ei->second[2], c);
				// fetch nodal displacements and element traction
				Eigen::Vector3d
					ua(u[3*(ei->second[0]-NODE_BASE_INDEX)], u[3*(ei->second[0]-NODE_BASE_INDEX)+1], u[3*(ei->second[0]-NODE_BASE_INDEX)+2]),
					ub(u[3*(ei->second[1]-NODE_BASE_INDEX)], u[3*(ei->second[1]-NODE_BASE_INDEX)+1], u[3*(ei->second[1]-NODE_BASE_INDEX)+2]),
					uc(u[3*(ei->second[2]-NODE_BASE_INDEX)], u[3*(ei->second[2]-NODE_BASE_INDEX)+1], u[3*(ei->second[2]-NODE_BASE_INDEX)+2]),
					qe(q[3*(ei->first-ELEM_BASE_INDEX)],q[3*(ei->first-ELEM_BASE_INDEX)+1],q[3*(ei->first-ELEM_BASE_INDEX)+2]);
				// compute area and normal
				n = (b-a).cross(c-a);
				ar= 0.5*n.norm();
				n.normalize();
				// compute distance
//				dsq = ((a+b+c)/3.0 - p[pi]).squaredNorm();
				double s_,t_; bool f;
				dsq = sqDistPtTri(p[pi],a,b,c,s_,t_,f);
				esq = std::max( (b-a).squaredNorm(), std::max( (c-a).squaredNorm(), (c-b).squaredNorm() ) );
				if( snap || dsq<0.01*esq ){ // too close, snap to element (0.09 --> 1/3 edge length)
					snap=true;
					if( dsq < clDistSq ) computeElementStresses(S, a,b,c, ua,ub,uc, qe, E,nu); // snap to the overall closest element
				}else
				{
					// choose integration order
					if(      dsq < esq   ) q_order = 20; // closer than edge length
					else if( dsq < 9*esq ) q_order = 16; // 1 to 3 edge lengths
					else if( dsq <25*esq ) q_order =  9; // 3 to 5 edge lengths
					else                   q_order =  3; // further than 5 edge lengths
					for(unsigned int qi=0; qi<getTriangleGaussNumPoints(q_order); ++qi){ // loop over integration points
						//fetch integration point
						Eigen::Vector3d qp, u_qp,dr;
						double S1,S2,r;
						getTriangleGaussPoint(q_order, qi, qp); // qp now has barycentric coords
						u_qp = ua*qp[0] + ub*qp[1] + uc*qp[2];
						qp.noalias() = a*qp[0] + b*qp[1] + c*qp[2]; // qp now has cartesian coords
						// evaluate stress kernels
						dr = (qp - p[pi]);
						r  = dr.norm();
						dr/= r;
						for(int i=0;i<3;++i) for(int j=i;j<3;++j) for(int k=0;k<3;++k){
							S1 = 1.0/(8.0*M_PI*(1.0-nu)*r*r)*( (1.0-2.0*nu)*( (k==j)*dr[i]+(k==i)*dr[j]-(i==j)*dr[k] ) + 3.0*dr[i]*dr[j]*dr[k] );
							S2 = E/(8.0*M_PI*(1.0-nu*nu)*r*r*r)*(
								 3.0*dr.dot(n)*( (1.0-2.0*nu)*(i==j)*dr[k]+nu*(j==k)*dr[i]+nu*(i==k)*dr[j]-5.0*dr[i]*dr[j]*dr[k] )
								 +3.0*n[k]*(1.0-2.0*nu)*dr[i]*dr[j]+n[i]*( (1.0-2.0*nu)*(j==k)+3.0*nu*dr[j]*dr[k] )
								 +n[j]*( (1.0-2.0*nu)*(i==k)+3.0*nu*dr[i]*dr[k] ) -n[k]*(i==j)*(1.0-4.0*nu) );
							// sum up results
							S(i,j) += ar*getTriangleGaussWeight(q_order,qi)*( S1*qe[k]-S2*u_qp[k] );
						}
					}
				}
				if( dsq < clDistSq ){
					clDistSq = dsq; // store squared distance to closest surface point
					clPt[pi] = (1.0-s_-t_)*a+s_*b+t_*c; // also store the location of this point
				}
			}
			// tension stresses (xx, yy, zz)
			t[pi][0]=S(0,0);
			t[pi][1]=S(1,1);
			t[pi][2]=S(2,2);
			// shear stresses (xy, xz, yz)
			s[pi][0]=S(0,1);
			s[pi][1]=S(0,2);
			s[pi][2]=S(1,2);
		}
		return 0;
	}

	int PostProcessor::findClosestPointInPlane(Eigen::Vector3d& cp,unsigned int node, const Eigen::Vector3d& n){
		// we want to find the closest point on the regular surface to the given node of the mesh
		// along the curve of the intersection of the regular surface with a plane specified
		// by the normal vector n and the given node
		double closestDistSq=DBL_MAX;
		Eigen::Vector3d a,b,c, p,q,r;
		double d, s, t;
		copyNode(node, p);

		for(elem_map::iterator ei=elems.begin(); ei!=elems.end(); ++ei)
		if( cracks.count(regions[ei->first])==0){ //only loop over non-crack elements
			copyNode(ei->second[0], a);
			copyNode(ei->second[1], b);
			copyNode(ei->second[2], c);
			if( (a-p).dot(n) * (b-p).dot(n) <= 0.0 ||
				(a-p).dot(n) * (c-p).dot(n) <= 0.0 ||
				(b-p).dot(n) * (c-p).dot(n) <= 0.0	){ // test if element intersects the given plane
				//re-order nodes such that ab and ac intersect the plane but bc does not
				if( std::abs((a-p).dot(n))<10*DBL_EPSILON ){ q=a; a=b; b=q; }
				if( std::abs((a-p).dot(n))<10*DBL_EPSILON ){ q=a; a=c; c=q; }
				if((a-p).dot(n) * (b-p).dot(n) > 0.0){ q=a; a=c; c=q; } else // a and b are on the same side --> swap a and c
				if((a-p).dot(n) * (c-p).dot(n) > 0.0){ q=a; a=b; b=q; }      // a and c are on the same side --> swap a and b
				// otherwise b and c are on the same side --> nothing to do
				// q = intersection of ab with plane (p,n)
				s = (p-b).dot(n) / (a-b).dot(n);
				q = s*a + (1.0-s)*b;

				// r = intersection of ac with plane (p,n)
				t = (p-c).dot(n) / (a-c).dot(n);
				r = t*a + (1.0-t)*c;

				// compute distance
				s = (p-r).dot(q-r) / (q-r).squaredNorm(); //dot(q-r,q-r);
				if( s > 1.0 ) s=1.0; else
				if( s < 0.0 ) s=0.0;
				q = s*q + (1.0-s)*r;
				d = (p-q).squaredNorm();
				if( d < closestDistSq ){
					cp = q;
					closestDistSq = d;
					//printf("\n%% nd %d: d=%.3lf , out-of-plane %.3lf", node, sqrt(d), (cp-p).dot(n));
				}
			}
		}
		return 0;
	}
}
