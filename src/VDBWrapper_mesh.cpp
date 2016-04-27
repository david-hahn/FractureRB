/* 
 * File:   VDBWrapper.cpp
 * Author: David
 * 
 * Created on 15. April 2014
 * 
 * contains the volume to mesh conversion method for VDBWrapper
 */

#include "VDBWrapper.h"
#include "distPtToTri.h"

//#define _DEBUG_SEGMENT_
#include "mySegment.h"

#include <openvdb/tools/MeshToVolume.h>
#include <openvdb/tools/VolumeToMesh.h>
#include <openvdb/tools/VolumeToSpheres.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/Composite.h>

//VCG library stuff ...
//#include <vcg/complex/complex.h>
//#include <vcg/complex/algorithms/smooth.h>
//#include <vcg/complex/algorithms/local_optimization.h>
//#include <vcg/complex/algorithms/local_optimization/tri_edge_collapse_quadric.h>
//#include <vcg/complex/algorithms/update/normal.h>
//#include <vcg/wrap/io_trimesh/export_obj.h>
#include "vcgHelper.h"
#include <vcg/wrap/io_trimesh/export_stl.h>
#include <vcg/complex/algorithms/hole.h>
#include <map>


namespace FractureSim{
	using namespace myVCG;
	using namespace std;

//	namespace localVCG{
//		// from VCG example tridecimator.cpp
//		class MyVertex;
//		class MyEdge;
//		class MyFace;
//		struct MyUsedTypes: public UsedTypes<Use<MyVertex>::AsVertexType,Use<MyEdge>::AsEdgeType,Use<MyFace>::AsFaceType>{};
//		class MyVertex  : public Vertex< MyUsedTypes,
//			vertex::VFAdj,
//			vertex::Coord3f,
//			vertex::Normal3f,
//			vertex::Mark,
//			vertex::BitFlags  >{
//		public:
//			vcg::math::Quadric<double> &Qd() {return q;}
//		private:
//			math::Quadric<double> q;
//		};
//		class MyEdge : public Edge< MyUsedTypes> {};
//		typedef BasicVertexPair<MyVertex> VertexPair;
//		class MyFace    : public Face< MyUsedTypes,
//			face::VFAdj,
//			face::VertexRef,
//			face::BitFlags > {};
//		class VcgMesh    : public vcg::tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace> > {};
//		class MyTriEdgeCollapse: public vcg::tri::TriEdgeCollapseQuadric< VcgMesh, VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex>  > {
//		public:
//			typedef  vcg::tri::TriEdgeCollapseQuadric< VcgMesh,  VertexPair, MyTriEdgeCollapse, QInfoStandard<MyVertex>  > TECQ;
//			typedef  VcgMesh::VertexType::EdgeType EdgeType;
//			inline MyTriEdgeCollapse(  const VertexPair &p, int i, BaseParameterClass *pp) :TECQ(p,i,pp){}
//		};
//	}
//	using namespace localVCG;
//
//	void copyVDBtoVCG(vdb::tools::VolumeToMesh& mesher, VcgMesh& mesh); // local fcn prototype

	int VDBWrapper::mesh(node_map& nodes, elem_map& elems, int nTris, double adaptive, double offsetVoxels){
		//printf(" params: nTris=%d, adaptive=%.3lf, offset=%.3lf\n", nTris, adaptive, offsetVoxels);
		vdb::tools::VolumeToMesh* mesher =
            new vdb::tools::VolumeToMesh(0.0,adaptive); // first param is isosurface (which level-set), second is adaptivity of mesh in [0,1]
        nodes.clear(); elems.clear();
        node_map tmpNodes;
        
        (*mesher)(*objectGrid);
		if(nTris<0){ // use original resolution
            for(int i=0; i<mesher->pointListSize(); ++i){
                tmpNodes[i+NODE_BASE_INDEX].push_back(mesher->pointList()[i][0]);
                tmpNodes[i+NODE_BASE_INDEX].push_back(mesher->pointList()[i][1]);
                tmpNodes[i+NODE_BASE_INDEX].push_back(mesher->pointList()[i][2]);
            }
            int k=ELEM_BASE_INDEX;
            for(int i=0; i<mesher->polygonPoolListSize(); ++i){
                vdb::tools::PolygonPool &pool = mesher->polygonPoolList()[i];
                for(int j=0; j<pool.numTriangles(); ++j){
                    //swap ordering of nodes to get outside-normals
                    elems[k].push_back(pool.triangle(j)[0]+NODE_BASE_INDEX);
                    elems[k].push_back(pool.triangle(j)[2]+NODE_BASE_INDEX);
                    elems[k].push_back(pool.triangle(j)[1]+NODE_BASE_INDEX);
                    ++k;
                }
                for(int j=0; j<pool.numQuads(); ++j){
                    //swap ordering of nodes to get outside-normals
                    elems[k].push_back(pool.quad(j)[0]+NODE_BASE_INDEX);
                    elems[k].push_back(pool.quad(j)[3]+NODE_BASE_INDEX);
                    elems[k].push_back(pool.quad(j)[2]+NODE_BASE_INDEX);
                    ++k;
                    elems[k].push_back(pool.quad(j)[0]+NODE_BASE_INDEX);
                    elems[k].push_back(pool.quad(j)[2]+NODE_BASE_INDEX);
                    elems[k].push_back(pool.quad(j)[1]+NODE_BASE_INDEX);
                    ++k;
                }
            }
            delete mesher;
        }else{ // use different - assumed coarser - resolution
            // create a VCG mesh from the VDB mesher's output
            // run topology-preserving decimation
            // copy result to output params (nodes,elems)
            VcgMesh mesh;
            copyVDBtoVCG(*mesher,mesh);
            delete mesher;

            TriEdgeCollapseQuadricParameter qparams;
            // defaults are:
                //BoundaryWeight=.5;
                //CosineThr=cos(M_PI/2);
                //FastPreserveBoundary=false;
                //NormalCheck=false;
                //NormalThrRad=M_PI/2;
                //OptimalPlacement=true;
                //PreserveBoundary = false;
                //PreserveTopology = false;
                //QuadricEpsilon =1e-15;
                //QualityCheck=true;
                //QualityQuadric=false;
                //QualityThr=.1;
                //QualityWeight=false;
                //SafeHeapUpdate =false;
                //ScaleFactor=1.0;
                //ScaleIndependent=true;
                //UseArea=true;
                //UseVertexWeight=false;
            qparams.PreserveTopology = true;
			qparams.QualityQuadric   = true;
			qparams.QualityThr       = DBL_MAX;			

            vcg::LocalOptimization<VcgMesh> DeciSession(mesh,&qparams);
            DeciSession.Init<MyTriEdgeCollapse>();
            DeciSession.SetTargetSimplices(nTris);
            DeciSession.DoOptimization();
			// don't do these, they might destroy manifoldness
			//vcg::tri::Clean<VcgMesh>::MergeCloseVertex(mesh,0.5*voxelSize);
			//vcg::tri::Clean<VcgMesh>::RemoveDuplicateFace(mesh);
			// instead we want to clean up self-intersections by deleting the 1-ring of all elements
			// that are involved in self-intersections and applying hole-filling afterwards ...
			int count = removeSelfIntersectingFaces(mesh);
			if( count>0) printf("\n%% ... deleted %d faces", count);
			if(count!=0){ //faces deleted, fill hole ...
				//vcg::tri::Clean<VcgMesh>::RemoveDegenerateEdge(mesh);
				//vcg::tri::Clean<VcgMesh>::RemoveDegenerateFace(mesh);
				vcg::tri::UpdateTopology<VcgMesh>::VertexFace(mesh);
				vcg::tri::UpdateTopology<VcgMesh>::FaceFace(mesh);
				vcg::tri::UpdateFlags<VcgMesh>::FaceBorderFromFF(mesh);
				vcg::tri::Clean<VcgMesh>::RemoveNonManifoldFace(mesh);
				
				vcg::tri::UpdateTopology<VcgMesh>::VertexFace(mesh);
				vcg::tri::UpdateTopology<VcgMesh>::FaceFace(mesh);
				vcg::tri::UpdateFlags<VcgMesh>::FaceBorderFromFF(mesh);
				vcg::tri::Clean<VcgMesh>::RemoveSmallConnectedComponentsSize(mesh,6);
				
				vcg::tri::Clean<VcgMesh>::RemoveUnreferencedVertex(mesh);
				
				vcg::tri::UpdateTopology<VcgMesh>::VertexFace(mesh);
				vcg::tri::UpdateTopology<VcgMesh>::FaceFace(mesh);
				vcg::tri::UpdateFlags<VcgMesh>::FaceBorderFromFF(mesh);
				vcg::tri::UpdateNormal<VcgMesh>::PerVertexPerFace(mesh);
				count=vcg::tri::Hole<VcgMesh>::EarCuttingIntersectionFill<vcg::tri::SelfIntersectionEar<VcgMesh> >(mesh,INT32_MAX,false);
				printf(" then filled %d holes",count);
				if(count==0) return -1; // failed to fill hole, assume the mesh is not good and abort
			}

			// offset along the outward normal so that we have more of the detailed surface inside the decimated mesh
			vcg::tri::UpdateNormal<VcgMesh>::PerVertex(mesh); // compute area-weighted normals
			vcg::tri::UpdateNormal<VcgMesh>::NormalizePerVertex(mesh); // make them unit-vectors
#pragma omp parallel for schedule(guided, mesh.vert.size()/32)
			for(int i=0; i<mesh.vert.size(); ++i) if(!mesh.vert[i].IsD()){
				Eigen::Vector3d p,n;
				mesh.vert[i].P().ToEigenVector(p);
				mesh.vert[i].N().ToEigenVector(n);
				p+=offsetVoxels*voxelSize*n; // how should we choose the offset distance?
				mesh.vert[i].P().FromEigenVector(p);
			}

            //vcg::tri::io::ExporterOBJ<VcgMesh>::Save(mesh,"test.obj",0);

            id_map vertexId;
            unsigned int k=NODE_BASE_INDEX;
            for(VcgMesh::VertexIterator vi=mesh.vert.begin(); vi!=mesh.vert.end(); ++vi) if(!vi->IsD()){
                tmpNodes[k].push_back(vi->P()[0]);
                tmpNodes[k].push_back(vi->P()[1]);
                tmpNodes[k].push_back(vi->P()[2]);
                vertexId[vi-mesh.vert.begin()]=k;
                ++k;
            }
            k=ELEM_BASE_INDEX;
            for(VcgMesh::FaceIterator fi=mesh.face.begin(); fi!=mesh.face.end(); ++fi) if(!fi->IsD()){
                elems[k].push_back(vertexId[vcg::tri::Index(mesh, fi->V(0))]);
                elems[k].push_back(vertexId[vcg::tri::Index(mesh, fi->V(1))]);
                elems[k].push_back(vertexId[vcg::tri::Index(mesh, fi->V(2))]);
                ++k;
            }
            mesh.Clear();
        }
        // finally renumber mesh to be HyENA compatible
        id_map renumber; unsigned int k=NODE_BASE_INDEX;
        for(elem_map::iterator it=elems.begin(); it!=elems.end(); ++it){
            for(unsigned int j=0; j<3; ++j){
                if( renumber.count( it->second[j] )==0 )
                    renumber[it->second[j]]=k++;
                it->second[j]=renumber[it->second[j]];
            }
        }
        for(node_map::iterator it=tmpNodes.begin(); it!=tmpNodes.end(); ++it){
            nodes[renumber[it->first]]=it->second;
        }
		return elems.size();
	}

	int VDBWrapper::writeVisualMesh(
			node_map& nodes, elem_map& elems, id_map& regions, id_set& cracks,
			const vector_type& u, const vector_type& u_c, std::string filename,
            double adaptiveVDB, double adaptiveQuadric, bool visDisplace, bool visCOD,
			bool visClose, bool visOBJ, double segment, bool writePerSegment, double postDecimation
	){
		if(!nearTriGridInitialized) return -1;
		//overview:
		//      intersect the objectGrid with all crackGrids
		//      -- maybe run segmentation here?
		//      triangulate the result
		//      -- or perhaps run segmentation here? (not sure we want segmentation at all atm.)
		//      voxelize the given mesh to create a closest primitive index (perhaps on coarser resolution)
		//      interpolate displacements from the given mesh onto the fine mesh
		//      -- when close to a crack check on which side the vertex is and store the sign and 1/2 of the opening displacement
		//		-- run over the elements and "diffuse" the sign in the 1-ring for each vertex and each crack to reduce errors when deciding the sign
		//      write the output file
		vdb::FloatGrid::Ptr fracturedObject = objectGrid->deepCopy();
		for(std::map<unsigned int, vdb::FloatGrid::Ptr>::iterator it=crackGrids.begin();
			it!=crackGrids.end(); ++it){
			vdb::tools::csgIntersection(*fracturedObject, *(it->second)->deepCopy());
		}

		VcgMesh mesh;
		if(segment>=0.0){
			vdb::FloatGrid::Ptr copyOfObject = fracturedObject->deepCopy();
			std::vector<vdb::FloatGrid::Ptr> segments =
				mySegment( *copyOfObject, segment*voxelSize,true );
			printf(" have %u segments ... ", segments.size());
			for(int i=0; i<segments.size(); ++i){
				vdb::tools::VolumeToMesh mesher(0.0,clamp(adaptiveVDB,0.0,1.0));
				mesher(*segments[i]);
				copyVDBtoVCG(mesher, mesh);
				if(writePerSegment) outputSubMesh( i, mesh, nodes,elems,regions,cracks,u,u_c,filename,adaptiveQuadric,visDisplace,visCOD,visClose,visOBJ,postDecimation);
			}
			if(!writePerSegment)    outputSubMesh(-1, mesh, nodes,elems,regions,cracks,u,u_c,filename,adaptiveQuadric,visDisplace,visCOD,visClose,visOBJ,postDecimation);
		}else{
			vdb::tools::VolumeToMesh mesher(0.0,clamp(adaptiveVDB,0.0,1.0)); // first param is isosurface (which level-set), second is adaptivity of mesh in [0,1]
			mesher(*fracturedObject);
			copyVDBtoVCG(mesher, mesh);
			outputSubMesh(-1, mesh, nodes,elems,regions,cracks,u,u_c,filename,adaptiveQuadric,visDisplace,visCOD,visClose,visOBJ,postDecimation);
		}

		return 0;
	}

	int VDBWrapper::outputSubMesh( int subMeshId, VcgMesh& mesh,
			node_map& nodes, elem_map& elems, id_map& regions, id_set& cracks,
			const vector_type& u, const vector_type& u_c, std::string filename,
            double adaptive, bool visDisplace, bool visCOD, bool visClose, bool visOBJ, double postDecimation
	){
        if(adaptive>0.0){
            TriEdgeCollapseQuadricParameter qparams;
            qparams.PreserveTopology = true;
            qparams.QualityQuadric   = true;
            qparams.QualityThr       = DBL_MAX;			
            vcg::LocalOptimization<VcgMesh> DeciSession(mesh,&qparams);
            DeciSession.Init<MyTriEdgeCollapse>();
            DeciSession.SetTargetMetric(adaptive);
            DeciSession.DoOptimization();
        }
		vcg::tri::UpdateNormal<VcgMesh>::PerVertex(mesh); // compute area-weighted normals
		vcg::tri::UpdateNormal<VcgMesh>::NormalizePerVertex(mesh); // make them unit-vectors

        vect3d_map crackClosingDisp;
        if(visClose){ // compute crack-closing displacement (undo "fattening" in the level-set needed for segmentation and intersections)
			typedef vdb::tools::GridSampler<vdb::FloatGrid::ConstAccessor, vdb::tools::BoxSampler> fltBoxSampler;			
			std::map<unsigned int, std::vector< vdb::Vec3R > > pointMap;
			std::map<unsigned int, std::vector< float > >   distanceMap;
			std::map<int, unsigned int> gridMap, positionMap;
			double tmpDist, val; unsigned int gridID;
            for(int i=0; i<mesh.vert.size(); ++i) if(!mesh.vert[i].IsD()){
				tmpDist=-voxelSize;
				for(std::map<unsigned int, vdb::FloatGrid::Ptr>::iterator it=
					crackGrids.begin(); it!=crackGrids.end(); ++it){
					fltBoxSampler sampler( it->second->getConstAccessor(), *resXform );
					val = sampler.wsSample( vdb::Vec3d( mesh.vert[i].P()[0], mesh.vert[i].P()[1], mesh.vert[i].P()[2] ) );
					if( val > tmpDist){
						tmpDist=val; gridID=it->first;
					}
				}
				if( tmpDist > -0.2*voxelSize ){
					positionMap[i]=pointMap[gridID].size();
					gridMap[i]=gridID;
					pointMap[gridID].push_back(vdb::Vec3d( mesh.vert[i].P()[0], mesh.vert[i].P()[1], mesh.vert[i].P()[2] ));
				}
            }
			vdb::tools::ClosestSurfacePoint<vdb::FloatGrid> cspSearcher;
			for(std::map<unsigned int, std::vector< vdb::Vec3R > >::iterator it =
				pointMap.begin(); it != pointMap.end(); ++it){
				cspSearcher.initialize<vdb::util::NullInterrupter>(*signedCrackGrids[it->first]);
				distanceMap[it->first].clear();
				cspSearcher.searchAndReplace(it->second,distanceMap[it->first]);
			}
			for(int i=0; i<mesh.vert.size(); ++i) if(!mesh.vert[i].IsD()){
				if(gridMap.count(i)>0){
					crackClosingDisp[i].setZero();
					crackClosingDisp[i][0] = pointMap[gridMap[i]][positionMap[i]][0] - mesh.vert[i].P()[0];
					crackClosingDisp[i][1] = pointMap[gridMap[i]][positionMap[i]][1] - mesh.vert[i].P()[1];
					crackClosingDisp[i][2] = pointMap[gridMap[i]][positionMap[i]][2] - mesh.vert[i].P()[2];
					//if( distanceMap[gridMap[i]][positionMap[i]] > 1.2*halfDiag*voxelSize )
					//	crackClosingDisp[i]*= 1.2*halfDiag*voxelSize/distanceMap[gridMap[i]][positionMap[i]];
					//if( crackClosingDisp[i].norm() > 1.2*halfDiag*voxelSize )
					//	crackClosingDisp[i]*= 1.2*halfDiag*voxelSize/crackClosingDisp[i].norm();
				}
			}

		}
        if(visDisplace){
            //vcg::tri::io::ExporterOBJ<VcgMesh> writerObj;
            //writerObj.Save(mesh,(filename+"_undef.obj").c_str(),vcg::tri::io::Mask::IOM_VERTNORMAL);
            // interpolate displacement at all nodes in the hi-res mesh
            // then offset all nodes to world space coordinates
            typedef pair<char, Eigen::Vector3d> codPair; // store sign and COD
            typedef map<unsigned int, codPair> codMap; // store region and COD with sign
            map<int, codMap> vertCOD, vertCODaggregate;
    #pragma omp parallel
            {
                int tri; bool isCrk; char sign;
                Eigen::Vector3d p,n,u_p,tmp;
                IdSetGrid::ConstAccessor nearTri = nearTriGrid->getConstAccessor();
                map<int, codMap> vertCOD_thread; // writing vertCOD directly is NOT THREAD-SAFE!
    #pragma omp for schedule(guided, mesh.vert.size()/32)
                for(int i=0; i<mesh.vert.size(); ++i) if(!mesh.vert[i].IsD()){
                    mesh.vert[i].P().ToEigenVector(p);
                    mesh.vert[i].N().ToEigenVector(n);

                    vdb::Coord c_p = nearTriXform->worldToIndexCellCentered(vdb::Vec3d(p[0],p[1],p[2]));
                    /**/ // using also neighboring voxels ...
                    id_set nearTris = nearTri.getValue(c_p).get();
                    for(int k=0; k<6; ++k){ // 6 or 26 neighbors?
                        const id_set& neighTris=nearTri.getValue(c_p+vdb::util::COORD_OFFSETS[k]).get();
                        nearTris.insert( neighTris.begin(), neighTris.end() );
                    } /*/
                    const id_set& nearTris = nearTri.getValue(c_p).get(); /**/
                    if(nearTris.empty()){
                        /*printf("x");/**/ continue; // no triangles found
                    } 
                    double d,s,t; bool up;
                    typedef multimap<double, unsigned int>::value_type pair_d_u;
                    typedef multimap<double, unsigned int>::iterator   iter_d_u;
                    multimap<double, unsigned int> sqDistPerTri;
                    value_map sPerTri, tPerTri;
                    id_set doneRegions;
                    for(id_set::const_iterator it = nearTris.begin(); it!=nearTris.end(); ++it){
                        Eigen::Vector3d
                            a(nodes[elems[*it][0]][0], nodes[elems[*it][0]][1], nodes[elems[*it][0]][2]),
                            b(nodes[elems[*it][1]][0], nodes[elems[*it][1]][1], nodes[elems[*it][1]][2]),
                            c(nodes[elems[*it][2]][0], nodes[elems[*it][2]][1], nodes[elems[*it][2]][2]);
                        d=sqDistPtTri(p,a,b,c,s,t,up);
                        sqDistPerTri.insert(pair_d_u(d,*it));
                        sPerTri[*it]=s; tPerTri[*it]=t;
                        //printf("p(%+.2lf,%+.2lf,%+.2lf) t%5u: sqDist=%+.4lf st(%+.3lf,%+.3lf)\n",p[0],p[1],p[2],*it,d,s,t);
                    }
                    // first, take the continous part of the displacement field from the overall closest triangle
                    // read from u_c if it belongs to a crack surface, u otherwise
                    // then add the opening displacement for all crack surfaces
                    tri  = sqDistPerTri.begin()->second; // maybe consider some blending between near triangles of different regions?
                    isCrk= cracks.count(regions[tri]) > 0;
                    if( tri>0){
                        interpDisplacement(u_p, p,n, isCrk? u_c : u, tri,sPerTri[tri],tPerTri[tri], nodes,elems); // when taking the continous part off a crack, look in u_c
                    }else{ u_p.setZero(); /*printf("x");/**/ }

                    // now add the opening displacement for all crack surfaces
                    if(visCOD) for(iter_d_u it=sqDistPerTri.begin(); it!=sqDistPerTri.end(); ++it){
                        tri=it->second;
                        if( tri>0 && cracks.count(regions[tri])>0 && doneRegions.count(regions[tri])==0){
                            sign=interpDisplacement(tmp, p,n, u, tri,sPerTri[tri],tPerTri[tri], nodes,elems, regions[tri], nearTriXform->voxelSize()[0]); // get the opening displacement from the nearest crack
                            //u_p+=sign*tmp;
                            if(vertCOD_thread.count(i)==0) vertCOD_thread[i].clear();
                            vertCOD_thread[i][regions[tri]]=codPair(2*sign,tmp);
                            doneRegions.insert(regions[tri]);
                            //if(tmp.norm() > nearTriXform->voxelSize()[0]) printf(" large cod interp!\n");
                        }
                    }
                    mesh.vert[i].P().FromEigenVector(p+u_p);
                }
    #pragma omp critical
                { // reduce vertex cod and sign maps over threads
                    vertCOD.insert(vertCOD_thread.begin(), vertCOD_thread.end());
                }
                vertCOD_thread.clear();
            }
            int m_steps=2; // use m steps of 1-ring voting on the sign
            for(int m=1;m<=m_steps;++m){ 
                vertCODaggregate=vertCOD; // make a copy
                for(int i=0; i<mesh.face.size(); ++i) if(!mesh.face[i].IsD()){
                    int v[3];
                    v[0] = vcg::tri::Index(mesh, mesh.face[i].V(0));
                    v[1] = vcg::tri::Index(mesh, mesh.face[i].V(1));
                    v[2] = vcg::tri::Index(mesh, mesh.face[i].V(2));
        //            printf(" %d %d %d (%d %d %d)\n", v0, v1, v2,
        //                vertCOD.count(v0),vertCOD.count(v1),vertCOD.count(v2));
                    for(int k=0; k<3; ++k){
                        if( vertCOD.count(v[k])>0){
                            for(codMap::iterator it=vertCOD[v[k]].begin(); it!=vertCOD[v[k]].end(); ++it){
                                for(int l=0; l<3; ++l){
                                    if( k!=l && vertCOD.count(v[l])>0){
                                        if( vertCOD[v[l]].count( it->first )>0){
                                            vertCODaggregate[v[l]][it->first].first += it->second.first/2;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if(m<m_steps){ // renormalize before next iteration
                    vertCOD=vertCODaggregate;
                    for(map<int, codMap>::iterator i1 = vertCOD.begin(); i1!=vertCOD.end(); ++i1){
                        for(codMap::iterator i2 = i1->second.begin(); i2!=i1->second.end(); ++i2){
                            if(i2->second.first >0) i2->second.first= 2;
                            if(i2->second.first <0) i2->second.first=-2;
                        }
                    }
                }
            }
            vertCOD.clear();
            
            for(map<int, codMap>::iterator i1 = vertCODaggregate.begin(); i1!=vertCODaggregate.end(); ++i1){
                Eigen::Vector3d p; //printf(" v %d has (region, sign)\n", i1->first);
                for(codMap::iterator i2 = i1->second.begin(); i2!=i1->second.end(); ++i2){
                    //printf(" (%d, %d)\n" ,i2->first, i2->second.first);
                    mesh.vert[i1->first].P().ToEigenVector(p);
                    if(i2->second.first >0) p+=i2->second.second;
                    if(i2->second.first <0) p-=i2->second.second;
                    mesh.vert[i1->first].P().FromEigenVector(p);
                }
            }
        }
		Eigen::Vector3d p;
		if(visClose) for(vect3d_map::iterator it=crackClosingDisp.begin(); it!=crackClosingDisp.end(); ++it){
            mesh.vert[it->first].P().ToEigenVector(p);
            p+=0.9*it->second; // apply only 90% of the closing displacement to prevent coincident surfaces in the visual mesh
            mesh.vert[it->first].P().FromEigenVector(p);
        }

		//tri::Smooth<VcgMesh>::VertexCoordPlanarLaplacian(mesh,3,math::ToRad(15.0));

		if( postDecimation>0.0 ){
            TriEdgeCollapseQuadricParameter qparams;
            qparams.PreserveTopology = true;
            qparams.QualityQuadric   = true;
            qparams.QualityThr       = DBL_MAX;			
            vcg::LocalOptimization<VcgMesh> DeciSession(mesh,&qparams);
            DeciSession.Init<MyTriEdgeCollapse>();
            DeciSession.SetTargetMetric(postDecimation);
            DeciSession.DoOptimization();
        }

        vcg::tri::UpdateNormal<VcgMesh>::PerVertex(mesh); // compute area-weighted normals
        vcg::tri::UpdateNormal<VcgMesh>::NormalizePerVertex(mesh); // make them unit-vectors

		if(subMeshId>=0){
			stringstream sstr; sstr << filename << "_s";
			sstr.fill('0'); sstr.width(4);
			sstr << subMeshId; filename = sstr.str();
		}
        if(visOBJ){
            vcg::tri::io::ExporterOBJ<VcgMesh> writerObj;
            writerObj.Save(mesh,(filename+".obj").c_str(),/**/vcg::tri::io::Mask::IOM_NONE/*/vcg::tri::io::Mask::IOM_VERTNORMAL/**/);
        }else{
            vcg::tri::io::ExporterSTL<VcgMesh> writer;
            writer.Save(mesh,(filename+".stl").c_str(),true,vcg::tri::io::Mask::IOM_NONE);
        }
		mesh.Clear();
		return 0;
	}

	char VDBWrapper::interpDisplacement(
		Eigen::Vector3d& u_p, const Eigen::Vector3d p,const Eigen::Vector3d n, const vector_type& u,
		unsigned int tri, double s, double t, node_map& nodes, elem_map& elems, int crackRegion, double bgval
	){
		Eigen::Vector3d
			u0(u[3*(elems[tri][0]-NODE_BASE_INDEX)  ],u[3*(elems[tri][0]-NODE_BASE_INDEX)+1],u[3*(elems[tri][0]-NODE_BASE_INDEX)+2]),
			u1(u[3*(elems[tri][1]-NODE_BASE_INDEX)  ],u[3*(elems[tri][1]-NODE_BASE_INDEX)+1],u[3*(elems[tri][1]-NODE_BASE_INDEX)+2]),
			u2(u[3*(elems[tri][2]-NODE_BASE_INDEX)  ],u[3*(elems[tri][2]-NODE_BASE_INDEX)+1],u[3*(elems[tri][2]-NODE_BASE_INDEX)+2]),
			v0(nodes[elems[tri][0]][0], nodes[elems[tri][0]][1], nodes[elems[tri][0]][2]);
		Eigen::Vector3d
			e0(nodes[elems[tri][1]][0], nodes[elems[tri][1]][1], nodes[elems[tri][1]][2]),
			e1(nodes[elems[tri][2]][0], nodes[elems[tri][2]][1], nodes[elems[tri][2]][2]);
		e0-=v0; // e0 is now the edge-vector b-a
		e1-=v0; // e1 is now the edge-vector c-a
		v0-=p;  // v0 is now the vector a-p

		u_p = (1-s-t)*u0 + s*u1 + t*u2; //printf("(%.2lf,%.2lf)", s,t);

		if(crackRegion>=0){ // we are interpolating crack opening displacement!
            char sign=1;
            double dist = ((s*e0 + t*e1)+v0).norm();
			// decide whether we are on the positive or the negative side of the crack
			//if we are close to the crack, do this based on level-set-detailed information!
			vdb::Vec3d p_vdb = vdb::Vec3d(p[0], p[1], p[2]);
			vdb::Vec3d n_vdb = vdb::Vec3d(n[0], n[1], n[2]); n_vdb.normalize();
			vdb::Coord p_crd = resXform->worldToIndexNodeCentered(p_vdb);
			float sdfValue, sdfDiff;
			bool isOnValue = signedCrackGrids[crackRegion]->getConstAccessor().probeValue(p_crd,sdfValue);
			vdb::tools::GridSampler<vdb::FloatGrid::ConstAccessor, vdb::tools::BoxSampler> sampler(
				signedCrackGrids[crackRegion]->getConstAccessor(), *resXform
			);
			if(isOnValue){
				//sdfDiff = sampler.wsSample(p_vdb - 0.5*voxelSize*n_vdb) - sampler.wsSample(p_vdb); // n_vdb is normalized, since it is a SDF grad-norm should be 1, we take a finite difference at 0.5*voxelSize which should be the max. difference
				//if(std::abs(sdfDiff) < 0.125*voxelSize && std::abs(sdfValue) < voxelSize ){
				//	p_vdb -= n_vdb*voxelSize;
				//	p_crd = resXform->worldToIndexNodeCentered(p_vdb);
				//}
				////printf("\n sdf grad %.3le", sdfDiff);
				//if( std::abs(sdfDiff) > 0.25*voxelSize && std::abs(sdfValue) < voxelSize){
				//	//printf(" n (%.3lf, %.3lf, %.3lf) d %.3le\n", n_vdb[0],n_vdb[1],n_vdb[2], sdfDiff);
				//	if( sdfDiff < 0.0 ) u_p*=-1.0; // use orientation of normals if close to the crack surface and normals are sufficiently aligned
				//}else{ // otherwise use SDF value
					if( std::abs(sdfValue) >= (NBHW-1)*voxelSize ){
						if( sdfValue < 0.0 ) sign=-1;//u_p*=-1.0; // easy case: in the narrow-band but not close to the crack-surface
						//printf("\n* sdf %+.3le",sdfValue);
					}else{ // close to the crack-surface: use box-sampling in the grid (don't do this close to the background value though!)
						sdfValue = sampler.wsSample(p_vdb); 
						if( sdfValue < 0.0 ) sign=-1;//u_p*=-1.0;
					}
				//}
			}else{ //not in the crack's narrow-band, use coarse mesh information
				Eigen::Vector3d tri_n = e0.cross(e1);
				if(v0.dot(tri_n) < 0.0) sign=-1;//u_p*=-1.0; // this is ok
			}

			u_p*=0.5; // now we have half-cod          
            // finally blend cod over dist
            dist-=halfDiag*voxelSize;
            if(dist<0) dist=0.0;
            if(dist>bgval) dist=bgval; // clamp
            //double  k = bgval*bgval; // n=2; k=bg^n;
            //double c0 = (k - sqrt(k*(k + 4)))/(2*k);
            //double c1 = (sqrt(k*(k + 4)) - k)/2;
			//double  f = c0 + 1/( c1 + dist*dist ); // inverse-quadratic blend
			double f = 1.0 - (dist/bgval); // linear blend
			if(u_p.squaredNorm() > bgval*bgval) u_p*=bgval/u_p.norm(); //clamp magnitude of COD to cutoff distance to avoid foldover
            u_p*=f; //printf("\nu_p blend f=%.3lf  d=%.5lf",f,dist);
			return sign;
		}
        return 0;
	}
    
//	void copyVDBtoVCG(vdb::tools::VolumeToMesh& mesher, VcgMesh& mesh){
//        // add vertices
//		unsigned int firstVert = mesh.vert.size();
//        VcgMesh::VertexIterator vi = vcg::tri::Allocator<VcgMesh>::AddVertices(mesh,mesher.pointListSize());
//        for(unsigned int i=0; i<mesher.pointListSize(); ++i,++vi){
//            vi->P()[0] = mesher.pointList()[i][0];
//            vi->P()[1] = mesher.pointList()[i][1];
//            vi->P()[2] = mesher.pointList()[i][2];
//        }
//        // count how many tri-faces we have
//        unsigned int nTrisOri = 0;
//        for(int i=0; i<mesher.polygonPoolListSize(); ++i){
//            vdb::tools::PolygonPool& pool = mesher.polygonPoolList()[i];
//            nTrisOri += pool.numTriangles();
//            nTrisOri += 2*pool.numQuads();
//        }
//        // add triangles
//        VcgMesh::FaceIterator fi = vcg::tri::Allocator<VcgMesh>::AddFaces(mesh,nTrisOri);
//        for(int i=0; i<mesher.polygonPoolListSize(); ++i){
//            vdb::tools::PolygonPool &pool = mesher.polygonPoolList()[i];
//            for(int j=0; j<pool.numTriangles(); ++j){
//                //swap ordering of nodes to get outside-normals
//                fi->V(0)=&(mesh.vert[firstVert+ pool.triangle(j)[0]]);
//                fi->V(1)=&(mesh.vert[firstVert+ pool.triangle(j)[2]]);
//                fi->V(2)=&(mesh.vert[firstVert+ pool.triangle(j)[1]]);
//                ++fi;
//            }
//            for(int j=0; j<pool.numQuads(); ++j){
//                //swap ordering of nodes to get outside-normals
//                fi->V(0)=&(mesh.vert[firstVert+ pool.quad(j)[0]]);
//                fi->V(1)=&(mesh.vert[firstVert+ pool.quad(j)[2]]);
//                fi->V(2)=&(mesh.vert[firstVert+ pool.quad(j)[1]]);
//                ++fi;
//                fi->V(0)=&(mesh.vert[firstVert+ pool.quad(j)[0]]);
//                fi->V(1)=&(mesh.vert[firstVert+ pool.quad(j)[3]]);
//                fi->V(2)=&(mesh.vert[firstVert+ pool.quad(j)[2]]);
//                ++fi;
//            }
//        }
//	}

	//****************************************************************************
	int VDBWrapper::addToNearTriGrid(node_map& nodes, elem_map& elems, double res){
        using namespace vdb;
		if(!nearTriGridInitialized){
			nearTriGrid = IdSetGrid::create();
			nearTriGridInitialized=true;
			//create the transform for the given resolution
            nearTriXform = vdb::math::Transform::createLinearTransform(res);
		}
		//"voxelize" ...
    	IdSetGrid::Accessor acc = nearTriGrid->getAccessor();
        Coord ijk, n_ijk; const Vec3d half(0.5,0.5,0.5);
        for(elem_map::iterator it=elems.begin(); it!=elems.end(); ++it){
            Vec3d
                a(nodes[it->second[0]][0],nodes[it->second[0]][1],nodes[it->second[0]][2]),
                b(nodes[it->second[1]][0],nodes[it->second[1]][1],nodes[it->second[1]][2]),
                c(nodes[it->second[2]][0],nodes[it->second[2]][1],nodes[it->second[2]][2]);
            a = nearTriXform->worldToIndex(a);
            b = nearTriXform->worldToIndex(b);
            c = nearTriXform->worldToIndex(c);
            BBoxd elemBBox; elemBBox.expand(a); elemBBox.expand(b); elemBBox.expand(c);
            ijk = util::nearestCoord(elemBBox.getCenter());
//            printf("bbox of tri %d has min (%.3lf,%.3lf,%.3lf) and max (%.3lf,%.3lf,%.3lf)\n", it->first,
//                elemBBox.min()[0],elemBBox.min()[1],elemBBox.min()[2], elemBBox.max()[0],elemBBox.max()[1],elemBBox.max()[2]);
//            printf("  ---  center is (%.3lf,%.3lf,%.3lf) with nearest voxel (%d,%d,%d)\n",
//                elemBBox.getCenter()[0],elemBBox.getCenter()[1],elemBBox.getCenter()[2], ijk[0],ijk[1],ijk[2]);

            deque<Coord> coordList;
            coordList.push_back(ijk);
            while(!coordList.empty()){
                ijk = coordList.back(); coordList.pop_back();
                if(!acc.isValueOn(ijk)){ //initialize
                    acc.setValue(ijk, IdSetAdapter());
                }
                acc.getValue(ijk).get().insert(it->first);
                for(int i=0; i<26; ++i){
                    n_ijk = ijk + util::COORD_OFFSETS[i];
                    BBoxd voxelBBox(n_ijk.asVec3d()-half, n_ijk.asVec3d()+half); // bounding box of the voxel (cell-centered: coord is center of a cube)
                    if(voxelBBox.hasOverlap(elemBBox)){
                        if(acc.isValueOn(n_ijk)){
                            if(acc.getValue(n_ijk).get().count(it->first)==0){
                                coordList.push_back(n_ijk);
                            }
                        }else{
                            coordList.push_back(n_ijk);
                        }
                    }
                }
            }
        }
		return 0;
	}
	
	// NEW FOR FractureRB
	vdb::BBoxd VDBWrapper::getNearTris(const vdb::Vec3d& p, id_set& tris){
		vdb::Coord c_p = nearTriXform->worldToIndexCellCentered(p);
		IdSetGrid::ConstAccessor nearTri = nearTriGrid->getConstAccessor();
		id_set& nearTris = nearTri.getValue(c_p).get();
		// copy to output param
		tris.clear();
		tris.insert(nearTris.begin(), nearTris.end());
		
		vdb::Vec3d
			p_c_p = nearTriXform->indexToWorld(c_p),
			halfVoxel = nearTriXform->voxelSize()*0.5;
		vdb::BBoxd voxelBBoxWorld(
			p_c_p-halfVoxel,
			p_c_p+halfVoxel
		);
		
//		//debug output
//		printf("near-tri lookup:\n"
//			   " :: input point \t(%.3lf, %.3lf, %.3lf)\n"
//			   " :: mapped coord\t(%d, %d, %d)\n"
//			   " :: in world    \t(%.3lf, %.3lf, %.3lf)\n"
//			   " :: bbox min    \t(%.3lf, %.3lf, %.3lf)\n"
//			   " :: bbox max    \t(%.3lf, %.3lf, %.3lf)\n",
//				p[0],p[1],p[2],
//				c_p[0],c_p[1],c_p[2],
//				p_c_p[0],p_c_p[1],p_c_p[2],
//				voxelBBoxWorld.min()[0],voxelBBoxWorld.min()[1],voxelBBoxWorld.min()[2],
//				voxelBBoxWorld.max()[0],voxelBBoxWorld.max()[1],voxelBBoxWorld.max()[2]
//		);
		return voxelBBoxWorld;
	}

	// NEW FOR FractureRB
	unsigned int VDBWrapper::findClosestSurfaceTri(const Eigen::Vector3d& p,
		node_map& nodes, elem_map& elems, id_map& regions, id_set& cracks
	){
	    vdb::Coord c_p = nearTriXform->worldToIndexCellCentered(vdb::Vec3d(p[0],p[1],p[2]));
		IdSetGrid::ConstAccessor nearTri = nearTriGrid->getConstAccessor();
		id_set nearTris = nearTri.getValue(c_p).get();
        for(int k=0; k<26; ++k){ // 6 or 26 neighbors?
            const id_set& neighTris=nearTri.getValue(c_p+vdb::util::COORD_OFFSETS[k]).get();
            nearTris.insert( neighTris.begin(), neighTris.end() );
        }
        if(nearTris.empty()){
			return ELEM_BASE_INDEX-1; // no triangles found
        } 
        double d,s,t, minDist=1.0;
		bool outside, first=true;
		unsigned int closestTri=ELEM_BASE_INDEX-1;
        for(id_set::const_iterator it = nearTris.begin(); it!=nearTris.end(); ++it){
			if( cracks.count(regions[*it])==0 ){ // not a crack element
				Eigen::Vector3d
					a(nodes[elems[*it][0]][0], nodes[elems[*it][0]][1], nodes[elems[*it][0]][2]),
					b(nodes[elems[*it][1]][0], nodes[elems[*it][1]][1], nodes[elems[*it][1]][2]),
					c(nodes[elems[*it][2]][0], nodes[elems[*it][2]][1], nodes[elems[*it][2]][2]);
				d=sqDistPtTri(p,a,b,c,s,t,outside);
				if(first || d < minDist){
					minDist=d;
					closestTri= *it;
					first=false;
				}
			}
        }
		return closestTri;
	}

	// NEW FOR FractureRB
	std::vector<vdb::FloatGrid::Ptr> VDBWrapper::getSegments(double handleThreshold, bool useTiles) const{
		std::vector<vdb::FloatGrid::Ptr> result;
		try{
			// need to work on copies, otherwise we'll lose data
			//printf("%% copying grid ...\n");
			vdb::FloatGrid::Ptr theGrid = objectGrid->deepCopy(), currentCrack;
			//printf("%% ... have %d active voxels\n",theGrid->activeVoxelCount());
			// intersect object with fractures
			for(std::map<unsigned int, vdb::FloatGrid::Ptr>::const_iterator it=
				crackGrids.begin(); it!=crackGrids.end(); ++it
			){
				currentCrack= it->second->deepCopy();
				//printf("%% intersecting with crack %d ...\n", it->first);
				for(vdb::FloatTree::RootNodeType::ChildOnIter cit = objectGrid->tree().beginRootChildren(); cit.test(); ++cit){
					currentCrack->getAccessor().touchLeaf(cit.getCoord()); // make sure the crack-grid covers the object grid
					currentCrack->getAccessor().setValueOn(cit.getCoord());
					// if we find a background value, the sign needs to be flipped
					if( currentCrack->getAccessor().getValue(cit.getCoord()) == currentCrack->background()/* > crackGrids[i]->voxelSize()[0]*/ ){
						currentCrack->getAccessor().setValue(cit.getCoord(), -currentCrack->background());
						// this should not happen, but if it did, we just fixed it
						// still, let's print a warning
						printf("\n%% !!! WARNING: VDBWrapper::getSegments: %s has inconsistent tile coverage!\n", currentCrack->getName().c_str());
					}
				}
				currentCrack->signedFloodFill();
				//Note: in hindsight it might have been more efficient to not flip signs on the crack-grids (ever) and instead use a modified csgIntersection method
				vdb::tools::csgIntersection(*theGrid, *currentCrack );
				//printf("%% ... have %d active voxels\n",theGrid->activeVoxelCount());
			}

			//printf("%% segmenting ...\n");
			result = mySegment(*theGrid,handleThreshold*voxelSize,useTiles);

			// idea: could dilate all segments by halfDiag*voxelSize then intersect with the original surface to close cracks (faster than projection?)
		}catch(std::exception& e){
			printf("exception in VDBWrapper::getSegments: %s\n", e.what());
			result.clear();
		}catch(...){
			printf("caught ellipsis in VDBWrapper::getSegments\n");
			result.clear();
		}
		//printf("%% segmentation done\n");
		return result;
	}
}
