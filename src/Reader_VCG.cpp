//VCG library stuff ...
#include <vcg/complex/complex.h>
#include <vcg/wrap/io_trimesh/import.h>
#include "types.h"

using namespace vcg;
using namespace tri;

namespace FractureSim{
	namespace localReaderVCG{
		class MyVertex;
		class MyEdge;
		class MyFace;
		struct MyUsedTypes: public UsedTypes<Use<MyVertex>::AsVertexType,Use<MyEdge>::AsEdgeType,Use<MyFace>::AsFaceType>{};
		class MyVertex  : public Vertex< MyUsedTypes,
			vertex::VFAdj,
			vertex::Coord3f,
			vertex::Normal3f,
			vertex::Mark,
			vertex::BitFlags  >{};
		class MyEdge : public Edge< MyUsedTypes> {};
		class MyFace    : public Face< MyUsedTypes,
			face::VFAdj,
			face::VertexRef,
			face::BitFlags > {};
		class MyMesh    : public vcg::tri::TriMesh<std::vector<MyVertex>, std::vector<MyFace> > {};
	}
	int readVCG(std::string filename, node_map& nodes, elem_map& elems){
		localReaderVCG::MyMesh mesh;
		nodes.clear(); elems.clear();
		int r = io::Importer<localReaderVCG::MyMesh>::Open(mesh,filename.c_str());

		id_map vertexId; unsigned int k=NODE_BASE_INDEX;
        for(localReaderVCG::MyMesh::VertexIterator vi=mesh.vert.begin(); vi!=mesh.vert.end(); ++vi) if(!vi->IsD()){
			nodes[k].assign(3,0.0);
			nodes[k][0]= vi->P()[0];
			nodes[k][1]= vi->P()[1];
			nodes[k][2]= vi->P()[2];
            vertexId[vi-mesh.vert.begin()]=k;
            ++k;
        }
        k=ELEM_BASE_INDEX;
        for(localReaderVCG::MyMesh::FaceIterator fi=mesh.face.begin(); fi!=mesh.face.end(); ++fi) if(!fi->IsD()){
			elems[k].assign(3,0);
			elems[k][0]=vertexId[tri::Index(mesh, fi->V(0))];
            elems[k][1]=vertexId[tri::Index(mesh, fi->V(1))];
            elems[k][2]=vertexId[tri::Index(mesh, fi->V(2))];
            ++k;
        }
		return r;
	}
}