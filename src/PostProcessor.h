#ifndef POSTPROCESSOR_H
#define	POSTPROCESSOR_H

#include "types.h"

namespace FractureSim{
	class FractureModel;

    enum VTK_CELL_TYPE{
        VTK_VERTEX          = 1,
        VTK_POLY_VERTEX     = 2,
        VTK_LINE            = 3,
        VTK_POLY_LINE       = 4,
        VTK_TRIANGLE        = 5,
        VTK_TRIANGLE_STRIP  = 6,
        VTK_POLYGON         = 7,
        VTK_PIXEL           = 8,
        VTK_QUAD            = 9,
        VTK_TET             = 10,
        VTK_VOXEL           = 11,
        VTK_HEXAHEDRON      = 12 
    };
    
    class PostProcessor{
    public:
        PostProcessor(
            node_map& nodes_, elem_map& elems_, id_map& regions_, id_set& cracks_,
            elem_map& crackTips_, elem_map& parents_, state_map& crackTipStates_,
			value_map& crackAreas_, double crackMeshSize_=0.0
        ): nodes(nodes_), elems(elems_), regions(regions_), cracks(cracks_),
           crackTips(crackTips_), parents(parents_), compressive(1.0),
           crackTipStates(crackTipStates_), crackAreas(crackAreas_), crackMeshSize(crackMeshSize_)
		{}
		
		PostProcessor(
			PostProcessor& ori
		): nodes(ori.nodes), elems(ori.elems), regions(ori.regions),
           cracks(ori.cracks), crackTips(ori.crackTips), parents(ori.parents), compressive(1.0),
           crackTipStates(ori.crackTipStates), crackAreas(ori.crackAreas), crackMeshSize(ori.crackMeshSize)
		{}

        virtual ~PostProcessor(){}
        
        /* Compute stress intensity factors in local coordinate system for all
         * crack-tip nodes, using displacement correlation technique
         */
		int computeNodeSIFs(
            const vector_type& displacements, double E, double nu,
			vect3d_map& sifs, vect3d_map& faceNormals, vect3d_map& tangents
        );
		
        /* Estimate stress intensity factors in local coordinate system for all
         * crack-tip nodes, using the regular stress field and the estimated
		 * fracture size (does not require an up-to-date BEM solution)
		 * only the regular surface displacements and tractions are used
         */
		int estimateNodeSIFs(
            const vector_type& displacements, const vector_type& tractions, double E, double nu,
			vect3d_map& sifs, vect3d_map& faceNormals, vect3d_map& tangents, state_map& nodeStates
        );

		/*
		 * When using INITIAL_BEM method, estimate COD by inverse DCT on estimated SIFs
		 */
		int estimateCOD(
			const elem_map& newElems, vect3d_map& nodeSIFs, const id_set& ctNodes,
			double E, double nu, vect3d_map& estCOD
		);

		/* Clears stored VTK data entries.
         * Data arrays pointed to will NOT be deleted!
         */
        inline void vtkClearData(){
            vtkData.clear();
            vtkDataDimension.clear();
            vtkDataIsCellData.clear();
            vtkDataNames.clear();
        }

		/* Add point or cell data to the VTK output
		 */
        int vtkAddData(std::string name, int dimension, bool isCellData, vector_type& data);

		/* Write mesh and all associated data to a VTK file
		 */
        int vtkWrite(std::string fileName, VTK_CELL_TYPE cellType);

		/* Compute stresses on triangles in the mesh based on nodal displacements
		 * and add to VTK data
         * output-param retValues       (terminology of "node_map" is misused here)
         * is a map<int, vector<double> > that maps element-IDs to 9 double values:
         * the first value is the max. principal stress value,
		 * the next 3 values are the plane-normal across which the principal stress is given;
		 * the following 4 values are magnitude & normal vector for the min. principal stress
		 * the last value is a flag: 0 for regular surface elements, >0 for fracture elements
		 * for fractures, we specify whether max. and min. principal stress reach their largest
		 * magnitude for the positive or negative side of the fracture (sign of applied COD):
		 * 1: both positive; 2: max. negative, min. positive; 3: max. pos., min. neg.; 4: both neg.
		 */
		int computeSurfaceStresses(
            node_map& retValues, const vector_type& displacements,
			const vector_type& crackBaseDisplacements,
			const vector_type& tractions,
            double E=1.0, double nu=0.0, bool addToVTK=false, int ignoreElems=-1
        );

		/* Compute interior stresses S (tensions holds S_xx, S_yy, S_zz, while
		 * shears holds S_xy, S_xz, and S_yz) at each location in points
		 * based on the boundary displacements and tractions using stress kernels.
		 * since this requires distance computation to all elements on the regular surface
		 * (not on fractures) we also report the closest point on the surface for each evaluation point
		 * See [Kielhorn2009] Appendix 1 for details.
		 */
		int computeInteriorStresses(
			vect3d_map& tensions, vect3d_map& shears, vect3d_map& closestPoint, const std::vector<Eigen::Vector3d> points,
			const vector_type& displacements, const vector_type& tractions,
			double E=1.0, double nu=0.0
		);
		
		/* In the triangle given by the node-IDs in el, 
		 * find the node-ID which is neither nd_a nor nd_b
		 */
        inline unsigned int findInteriorNode(
            unsigned int nd_a, unsigned int nd_b,
            const std::vector<unsigned int>& el
        ){
            if(el[0]!=nd_a && el[0]!=nd_b)
                return el[0];
            else if(el[1]!=nd_a && el[1]!=nd_b)
                return el[1];
            else //if(el[2]!=nd_a && el[2]!=nd_b)
                return el[2];
        }

		/* Copy the coordinates of the specified node to coords
		 */
        inline void copyNode(unsigned int node, Eigen::Vector3d& coords){
            coords[0]=nodes[node][0];
            coords[1]=nodes[node][1];
            coords[2]=nodes[node][2];
        }

		inline void setCompressiveFactor(double c){ compressive=c; }

    protected:
        node_map& nodes;
        elem_map& elems;
        id_map& regions;
        id_set& cracks;
        elem_map& crackTips;
        elem_map& parents;
        state_map& crackTipStates;
		value_map& crackAreas;
        double crackMeshSize;
		double compressive; // compressive toughness = this factor * tensile toughness

        std::vector<std::string>    vtkDataNames;
        std::vector<int>            vtkDataDimension;
        std::vector<bool>           vtkDataIsCellData;
        std::vector<vector_type*>	vtkData;
		vector_type s_xx, s_yy, s_zz, s_xy, s_xz, s_yz; // cartesian stresses per node/element

		vect3d_map lastCrackTipFaceNormals, lastCrackTipTangents; // store the local coordinate system for the last known crack-tip (either computeSIFs or estimateSIFs called)
		
		void vtkWriteData(std::ofstream& out, int i, int n);
        
        /* Build a local coordinate frame used to compute SIFs
         * Input are 3 points a,b,c
         * Output are 3 unit vector_types n1,n2,n3, where
         * n1 is the face normal of the triangle (a,b,c),
         * n2 is the edge normal of the edge (a,b) in-plane, and
         * n3 is the tangent unit vector along the edge (a,b)
         */
        inline void getLocalCoordFrame(
            const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c,
            Eigen::Vector3d& n1, Eigen::Vector3d& n2, Eigen::Vector3d& n3
        ){
            // n3 is the edge unit vector_type --> mode III direction (sliding)
            n3 = (b-a); n3.normalize();
            // n1 is the face normal --> mode I direction (opening)
            n1 = n3.cross(c-a); n1.normalize();
            // n2 is the in-plane edge normal --> mode II direction (shear)
            n2 =-n1.cross(n3); n2.normalize(); // normalization should be obsolete
            // now we have a local coordinate system for the edge a-b
        }
		
		// get averaged tangents and out-of-plane normals at crack-tip nodes
		int getCrackTipAxis(vect3d_map& normals, vect3d_map& tangents);
		
		// computes the deformation gradient F = X B^-1 of a surface element
		int computeElementDefGrad(
			Eigen::Matrix3d& F,
			const Eigen::Vector3d&  a, const Eigen::Vector3d&  b, const Eigen::Vector3d&  c,
			const Eigen::Vector3d& ua, const Eigen::Vector3d& ub, const Eigen::Vector3d& uc
		);
		// computes cartesian stresses based on the deformation gradient F
		// and the linear material model S = 2 mu (F-I) + lambda tr(F-I);
		int computeElementStresses(
			Eigen::Matrix3d& S,
			const Eigen::Vector3d&  a, const Eigen::Vector3d&  b, const Eigen::Vector3d&  c,
			const Eigen::Vector3d& ua, const Eigen::Vector3d& ub, const Eigen::Vector3d& uc,
			const Eigen::Vector3d&  q, double E, double nu
		);
		// computes SVD of the deformation gradient F = U*S*Vt
		// and principal stresses P = 2 mu (S-I) + lambda tr(S-I);
		int computeElementPrincipalStresses(
			Eigen::Matrix3d& U, Eigen::Matrix3d& S, Eigen::Matrix3d& Vt, Eigen::Matrix3d& P,
			const Eigen::Vector3d&  a, const Eigen::Vector3d&  b, const Eigen::Vector3d&  c,
			const Eigen::Vector3d& ua, const Eigen::Vector3d& ub, const Eigen::Vector3d& uc,
			const Eigen::Vector3d&  q, double E, double nu
		);
		// find the closest point to the specified node on the regular surface within a plane
		// defined by the normal vector n and the given node
		int findClosestPointInPlane(Eigen::Vector3d& cp,unsigned int node, const Eigen::Vector3d& n);
    };   
}

#endif
