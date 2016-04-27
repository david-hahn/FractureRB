/* 
 * File:   SubsampledCrackTip.h
 * Author: David
 *
 * Created on 24. Feb. 2014, 13:49
 */


#ifndef SUBSAMPLEDCRACKTIP_H
#define	SUBSAMPLEDCRACKTIP_H

#include "types.h"
#include <string>

namespace FractureSim{
	class VDBWrapper;
	class FractureModel;

	class SubsampledCrackTip{
	public:
		SubsampledCrackTip(
			node_map& nodes_, elem_map& elems_, id_map& regions_,
			elem_map& crackTips_, elem_map& crackTipParents_, state_map& crackTipStates_,
			VDBWrapper* levelSet_, FractureModel* fractureModel_
		):  nodes(nodes_), elems(elems_), regions(regions_),
            crackTips(crackTips_), crackTipParents(crackTipParents_),
            crackTipStates(crackTipStates_), startScale(1.0)
		{
			levelSet=levelSet_;
			fractureModel=fractureModel_;
            outFileName.clear();
			sifOut=NULL;
		}

		// initialize from BEM mesh
		int init(double targetMeshSize, unsigned int nSamples, bool substepOutput=false);

        /* Propagate the crack based on SIF values at the crack tip
         * Returns additional nodes and elements in the mesh.
         * Also returns the new crack-tip and parent elements maps.
         * newCrackTips also contains the old nodes of the crack tip element
         * as additional entries in the value vector - should be truncated when
         * not needed anymore. Format is (new_a, new_b, old_a, old_b), even if
         * new_a==old_a or  new_b==old_b
         * The crack-tip states will be UPDATE for all propagated tip-elements,
         * ACTIVE_A or ACTIVE_B if they were in one of these states before,
         * ACTIVE for all not-propagated but ACTIVE tip-elements, and
         * INACTIVE for all tip-elements that were INACTIVE before
         */
		int propagate(
			node_map& newNodes, elem_map& newElems, id_map& newRegions,
			elem_map& newCrackTips, elem_map& newParents, state_map& newCrackTipStates,
			vect3d_map& crackTipFaceNormals, vect3d_map& crackTipTangents,
			vect3d_map& nodeSIFs, unsigned int nSubSteps
		);
        /* for debugging: write vertices to vtk file
         */
        int vtkWrite(std::string fileName);

		/* Create a new crack anchored at p with face normal n1 and in-plane "forward" normal n2
		 * newElems will be cleared and filled with additional elements for the caller's information
		 * the new geometry will be added to the mesh
		 */
		int startCrack(
			elem_map& newElems, Eigen::Vector3d p, Eigen::Vector3d n1,
			Eigen::Vector3d n2, unsigned int region
		);
        
        /* set the filename for substep VDB output
         * only has an effect if initialized with substep output
         * sequence number will be appended to the given name
         */
        inline void setFilename(std::string name){
            outFileName=name;
        }
		
		/* can be used to make crack seeds smaller or larger, default is 1.0
		 * set this before calling startCrack
		 */
		inline void setStartScale(double s){startScale=s;}

		/* Old version:
		 * before advancing vertices, make sure that states (of vertices)
         * has data consistent with crackTipStates
         */
        void updateVertexStates();

		inline unsigned int getVertsPerSegment(){return vertsPerSegment;}

		/* NEW FOR estimated SIF fractures
		 * Run through the high-res crack-front states and coarsify to mesh resolution
		 * A mesh node is active if there exists an active vertex on any of it's adjacent edges
		 * Parameter map will be cleared first
		 * Returns the number of active nodes (coarse) or negative on error
		 */
		int getCoarseCrackFrontStates(state_map& coarseStates);
		
		/* NEW FOR FractureRB:
		 * Copy marker data from source to this object
		 * sourceEdges lists which edges (wrt. the source's node numbering)
		 * should be copied (if found in the source object)
		 * targetNodes maps the source's to this object's node numbering
		 * if resetState==true all marker states will be decided based on
		 * whether they are inside of the level set's object grid
		 * otherwise the state from the source will be copied
		 */
		int copyEdges(SubsampledCrackTip& source, edge_imap sourceEdges, id_map targetNodes, bool resetState=false);
		
		/* If set, SIFs of node-vertices will be written at the beginning of every step
		 * SIFs are assumed to be piecewise constant over time and piecewise linear over edges
		 * set to "" to suppress output (default state)
		 */
		void setSIFsOutputFile(std::string filename);
		/* for debugging: add additional SIF data to the SIF output file
		 */
		void addToSIFsOutput(vect3d_map& nodeSIFs, std::string prefix);

		virtual ~SubsampledCrackTip();
	protected:
        static const double collapseThr; // collapse if edge is shorter than collapseThr*meshSize
        static const double subdivThr  ; // subdiv   if edge is longer  than subdivThr*meshSize
		double startScale; // can be used to make crack seeds smaller or larger, default is 1.0
		double meshSize;
        bool substepOutput;
        std::string outFileName;
        unsigned int vertsPerSegment;
		vect3d_map verts; // subsampled vertices along the crack-tip
        state_map  states; // store whether a vertex can propagate
        elem_map   neighbors; // store left and right neighbors of verts
        elem_map   vertNodes; // store which ct nodes are used to interpolate data
        node_map   nodeWeights; // store weights of nodes for interpolation
		id_map	   ctRegions; // store the region-ID from the parent triangle of the crack-tip element
        id_map     nodeVerts; //map ct nodes (old-id) to corresponding verts
        id_set     movedVerts; //store which verts have moved in the last step
		vect3d_map direction; //remember the most recent propagation direction of each vertex
		vect3d_map dir_min, dir_max; //new for FractureModel::cpVersion==2: also store min and max valid directions

		//references to the BEM mesh
		node_map&  nodes;
		elem_map&  elems;
		id_map&    regions;
		elem_map&  crackTips;
		elem_map&  crackTipParents;
		state_map& crackTipStates;

		//pointers to objects managed by FractureBEM
		VDBWrapper*    levelSet;
		FractureModel* fractureModel;
        
        /* Initialize vertices along the given crack-tip elements
         * The given set of crack tips MUST form closed loops,
         * otherwise the neighbor-mapping will be inconsistent.
         */
        int generateVertices(elem_map& theseCrackTips, unsigned int nextVert, int sourceNode=-1);
        
        /* once all active vertices have ben advanced for a certain number of
         * substeps, construct the update to the (coarser) BEM mesh
         */
        int buildMeshUpdate(
            node_map& newNodes, elem_map& newElems, id_map& newRegions,
            elem_map& newCrackTips, elem_map& newParents, state_map& newCrackTipStates
        );
        /* collapse edges whose nodeVerts have moved too close to each other
         * this is done on all edges but only modifies vertices
         * and updates the nodeVerts map. No new nodes or elements are added to
         * the crack-tip or crack-surface at this point. These tasks are
         * done while building the mesh update.
         */
        int collapseShortEdges();
        /* check the length of the edge (nd_a,nd_b) in terms of their mapped
         * vertices, if it is too long, subdivide it. All meshing operations are
         * completed in place.
         * Returns whether subdivision was done.
         * newA, newB, nextNode and nextElem are byRef params
         * and will be updated, as well as all mesh-maps.
         * However, nextCtElem can't be modified right away, as the subdivision
         * might happen before the original ctElem was completely processed.
         * If a subdivision happened, nextCtElem must be incremented once later!
         */
        bool checkSubdivideEdge(
            unsigned int nd_a, unsigned int nd_b,
            int& newA, int& newB, unsigned int nextCtElem,
            unsigned int& nextNode, unsigned int& nextElem,
            node_map& newNodes, elem_map& newElems, id_map& newRegions,
            elem_map& newCrackTips, elem_map& newParents,
            state_map& newCrackTipStates, id_map& addNodeVerts
        );
        /* Subdivisions need some cleanup work once the mesh-update has been
         * done, such as to make sure nodes are numbered in the same order in 
         * which they appear in the elements array.
         */
        void subdivCleanup(
            unsigned int& nextCtElem, unsigned int& nextElem,
            node_map& newNodes, elem_map& newElems, elem_map& newCrackTips,
            id_map& newNodeIDs, id_map& newVertNodes, id_map& addNodeVerts
        );
        /* add the triangle (nd_a,nd_b,nd_c) to newElems with ID nextElem
         * add regionID in newRegions
         * if ctElem>=0 this triangle will be added as the sole parent of
         * the crack-tip element ctElem
         */
        int addTri(
            elem_map& newElems, id_map& newRegions, unsigned int regionID,
            elem_map& newParents, int nextElem, int ctElem,
            unsigned int nd_a, unsigned int nd_b, unsigned int nd_c
        );
        /* double the sampling density along the crack-tip element (nd_a,nd_b)
         * returns the vertex ID of the central vertex along the element
         */
        unsigned int increaseSampling(unsigned int nd_a, unsigned int nd_b);
        /*
         */
        void redistributeVertices(unsigned int vertA, unsigned int vertB);
        void crackTipSmoothing();
        /* Update the interpolation weights on a segment
         */
        void updateWeights(
            unsigned int vertA, unsigned int vertB, 
            unsigned int nodeA, unsigned int nodeB 
        );
        /* Check if vertex data structures are consistent
         */
        void consitencyCheck();
		/* Remove vertices along a crack-tip edge element (node-vertices are kept and will become direct neighbors)
		 */
		void removeEdgeVerts(unsigned int nd_a, unsigned int nd_b);
		void deleteVertex(unsigned int vertex);
		inline CRACK_STATE updateState(CRACK_STATE s){
			if( s==ACTIVE   ) return UPDATE  ;
			if( s==ACTIVE_A ) return UPDATE_A;
			if( s==ACTIVE_B ) return UPDATE_B;
			return s;
		}
		void updateCrackTipStates(elem_map& newCrackTips, state_map& newCrackTipStates);
        void setNodeCoords(std::vector<double>& c, unsigned int vertex);

		// some stuff for debug-output of applied SIFs
		std::ofstream* sifOut;
		std::map<std::string, unsigned int> sifOutSteps, sifOutValues;
	};
}

#endif
