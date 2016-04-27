#ifndef HYENAWRAPPER_H
#define	HYENAWRAPPER_H

#include "types.h"

#ifdef _MSC_VER
	#ifndef HYENA_DLL
		#define HYENA_DLL __declspec(dllimport)
	#endif
#else
	// leave empty as GCC exports all symbols by default
	#define HYENA_DLL
#endif

namespace FractureSim{

    class HYENA_DLL HyENAWrapper{
    public:
        static HyENAWrapper* constructHyENAWrapper(
            node_map& nodes_, elem_map& elems_, id_map& regions_, id_set& cracks_,
            double youngsModulus, double poissonsRatio, bool is2D
        );
		virtual ~HyENAWrapper(){}

		/* Initialize the BEM matrices:
		 * first create a HyENAWrapper object
		 * - the mesh data is stored by ref! make sure the maps contain valid data
		 * than specify your boundary conditions using buildBoundaryCnd
		 * finally call init, nodes on the crack-tip are treated with special boundary conditions
		 * ToDo: (?) store (by ref) crack_tip_nodes as well to make things a bit more efficient and less confusing
		 */
        virtual int init(id_set& crack_tip_nodes, bool ignoreCracks=false) =0;
		/* Update the BEM matrices:
		 * first add (do not remove or change existing data!) new geometry to the mesh maps
		 * then call addCrack, again as in init the set of nodes on the crack-tip must be given
		 */
        virtual int addCrack(id_set& crack_tip_nodes) =0;
		/* Solve the BEM system and merge known and unknown data
		 * use getDisplacements and getTractions to access the result
		 */
        virtual int compute() =0;
		// just for testing ...
		virtual int innerEvalTest(std::string filename) =0;
        virtual void innerEval(vector_type& result, const vector_type& coords) =0;
		/* The standard displacement data (u) contain absolute displacement only for non-crack boundaries
		 * whereas we solve for crack-opening-displacements on fracture surfaces.
		 * This function evaluates the continous part of the displacement (u_c) field on cracks.
		 * The two sides of an opening crack can then be computed as u_c+1/2u and u_c-1/2u respectively.
		 */
		virtual vector_type& computeCrackBaseDisplacements() =0;
		
		/*
         * Build map of boundary conditions from a vector of strings, each entry follows the format
         * "boundary-id(f_x,f_y,f_z)", where f_i are forces along x,y,z direction (double) or
         * "boundary-id[d_x,d_y,d_z]", where d_i are displacements along x,y,z
         * or "boundary-id(fixed)" locks all dofs of nodes on boundary
         * or "boundary-id(crack)" applies crack-surface boundary conditions (COD computed)
         * boundary-ids will be read from file with same name as input-mesh and extension .boundary
         * i.e. mesh.boundary, which is expected to be in Elmer-mesh format.
         */
		virtual int buildBoundaryCnd(const std::vector<std::string>& specs) =0;

        /* Add a crack boundary condition for the given boundary-id
         * Returns 0 on success
         * or -1 if the boundary-id already had a BC (no change/overwrite)
         */
        virtual int addCrackBoundaryCnd(unsigned int bnd_id) =0;

		/* NEW FOR FractureRB:
		 * re-assemble the right hand side with updated boundary data
		 * the type of the BC can't be changed though
		 */
		virtual int updateRHS(vect3d_map& newBcData, id_set& crack_tip_nodes) =0;

		/* NEW FOR method INIT_BEM:
		 * update the right hand side with estimated opening displacements
		 */
		virtual int updateRHSwithCOD(vect3d_map& newCOD, elem_map& support_elems, id_set& support_crack_tip_nodes) =0;

		//getters
		inline vector_type& getDisplacements(){ return displacements; }
		inline vector_type& getTractions(){ return tractions; }
		inline node_map& getNodes(){ return nodes; }
        inline elem_map& getElems(){ return elems; }
		inline id_map& getRegions(){ return regions; }
        inline id_set& getCracks(){ return cracks; }
		inline double getYoungsModulus(){ return youngsModulus; }
		inline double getPoissonsRatio(){ return poissonsRatio; }

	protected:
		HyENAWrapper(
            node_map& nodes_, elem_map& elems_, id_map& regions_, id_set& cracks_,
            double youngsModulus, double poissonsRatio, bool is2D
        ) : nodes(nodes_), elems(elems_),
            regions(regions_), cracks(cracks_)
        {
            this->youngsModulus=youngsModulus;
            this->poissonsRatio=poissonsRatio;
			this->is2D=is2D;
        }
        // member variables
        double youngsModulus, poissonsRatio;
		bool is2D; //no longer supported

		// these maps will be dynamically updated as the geometry changes
        node_map& nodes;
        elem_map& elems;
        id_map& regions;
        id_set& cracks; //set of region IDs which are cracks
		vector_type displacements, tractions; // solution data

    };   
}

#endif
