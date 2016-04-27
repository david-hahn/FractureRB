#ifndef HYENAWRAPPERI_H
#define	HYENAWRAPPERI_H

#ifdef _MSC_VER
	// set this before including HyENAWrapper.h
	#define HYENA_DLL __declspec(dllexport)
#else
	// leave empty as GCC exports all symbols by default
	#define HYENA_DLL
#endif
#include "HyENAWrapper.h"

#include <set>

#include "hyena/core/common/enumerators.H"
#include "hyena/core/io/bcdata.hpp"

// fwd decl
namespace hyena{
    //struct BCData;
    template<ELEMENT_SHAPE SHAPE, APPROXIMATION ORDER>
	class Mesh;
	template<ELEMENT_SHAPE E, APPROXIMATION G, APPROXIMATION F, PROBLEM P, SPACE_TYPE S>
	class DofHandler;
	template<ELEMENT_SHAPE E, APPROXIMATION G, APPROXIMATION F, PROBLEM P,  SPACE_TYPE S>
	class LDof;
	template<ELEMENT_SHAPE E, APPROXIMATION G, APPROXIMATION F, PROBLEM P, SPACE_TYPE S>
	class Extractor;
}

namespace FractureSim{
    using namespace hyena;

    class HYENA_DLL HyENAWrapperImpl : public HyENAWrapper {
    public:
		HyENAWrapperImpl(
            node_map& nodes_, elem_map& elems_, id_map& regions_, id_set& cracks_,
            double youngsModulus, double poissonsRatio, bool is2D
        ) : HyENAWrapper(nodes_,elems_,regions_,cracks_,youngsModulus,poissonsRatio,is2D),
			mesh(NULL), uh(NULL), th(NULL)
        {}
		virtual ~HyENAWrapperImpl();

		/* Initialize the BEM matrices:
		 * first create a HyENAWrapper object
		 * - the mesh data is stored by ref! make sure the maps contain valid data
		 * than specify your boundary conditions using buildBoundaryCnd
		 * finally call init, nodes on the crack-tip are treated with special boundary conditions
		 * ToDo: (?) store (by ref) crack_tip_nodes as well to make things a bit more efficient and less confusing
		 */
        virtual int init(id_set& crack_tip_nodes, bool ignoreCracks=false);
		/* Update the BEM matrices:
		 * first add (do not remove or change existing data!) new geometry to the mesh maps
		 * then call addCrack, again as in init the set of nodes on the crack-tip must be given
		 */
        virtual int addCrack(id_set& crack_tip_nodes);
		/* Solve the BEM system and merge known and unknown data
		 * use getDisplacements and getTractions to access the result
		 */
        virtual int compute();
		// just for testing ...
		virtual int innerEvalTest(std::string filename);
        virtual void innerEval(vector_type& result, const vector_type& coords);
		/* The standard displacement data (u) contain absolute displacement only for non-crack boundaries
		 * whereas we solve for crack-opening-displacements on fracture surfaces.
		 * This function evaluates the continous part of the displacement (u_c) field on cracks.
		 * The two sides of an opening crack can then be computed as u_c+1/2u and u_c-1/2u respectively.
		 */
		virtual vector_type& computeCrackBaseDisplacements();
		
		/*
         * Build map of boundary conditions from a vector of strings, each entry follows the format
         * "boundary-id(f_x,f_y,f_z)", where f_i are forces along x,y,z direction (double) or
         * "boundary-id[d_x,d_y,d_z]", where d_i are displacements along x,y,z
         * or "boundary-id(fixed)" locks all dofs of nodes on boundary
         * or "boundary-id(crack)" applies crack-surface boundary conditions (COD computed)
         * boundary-ids will be read from file with same name as input-mesh and extension .boundary
         * i.e. mesh.boundary, which is expected to be in Elmer-mesh format.
         */
		virtual int buildBoundaryCnd(const std::vector<std::string>& specs);

        /* Add a crack boundary condition for the given boundary-id
         * Returns 0 on success
         * or -1 if the boundary-id already had a BC (no change/overwrite)
         */
        virtual int addCrackBoundaryCnd(unsigned int bnd_id);

		/* NEW FOR FractureRB:
		 * re-assemble the right hand side with updated boundary data
		 * the type of the BC can't be changed though
		 */
		virtual int updateRHS(vect3d_map& newBcData, id_set& crack_tip_nodes);

		/* NEW FOR method INIT_BEM:
		 * update the right hand side with estimated opening displacements
		 */
		virtual int updateRHSwithCOD(vect3d_map& newCOD, elem_map& support_elems, id_set& support_crack_tip_nodes);
		
        // HyENA configuration
        static const ELEMENT_SHAPE  E  = TRIANGLE;
        static const APPROXIMATION  G  = LINEAR;
        static const APPROXIMATION  DA = LINEAR;
        static const APPROXIMATION  NA = CONSTANT;
        static const PROBLEM        P  = REAL3_3D;
        static const SPACE_TYPE     DS = CONTINUOUS;
        static const SPACE_TYPE     NS = DISCONTINUOUS;

	protected:
		/*
         * Parse boundary condition
         * "boundary-id(f_x,f_y,f_z)", where f_i are forces along x,y,z direction (double)
         * or 'fixed', i.e. '1(fixed) locks all dofs of nodes on boundary 1
         * and 2(0,0,1.0) applies unit force in z-direction on boundary 2
         * boundary-ids will be read from file with same name as input-mesh and extension .bnd
         * i.e. mymesh.bnd, which is expected to be in Elmer-mesh format.
         */
        int parseBoundaryCnd(const char* spec, BCData& bc);
		std::map<unsigned int, BCData> boundaryConditions;

		// the mesh object hold geometry information (perhaps we don't need it here?)
        Mesh<E,G>* mesh; // using pointer here as the mesh object must be reconstructed for each geometry update
		// the dof handlers hold Dirichlet and Neumann data, indices, etc.
		DofHandler<E,G,DA,P,DS>* uh;
		DofHandler<E,G,NA,P,NS>* th;

		// HyENA LDof handling ...
		std::vector<const LDof<E,G,DA,P,DS>*> dk_ldofs, du_ldofs, cu_ldofs;
		std::vector<const LDof<E,G,NA,P,NS>*> nk_ldofs, nu_ldofs;
        id_map cu_IDtoIDX;

		matrix_type V, K, D, Kc, Dc, Dcc; // system matrix blocks
        vector_type fD, fN;  // right hand side blocks
		vector_type diri_data, neum_data; // known boundary data
		vector_type crackInnerEvals;
		vector_type test1,test2; // just for experimenting ...
        unsigned int M, N, L; // system size
		
		Eigen::PartialPivLU<matrix_type> luV, luS; // store factorizations in init()
    };   
}

#endif
