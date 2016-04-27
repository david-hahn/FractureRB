// Copyright (C) 2009-2010 Matthias Messner, Michael Messner, Franz
// Rammerstorfer, Peter Urthaler
//
// This file is part of HyENA - a C++ boundary element methods library.
//
// HyENA is free software: you can redistribute it and/or modify it under the
// terms of the GNU Lesser Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.
//
// HyENA is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU Lesser Public License for more
// details.
//
// You should have received a copy of the GNU Lesser Public License along with
// HyENA. If not, see <http://www.gnu.org/licenses/>.

/**
 * based on:
 * @file    innerevaluation.hpp
 * @author  Michael
 * @date    created:     12.07.11
 *          last change: 11.12.12
 * modified by DH
 */
#ifndef INNEREVALUATION_H_INCLUDED
#define INNEREVALUATION_H_INCLUDED

#include <float.h>

#include "hyena/mod/innerintegrator3d.hpp"
#include "src/distPtToTri.h"

#include "hyena/core/collocation/collokernel.hpp"
#include "hyena/core/traits/traits.H"


namespace hyena
{

  // forward declarations
  template<typename>  class InnerAssembler;
  template<typename, typename, typename, typename>  class ColloKernel;

  

  /**
   * InnerEvaluation provides the evaluation of an boundary integral operator
   * at an interior point. This is needed for either an indirect ansatz or for
   * the evaluation of the representation formula in the context of a direct
   * ansatz.
   *
   * @tparam row_traits_type
   * @tparam row_traits_type
   * @tparam col_traits_type
   * @tparam fsol_type
   * @tparam kernel_tag_type
    */
  template <typename row_traits_type,
            typename col_traits_type,
            typename fsol_type,
            typename kernel_tag_type>
  class InnerEvaluation : boost::noncopyable
  {
    // typedefs --------------------------------------------------------
    // public:
    // protected:
  private:
    typedef typename row_traits_type::dof_handler_type row_handler_type;
    typedef typename col_traits_type::dof_handler_type col_handler_type;

    typedef typename row_traits_type::const_ldof_array_type row_ldof_array_type;
    typedef typename col_traits_type::const_ldof_array_type col_ldof_array_type;

    typedef ColloKernel<row_traits_type, col_traits_type, fsol_type, 
                        kernel_tag_type> kernel_type;


    // member variables ------------------------------------------------
    // public:
    // protected:
  private:
    fsol_type& fsol_;
    const InnerIntegrator3d<kernel_type> integrate_;
    InnerAssembler<col_traits_type> assemble_;

    static const unsigned int nldgd_
    = ProblemTraits<col_traits_type::problem>::num_ldofs_per_gdof;



    // ctors -----------------------------------------------------------
  public:
    /**
     * ctor initializes the integration objects required for the evalutation of
     * the representation formula.
     *
     * @param[in] fsol fundamental solution
     */
    InnerEvaluation(fsol_type & fsol)
      : fsol_(fsol),
        integrate_(fsol_),
        assemble_()
    {  }

    // protected:
    // private:




    // member functions ------------------------------------------------
  public:
    template< typename DATA_VECTOR,
              typename POINT_VECTOR>
    void operator()(const row_ldof_array_type &ldofs,
                    const DATA_VECTOR &data,
                    const POINT_VECTOR &inner_points,
                    const double& sign,
                    DATA_VECTOR &res)
    {
      typedef typename POINT_VECTOR::const_iterator const_ip_iterator_type;

      const_ip_iterator_type ip_start = inner_points.begin();
      const_ip_iterator_type ip_end   = inner_points.end();

      res.resize(inner_points.size()*nldgd_);
      res.setZero();

      assemble_(ip_start, ip_end, ldofs, data, integrate_, res, sign);
    }

    //protected:
    //private:

  };




  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////


  /**
   * @tparam traits_type
   */
  template<typename traits_type>
  class InnerAssembler
  {
    // typedefs ---------------------------------------------------------
    // public:
    // protected:
  private:

    typedef typename traits_type::const_ldof_array_type y_ldof_array_type;
    typedef typename traits_type::superelement_extractor_type
    y_superelement_extractor_type;
    typedef typename traits_type::const_superelement_array_type
    y_const_superelement_array_type;
    typedef typename traits_type::gdof_type y_gdof_type;



    // static const & enums --------------------------------------------
    // public:
    // protected:
  private:
    static const unsigned int num_ldofs_per_gdof_
    = ProblemTraits<traits_type::problem>::num_ldofs_per_gdof;



    // ctors -----------------------------------------------------------
  public:
    /**
     * constructor initializes private superelement extractor
     */
    InnerAssembler()
      : y_se_extractor_()
    { }

    // protected:
    // private:



    // member functions ------------------------------------------------
  public:
    /**
     *
     */
    template< typename IP_ITER,
              typename DATA_VECTOR,
              typename INTEGRATOR >
    void operator()(const IP_ITER x_begin,
                    const IP_ITER x_end,
                    const y_ldof_array_type &y_ldofs,
                    const DATA_VECTOR &y_data,
                    const INTEGRATOR& integrate,
                    DATA_VECTOR &res,
                    const double sign=1.0) const
    {
      typedef typename INTEGRATOR::kernel_type kernel_type;

      typedef ReferenceElement<kernel_type::shape,
                               kernel_type::approx_y> y_type;

      y_const_superelement_array_type y_se;
      y_se_extractor_(y_ldofs, y_se);


#ifdef ASSEMBLE_INNER_OMP
#pragma omp parallel for schedule(guided)
#endif
      for (IP_ITER x_iter=x_begin; x_iter<x_end; ++x_iter) {

        typename kernel_type::kernel_block_type block;

        typename y_const_superelement_array_type::const_iterator
            y_se_it = y_se.begin(),
            closestTri;

		double       closestDist = DBL_MAX, cl_s, cl_t, tmp_s, tmp_t;
		bool         isSingular = false, cl_out, tmp_out;

        for (; y_se_it != y_se.end(); ++y_se_it) {

			// set block to zero
			for(unsigned int j=0; j < y_type::num_gdofs_per_element; ++j)
				block(0,j).zeros();

			//DH: estimate integration order based on proximity
			double distance_x2y = //( (*x_iter) - (*y_se_it)->getMidPoint() ).norm();
				sqDistPtTri(
                    *x_iter,
                    (*(*y_se_it)->getNode(0)),(*(*y_se_it)->getNode(1)),(*(*y_se_it)->getNode(2)),
                    tmp_s, tmp_t, tmp_out
                ); // assumes that we are using triangles in 3d

			//printf("\n%%   *  distance tri to point (%.4lf, %.4lf, %.4lf) for %d is %.6lf",
			//	(*x_iter)[0],(*x_iter)[1],(*x_iter)[2], (*y_se_it)->getElement()->getInputId(), sqrt(distance_x2y));
			
			if( distance_x2y < closestDist){
				closestDist = distance_x2y;
                cl_s = tmp_s; cl_t = tmp_t; cl_out = tmp_out;
				closestTri  = y_se_it;//->getElement()->getInputId();
			}
			

			// integration order chosen based on the proximity of the two superelements
			double h_y = (*y_se_it)->getMeshSize();
			distance_x2y = distance_x2y / (h_y*h_y);
			unsigned int order=10;
			if(      distance_x2y > 5.0*5.0 )	order = 3;
			else if( distance_x2y > 1.0     )	order = 6;
			else if( distance_x2y > 0.1*0.1 )	order = 9;
			else isSingular = true; //detect quasi-singular case (eval pt (*x_iter) is too close to triangle (*y_se_it) )
			if(isSingular) continue; // skip integration if we have a quasi-singular case
            // note that we do not want to break the loop, as we need to find the closest triangle as well (which might not be the first quasi-singular one)
            
			integrate(*x_iter, *y_se_it, block, order);

			// loop over all shape_functions on tau_y
			for(unsigned int ngdy=0; ngdy<y_type::num_gdofs_per_element; ++ngdy){
			const y_gdof_type* const y_gdof = (*y_se_it)->getGDof(ngdy);
			// loop over all num_ldofs_per_gdof on collo_x
			for(unsigned int nldx=0; nldx<num_ldofs_per_gdof_; ++nldx){
				const unsigned int i = distance(x_begin, x_iter);
				// loop over all num_ldofs_per_gdof on tau_y
				for(unsigned int nldy=0; nldy<num_ldofs_per_gdof_; ++nldy){
                    const unsigned int j = y_gdof->getLDof(nldy)->getID(); //DH ori was y_gdof->getLDof(nldy)->getIDX();
                    //printf("\n y_data[%d]=%.3le (idx %d)",j,y_data[j],y_gdof->getLDof(nldy)->getIDX());
                    res[i*num_ldofs_per_gdof_ + nldx]
                        += sign*block(0,ngdy)(nldx, nldy)*y_data[j];
				}
			  }
		    }
        } // end y_se_it

		//printf("\n%% closest tri to point (%.4lf, %.4lf, %.4lf) is %d with distance %.6lf", (*x_iter)[0],(*x_iter)[1],(*x_iter)[2], closestTri, sqrt(closestDist));
        if(isSingular || cl_out){
			if( traits_type::field_approximation == LINEAR ){ //interpolate PL in tri
				// this gives the correct result only when used for displacements as y_data! should check for DLP_tag here as well

				// inner integration ran into a quasi-singular case
				// or the evaluation point is outside (wrt the closest triangle's normal)
				// in which case the inner integration result would be 0
				// replace the result for this point by linear interpolation on the closest triangle
				// here we assume that y_type::num_gdofs_per_element == 3 (linear triangle)
				// and that num_ldofs_per_gdof_ == 3 (in 3d space)
				//printf("\n cl_tri interp for pt (%.4lf, %.4lf, %.4lf), dist %.3le, s,t (%.3lf, %.3lf)",
				//    (*x_iter)[0],(*x_iter)[1],(*x_iter)[2], sqrt(closestDist), cl_s,cl_t);
				Eigen::Vector3d u_x = interpInTri( cl_s, cl_t,
					y_data[(*closestTri)->getGDof(0)->getLDof(0)->getID()],
					y_data[(*closestTri)->getGDof(0)->getLDof(1)->getID()],
					y_data[(*closestTri)->getGDof(0)->getLDof(2)->getID()],
					y_data[(*closestTri)->getGDof(1)->getLDof(0)->getID()],
					y_data[(*closestTri)->getGDof(1)->getLDof(1)->getID()],
					y_data[(*closestTri)->getGDof(1)->getLDof(2)->getID()],
					y_data[(*closestTri)->getGDof(2)->getLDof(0)->getID()],
					y_data[(*closestTri)->getGDof(2)->getLDof(1)->getID()],
					y_data[(*closestTri)->getGDof(2)->getLDof(2)->getID()]
				);
				const unsigned int i = distance(x_begin, x_iter);
				res[3*i  ]= u_x[0]*=-sign;
				res[3*i+1]= u_x[1]*=-sign;
				res[3*i+2]= u_x[2]*=-sign;

				if(sign < 0.6 && sign > -0.6){ // we use |sign|=0.5 for crack opening displacements
					// need to check if the evaluation point is up or down of the element we snapped to
					double projection = (
                        (*(*closestTri)->getNode(1)) - (*(*closestTri)->getNode(0))
                    ).crossProduct(
                        (*(*closestTri)->getNode(2)) - (*(*closestTri)->getNode(0))
                    ).dotProduct(
                        (*x_iter)                    - (*(*closestTri)->getNode(0))
                    );
					if( projection < 0.0 ){
                        res[3*i  ]*=-sign;
                        res[3*i+1]*=-sign;
                        res[3*i+2]*=-sign;
                    }else{
                        res[3*i  ]*=sign;
                        res[3*i+1]*=sign;
                        res[3*i+2]*=sign;
                    }
				}
			}
			else{ //if( traits_type::field_approximation == CONSTANT ){
				const unsigned int i = distance(x_begin, x_iter);
				res[3*i  ]=0; res[3*i+1]=0; res[3*i+2]=0;
				//// this is wrong:
				//const unsigned int i = distance(x_begin, x_iter);
				//res[3*i  ]= y_data[(*closestTri)->getGDof(0)->getLDof(0)->getID()];
				//res[3*i+1]= y_data[(*closestTri)->getGDof(0)->getLDof(1)->getID()];
				//res[3*i+2]= y_data[(*closestTri)->getGDof(0)->getLDof(2)->getID()];
			}
        }

      } // end x_gdof_it
    }



  protected:
    // interpolate displacements in triangle with given (s,t) coords,
    // where (0,0) is node a, (1,0) is node b and (0,1) is node c.
	inline Eigen::Vector3d interpInTri( double s, double t,
        double u_a_x, double u_a_y, double u_a_z,
        double u_b_x, double u_b_y, double u_b_z,
        double u_c_x, double u_c_y, double u_c_z
    ) const{
		Eigen::Vector3d
            u_a(u_a_x, u_a_y, u_a_z),
            u_b(u_b_x, u_b_y, u_b_z),
            u_c(u_c_x, u_c_y, u_c_z);
        return (1.0-s-t)*u_a + s*u_b + t*u_c;
    }
    
  private:
    const y_superelement_extractor_type y_se_extractor_;

  };

} // end namespace hyena

#endif // INNEREVALUATION_H_INCLUDED
