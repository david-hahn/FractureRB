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
 * @file    innerintegrator3d.hpp
 * @author  Michael
 * @date    created:     12.07.11
 *          last change:
 * modified by DH
 */
#ifndef INNERINTEGRATOR_H_INCLUDED
#define INNERINTEGRATOR_H_INCLUDED




// own includes
#include "hyena/core/tria/element.hpp"
#include "hyena/core/tria/superelement.hpp"
#include "hyena/core/dofs/gdof.hpp"
#include "hyena/core/be/referenceelements.hpp"
#include "hyena/core/quadrature/quadrature.hpp"
#include "hyena/core/collocation/collokernel.hpp"
#include "hyena/core/common/enumerators.H"
#include "hyena/core/common/mat.hpp"
#include "hyena/core/traits/traits.H"



namespace hyena
{

  /**
   * @ingroup io
   *
   * Integration routine for inner point evaluation
   *
   * @tparam KERNEL kernel type
   */
  template<typename KERNEL>
  class InnerIntegrator3d : boost::noncopyable
  {
    // static const & enums --------------------------------------------
    // public:
    // protected:
    // private:



    // typedefs --------------------------------------------------------
  public:
    /**
     * Define the @p Kernel type which is given as template parameter at compile
     * time.
     */
    typedef KERNEL  kernel_type;

    // protected:
    // private:



    // member variables ------------------------------------------------
    // public:
    // protected:
  private:
    //!private kernel object
    const kernel_type kernel_;

    //! private quadrature object
    const QuadratureRule<kernel_type::shape, GAUSS> reg_quad_;




    // ctors -----------------------------------------------------------
  public:
    /**
     * This constructor initializes reg_quad_ rule with a dummy order (10) and
     * the kernel with fundsol. Finally the integration order is computed
     * depending on the distance of the inner point from the actual superelement
     * (which is numerically integrated over).
     *
     * @param[in] fundsol fundamental solution
     */
    template<typename FSOL>
    InnerIntegrator3d(FSOL& fsol)
      : kernel_(fsol),
        reg_quad_(10) // order specified here is never used, overridden by operator() param.
    { }


    // protected:
    // private:



    // member functions ------------------------------------------------
  public:
    /**
     * The overloaded operator() allows to use the InnerIntegrator3d in a
     * functor like manner. This is the core of the object. From here
     * singularityCheck, the right integration routine and the mapping to the
     * setting as defined by the duffy transformation and back are called.
     *
     * @tparam superelement_type Element type superelement_y is pointing to
     * @param[in] innerpoint_x innerpoint_point
     * @param[in] superelement_y support for ansatz functions
     * @param[in,out] result resulting matrix block
     */
    template<typename superelement_type>
    void operator()(const Point3  innerpoint_x,
                    const superelement_type *const superelement_y,
                    typename kernel_type::kernel_block_type& result,
					unsigned int order=10) const //caller should specify the actual integration order here
    {
      const ELEMENT_SHAPE shape_y = superelement_type::element_type::shape;
      const APPROXIMATION geom_y  = superelement_type::element_type::order;
	  
      // inner loop goes over all quadrature points
      for(unsigned int k=0; k<reg_quad_.getNumPoints(order); ++k)
        kernel_(reg_quad_.getPoint(order,k),
                innerpoint_x,
                superelement_y,
                ApproxTraits<shape_y, geom_y>::mapped_pos,
                reg_quad_.getWeight(order,k),
                result);
    }


    // protected:
    // private:

  };


} // end namespace hyena

#endif // INNERINTEGRATOR_H_INCLUDED
