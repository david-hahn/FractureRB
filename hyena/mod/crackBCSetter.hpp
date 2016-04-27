#ifndef CRACKBCSETTER_HPP_
#define CRACKBCSETTER_HPP_

#include <set>
#include <map>

//#include "hyena/core/dofs/bcsetter.hpp"
#include "hyena/core/io/bcdata.hpp"
#include "hyena/core/traits/traits.H"
#include "hyena/core/dofs/dofhandleraccessor.hpp"

namespace hyena
{

  //struct BCData;

  template<ELEMENT_SHAPE E,
           APPROXIMATION G,
           APPROXIMATION F,
           PROBLEM P,
           SPACE_TYPE S>
  class CrackBCSetter
    : public DofHandlerAccessor<E,G,F,P,S>
  {
    typedef DofTraits<E,G,F,P,S>   dof_traits_type;
    typedef TriaTraits<E,G>       tria_traits_type;

    typedef typename dof_traits_type::dof_handler_type dof_handler_type;
    typedef typename dof_traits_type::gdof_type  gdof_type;
    typedef typename dof_traits_type::ldof_type  ldof_type;
    typedef typename dof_traits_type::superelement_type  superelement_type;
    typedef typename dof_traits_type::gdof_array_type gdof_array_type;
    typedef typename dof_traits_type::ldof_array_type ldof_array_type;
    typedef typename dof_traits_type::superelement_array_type
    superelement_array_type;

//    typedef BCData  bc_data_type;
    
    typedef typename ProblemTraits<P>::value_type value_type;

    typedef typename tria_traits_type::element_type element_type;
    typedef typename tria_traits_type::node_type node_type;
    enum {num_ldofs_per_gdof = ProblemTraits<P>::num_ldofs_per_gdof,
          num_gdofs_per_element =
          ReferenceElement<E, F>::num_gdofs_per_element};

  private:
    const std::map<unsigned int, BCData>& bc_input_;
    const std::map<unsigned int, unsigned int>& region_id_;
    const std::set<unsigned int>& crack_tip_nodes_;
    const BC_TYPE bc_type_;
	bool is2d_;
    BCData defaultBC;
    std::map<unsigned int, unsigned int>& cu_ldofs_IDtoIDX_;
    
  public:
    CrackBCSetter(const dof_handler_type& dof_handler,
        const std::map<unsigned int, BCData>& bc_input,
        const BC_TYPE bc_type,
        const std::map<unsigned int, unsigned int>& region_id,
        const std::set<unsigned int>& crack_tip_nodes,
        std::map<unsigned int, unsigned int>& cu_ldofs_IDtoIDX, bool is2d=false
    ) : DofHandlerAccessor<E,G,F,P,S>(dof_handler),
        bc_input_(bc_input),
        bc_type_(bc_type),
        region_id_(region_id),
        crack_tip_nodes_(crack_tip_nodes),
        cu_ldofs_IDtoIDX_(cu_ldofs_IDtoIDX)
    {
        is2d_=is2d;
        // default BC is homogeneous Neumann
        defaultBC.bc_type=NEUMANN;
        defaultBC.homogenous=true;
    }


    template<typename VECTOR>
    void operator()(VECTOR& known_data) const
    {
      VECTOR dummy_1, dummy_2;
      setBC(known_data, dummy_1, dummy_2);
    }
/**/
    template<typename VECTOR>
    void operator()(VECTOR& known_data, VECTOR& robin_data_1,
                    VECTOR& robin_data_2) const
    {
      setBC(known_data, robin_data_1, robin_data_2);
    }

    // Initializes all ldofs to ANY_LDOF_TYPE
    void operator()() const
    {
      unsigned int ldof_cnt=0;
      typename ldof_array_type::const_iterator ldof_it =
           this->getLDofs().begin();
      for(; ldof_it !=  this->getLDofs().end(); ++ldof_it, ++ldof_cnt)
      {
        (*ldof_it)->setIDX(ldof_cnt);
        (*ldof_it)->setType(ANY_LDOF_TYPE);
      }
    }
/**/
  private:
 
    template<typename VECTOR>
    void setBC(VECTOR& known_data, VECTOR& robin_data_1,
               VECTOR& robin_data_2) const;

  };

#include "crackBCSetter.tpl"

} // end namespace hyena

#endif /* CRACKBCSETTER_HPP_ */
