
#ifndef CRACKBCSETTER_TPL
#define CRACKBCSETTER_TPL

#include<vector>


template<ELEMENT_SHAPE E,
         APPROXIMATION G,
         APPROXIMATION F,
         PROBLEM P,
         SPACE_TYPE S>
template<typename VECTOR>
void CrackBCSetter<E,G,F,P,S>::
 setBC(VECTOR& known_data, VECTOR& robin_data_1,
      VECTOR& robin_data_2) const
{

  std::vector<double> known_data_tmp;
//  std::vector<double> robin_data_1_tmp;
//  std::vector<double> robin_data_2_tmp;
  
  typename ldof_array_type::const_iterator ldof_it =
     this->getLDofs().begin();
  
  bool known = false;
  unsigned known_idx_cnt=0;
  unsigned unknown_idx_cnt=0;
  unsigned cunknown_idx_cnt=0;
  unsigned homo_idx_cnt=0;
//  unsigned robin_idx_cnt=0;
  unsigned int ldof_cnt = 0;
  BCData bcData;
  
  /**/
  if(bc_type_==DIRICHLET){
    for( std::map<unsigned int, unsigned int>::iterator it = cu_ldofs_IDtoIDX_.begin(); it!=cu_ldofs_IDtoIDX_.end(); ++it ){
        if(it->second >= cunknown_idx_cnt) cunknown_idx_cnt=it->second +1;
    }
  }/**/
  
  for(; ldof_it !=  this->getLDofs().end(); ++ldof_it, ++ldof_cnt)
  {
    const gdof_type *const gdof = (*ldof_it)->getGDof();
    unsigned int num_support = gdof->getNumberSupportElements();
    unsigned int support_cnt = 0;
    known = false;
    for(;support_cnt<num_support;++support_cnt) {
      const superelement_type *const se = gdof->getSuperElement(support_cnt);
      unsigned int input_id = se->getElement()->getInputId();

      std::map<unsigned int, unsigned int>::const_iterator sit=region_id_.find(input_id);
      if(sit==region_id_.end()){ // no specification found use default boundary condition
          bcData=defaultBC;
      }else{
        unsigned int region_id = sit->second;
        std::map<unsigned int, BCData>::const_iterator cit=bc_input_.find(region_id);
        if(cit==bc_input_.end()){ // no specification found use default boundary condition
            bcData=defaultBC;
        }else{
            bcData = bc_input_.find(region_id)->second;
        }
      }
        // if we do a pseudo-2d computation, lock z-displacements
        if (is2d_ && bc_type_==DIRICHLET && (*ldof_it)->getLIDX()==2 ) {
            //printf("%% !!! 2D !!!\n");
            known = true;
            (*ldof_it)->setIDX(homo_idx_cnt++);
            (*ldof_it)->setType(HOMOGEN);
//            (*ldof_it)->setIDX(known_idx_cnt++);
//            (*ldof_it)->setType(KNOWN);
//            known_data_tmp.push_back(0.0);
            break;
        }

      
      if(bcData.bc_type==CRACK){ // this is a crack ldof
          // for cracks, set Neumann ldofs to known with 0.0 traction
          if(bc_type_==NEUMANN){
            (*ldof_it)->setIDX(homo_idx_cnt++);
            (*ldof_it)->setType(HOMOGEN);
//            (*ldof_it)->setIDX(known_idx_cnt++);
//            (*ldof_it)->setType(KNOWN);
//            known_data_tmp.push_back(0.0);
          }
          // set Dirichlet ldofs to unknown
          // except at crack tips (0.0 opening displacement)
          if(bc_type_==DIRICHLET){
			unsigned int nd_input_id = se->getElement()->getNode((*ldof_it)->getGDof()->getReferenceElementIdx(support_cnt))->getInputId();
			//printf("%% diri ldof %d has node input id %d --- ",(*ldof_it)->getID(),nd_input_id);
			if(crack_tip_nodes_.count(nd_input_id)){ // this is a crack tip ldof (0.0 opening displacement)
				(*ldof_it)->setIDX(homo_idx_cnt++);
				(*ldof_it)->setType(CHOMOGEN);
				//(*ldof_it)->setIDX(known_idx_cnt++);
				//(*ldof_it)->setType(KNOWN);
				//known_data_tmp.push_back(0.0);
				//printf("\nct bc IDX %d ID %d nd %d", homo_idx_cnt-1, (*ldof_it)->getID(), nd_input_id);
			}else{ // this is an unknown cod-dof
                /**/
                unsigned int nextIDX;
                if( cu_ldofs_IDtoIDX_.count((*ldof_it)->getID())==0 ){
                    nextIDX=cunknown_idx_cnt++; // post-inc!
                    cu_ldofs_IDtoIDX_[ (*ldof_it)->getID() ] = nextIDX;
                }else
                    nextIDX=cu_ldofs_IDtoIDX_[ (*ldof_it)->getID() ];
                
				(*ldof_it)->setIDX(nextIDX);
				(*ldof_it)->setType(CUNKNOWN);
                /*/
				(*ldof_it)->setIDX(cunknown_idx_cnt++);
				(*ldof_it)->setType(CUNKNOWN);
                /**/
			}
          }
//          std::cout << "% crack " << bc_type_ << " ldof " << (*ldof_it)->getID()
//                    << " set to type " << (*ldof_it)->getType()
//                    << " because element_id is " << input_id 
//                    << " and region_id is " << region_id <<std::endl;
          known=true; // avoid defaulting to UNKNOWN below
          break;
          
      }else{ // not a crack ldof, this is based on bcsetter.tpl
//        typename bc_input_type::const_iterator ibc = bc_input_.find(input_id);
//        if (ibc==bc_input_.end())
//          continue;
//  //        HYENA_ERROR_MSG("Element BC Data missing");
//
//        bc_vector_type current_bc_vector=ibc->second;
//
//        if (num_ldofs_per_gdof != current_bc_vector.size())
//          HYENA_ERROR_MSG("size inp BC data doesn't match num_ldofs_per_gedof!");
//
//        bc_data_type current_bc = current_bc_vector[(*ldof_it)->getLIDX()];
//
        if (bcData.bc_type == bc_type_) {
          known = true;
          if (bcData.homogenous == false) {
            (*ldof_it)->setIDX(known_idx_cnt++);
            (*ldof_it)->setType(KNOWN);
            known_data_tmp.push_back(bcData.value[(*ldof_it)->getLIDX()]);
          }
          else if (bcData.homogenous == true){
            (*ldof_it)->setIDX(homo_idx_cnt++);
            (*ldof_it)->setType(HOMOGEN);
          }
          break;
        }
//        if (current_bc.bc_type == ROBIN && bc_type_==NEUMANN) {
//          known = true;
//          (*ldof_it)->setIDX(robin_idx_cnt++);
//          (*ldof_it)->setType(RUNKNOWN);
//          robin_data_1_tmp.push_back(current_bc.value[0]);
//          robin_data_2_tmp.push_back(current_bc.value[1]);
//          break;
//        }
//        // Attention: for interface dofs unknown_idx_counter is used, to
//        // assemble both UNKNOWN and IUNKNOWN in just one matrix block
//        if (current_bc.bc_type == INTERFACE) {
//          known = true;
//          (*ldof_it)->setIDX(unknown_idx_cnt++);
//          (*ldof_it)->setType(IUNKNOWN);
//          break;
//        }
      }
    }
    
    if(known == false){
        (*ldof_it)->setIDX(unknown_idx_cnt++);
        (*ldof_it)->setType(UNKNOWN);
    }

	//printf("%% bc_type %d: ldof %d (lidx %d) has type %d",bc_type_, (*ldof_it)->getID(), (*ldof_it)->getLIDX(), (*ldof_it)->getType());
	//if(known) printf(" and value %lf\n", *(known_data_tmp.end()));
	//else printf("\n");
  }

  known_data.resize(known_data_tmp.size());
  for (unsigned int cnt = 0; cnt < known_data_tmp.size();++cnt) {
    known_data(cnt)= known_data_tmp[cnt];
  }
//  robin_data_1.resize(robin_data_1_tmp.size());
//  for (unsigned int cnt = 0; cnt < robin_data_1_tmp.size();++cnt) {
//    robin_data_1(cnt)= robin_data_1_tmp[cnt];
//  }
//  robin_data_2.resize(robin_data_2_tmp.size());
//  for (unsigned int cnt = 0; cnt < robin_data_2_tmp.size();++cnt) {
//    robin_data_2(cnt)= robin_data_2_tmp[cnt];
//  }
}

#endif //CRACKBCSETTER_TPL
