#include "HyENAWrapperImpl.h"
#include "hyena/mod/crackBCSetter.hpp"
#include "hyena/mod/innerevaluation.hpp"

#include "hyena/core/tria/mesh.hpp"
#include "hyena/core/tria/surfacecurlelement.hpp"
#include "hyena/core/dofs/dofhandler.hpp"
#include "hyena/core/dofs/extractor.hpp"
#include "hyena/core/galerkin/galerkinkernel.hpp"
#include "hyena/core/galerkin/galerkinintegrator3d.hpp"
#include "hyena/core/galerkin/singleintegrator3d.hpp"
#include "hyena/core/galerkin/galerkinassembler.hpp"
#include "hyena/core/galerkin/massassembler.hpp"
#include "hyena/core/fundsol/fundamentalsolutions.hpp"
#include "hyena/core/utilities/mergedata.hpp"


namespace FractureSim{
	HyENAWrapper* HyENAWrapper::constructHyENAWrapper(
		node_map& nodes_, elem_map& elems_, id_map& regions_, id_set& cracks_,
		double youngsModulus, double poissonsRatio, bool is2D
	){
		return new HyENAWrapperImpl(nodes_,elems_,regions_,cracks_,youngsModulus,poissonsRatio,is2D);
	}


    const unsigned int iorder = 4;
    using namespace std;
    using namespace hyena;

    // HyENA typedefs
    typedef Elastostatic3D fsol_type;
    typedef TriaTraits<HyENAWrapperImpl::E, HyENAWrapperImpl::G>::mesh_type Mesh;
    typedef DofTraits<HyENAWrapperImpl::E, HyENAWrapperImpl::G, HyENAWrapperImpl::DA, HyENAWrapperImpl::P, HyENAWrapperImpl::DS> diri_traits_type;
    typedef DofTraits<HyENAWrapperImpl::E, HyENAWrapperImpl::G, HyENAWrapperImpl::NA, HyENAWrapperImpl::P, HyENAWrapperImpl::NS> neum_traits_type;
    typedef diri_traits_type::dof_handler_type diri_handler_type;
    typedef neum_traits_type::dof_handler_type neum_handler_type;
    typedef GalerkinKernel<neum_traits_type, neum_traits_type, fsol_type,
                           SLP_tag> slp_kernel_type;
    typedef GalerkinKernel<neum_traits_type, diri_traits_type, fsol_type, 
                           DLP_tag> dlp_kernel_type;
    typedef GalerkinKernel<diri_traits_type, neum_traits_type, fsol_type, 
                           ADLP_tag> adp_kernel_type;
    typedef GalerkinKernel<diri_traits_type, diri_traits_type, fsol_type, 
                           HSO_tag> hso_kernel_type;

    int HyENAWrapperImpl::init(id_set& crack_tip_nodes, bool ignoreCracks){
		//for(map<unsigned int, BCData>::iterator it = boundaryConditions.begin(); it != boundaryConditions.end(); ++it)
		//	cout << "% boundary " << it->first << " is " << it->second.bc_type << endl;
		//cout << "% crack tip nodes: ";
		//for( set<unsigned int>::iterator it = crack_tip_nodes.begin(); it != crack_tip_nodes.end(); ++it )
		//	cout << *it << " ";
		//cout << endl;
		//return -1;

        cu_IDtoIDX.clear();
        
        fsol_type fsol(youngsModulus,poissonsRatio);
        
		if(mesh!=NULL) delete mesh;
		mesh = new Mesh(nodes,elems);
		if(uh!=NULL) delete uh;
		uh = new diri_handler_type(*mesh);
        if(th!=NULL) delete th;
		th = new neum_handler_type(*mesh);

        CrackBCSetter<E,G,DA,P,DS> diri_set_bc(*uh, boundaryConditions, DIRICHLET, regions, crack_tip_nodes, cu_IDtoIDX, is2D);
        diri_set_bc(diri_data);

        CrackBCSetter<E,G,NA,P,NS> neum_set_bc(*th, boundaryConditions, NEUMANN, regions, crack_tip_nodes, cu_IDtoIDX);
        neum_set_bc(neum_data);

//        printf("%% cu_IDtoIDX\n");
//        for(id_map::iterator it=cu_IDtoIDX.begin(); it!=cu_IDtoIDX.end(); ++it)
//            printf("%% ID %d --> IDX %d\n",it->first, it->second);
        
        diri_traits_type::extractor_type diri_extractor(*uh);
		dk_ldofs.clear(); du_ldofs.clear(); cu_ldofs.clear();
		diri_extractor(dk_ldofs, KNOWN);
        diri_extractor(du_ldofs, UNKNOWN);
        if(!ignoreCracks) diri_extractor(cu_ldofs, CUNKNOWN);

        neum_traits_type::extractor_type neum_extractor(*th);
		nk_ldofs.clear(); nu_ldofs.clear();
        neum_extractor(nk_ldofs, KNOWN);
        neum_extractor(nu_ldofs, UNKNOWN);

		M = du_ldofs.size();
        N = nu_ldofs.size();
		L = cu_ldofs.size();

        GalerkinAssembler<slp_kernel_type> assemble_slp(fsol, iorder);
        GalerkinAssembler<dlp_kernel_type> assemble_dlp(fsol, iorder);
        GalerkinAssembler<adp_kernel_type> assemble_adp(fsol, iorder);
        GalerkinAssembler<hso_kernel_type> assemble_hso(fsol, iorder);
        MassAssembler<neum_traits_type, diri_traits_type> assemble_mass_dlp(iorder);
        MassAssembler<diri_traits_type, neum_traits_type> assemble_mass_adp(iorder);

        V.resize(N, N); V.setZero();
        K.resize(N, M); K.setZero();
		Kc.resize(N,L);Kc.setZero();
        D.resize(M, M); D.setZero();
		Dc.resize(M,L);Dc.setZero();
		Dcc.resize(L,L); Dcc.setZero();

        fD.resize(N); fD.setZero(); test1.resize(N); test1.setZero();
        fN.resize(M); fN.setZero(); test2.resize(M); test2.setZero();
		// we only treat traction-free cracks, so there is no rhs vector for crack dofs
        
        assemble_slp(nu_ldofs, V);
        assemble_dlp(nu_ldofs, du_ldofs, K);
        assemble_hso(du_ldofs, D);

		assemble_slp(nu_ldofs, nk_ldofs, neum_data, -1, fD);// fD += -V * gN
        assemble_dlp(nu_ldofs, dk_ldofs, diri_data, 1, fD);// fD += (1/2I + K) * gD
        assemble_mass_dlp(nu_ldofs, dk_ldofs, diri_data, 0.5, fD);

		assemble_adp(du_ldofs, nk_ldofs, neum_data, -1, fN);// fN += (1/2I - K') * gN
        assemble_mass_adp(du_ldofs, nk_ldofs, neum_data, 0.5, fN);
        assemble_hso(du_ldofs, dk_ldofs, diri_data, -1, fN);// fN += -D * gD

		assemble_dlp(nu_ldofs, cu_ldofs, Kc);
		assemble_hso(du_ldofs, cu_ldofs, Dc);
		assemble_hso(cu_ldofs, cu_ldofs, Dcc); // don't use symmetric assembly for the crack dofs

		 //cout << endl;
		 //cout << "V  =[" << endl << V   << "];" << endl;
		 //cout << "K  =[" << endl << K   << "];" << endl;
		 //cout << "D  =[" << endl << D   << "];" << endl;
		 //cout << "fD =[" << endl << fD  << "];" << endl;
		 //cout << "fN =[" << endl << fN  << "];" << endl;
		 //cout << "Kc =[" << endl << Kc  << "];" << endl;
		 //cout << "Dc =[" << endl << Dc  << "];" << endl;
		 //cout << "Dcc=[" << endl << Dcc << "];" << endl;
		 //cout << endl << "%";

		if(N>0){
			// invert A=[V -K; K' D] into A^-1 = [X -Y; Y' Z] and
			// store X-->V, Y-->K, and Z-->D
			// Dinv = inv(D);
			// Vinv = inv(V);
			// V = inv(V+K*Dinv*K');
			// D = inv(D+K'*Vinv*K);
			// K = -V*K*Dinv;
			// Eigen::ColPivHouseholderQR<matrix_type> qrD, qrV;
			// qrD.compute(D);
			// qrV.compute(V);
			//if(qrV.rank()<N) printf("\n!!! V might be rank deficient\n");
			//if(qrD.rank()<M) printf("\n!!! D might be rank deficient\n");
			// the partial piv. LU seems to be a bit faster than the QR, and these matrices should all be nice enough to solve
			luV.compute(V);
			V = (V+K*D.partialPivLu().inverse()*K.transpose()).partialPivLu().inverse();
			D = (D+K.transpose()*luV.inverse()*K).partialPivLu().inverse();
			K = -(D*K.transpose()*luV.inverse()).transpose();
			// cout << endl;
			// cout << "Vi  =[" << endl << V   << "];" << endl;
			// cout << "Ki  =[" << endl << K   << "];" << endl;
			// cout << "Di  =[" << endl << D   << "];" << endl;
			// cout << endl << "%";
			
			// // newer 2-step Schur complement pre-factorization
			// qrV.compute(V);
			// qrS.compute(K.transpose()*qrV.inverse()*K+D); // Schur complement of the first step K'*V^-1*K + D

		}else{
			printf("\n%% ... regularizing Neumann problem");
			// build a regularizer T for D, compute D'*D + T'*T and D'*fN
			matrix_type T(6,M);
			T.setZero(); double area; unsigned int i,j,k;
			for(elem_map::iterator el=elems.begin(); el!=elems.end(); ++el)
			if( cracks.count( regions[el->first] ) == 0 ){
				// every node has 3 displacement-dofs and gets a 3x3 block in T where we measure rigid translations
				// as well as a 3x3 block (in rows 4-6) where we measure rigid rotation
				// we treat rotations around the coordinate origin
				// ideally the mesh should already be specified in such a way that this makes sense (i.e. center of mass close to [0,0,0])
				Eigen::Vector3d r,
					a (nodes[el->second[0]][0], nodes[el->second[0]][1], nodes[el->second[0]][2]),
					b (nodes[el->second[1]][0], nodes[el->second[1]][1], nodes[el->second[1]][2]),
					c (nodes[el->second[2]][0], nodes[el->second[2]][1], nodes[el->second[2]][2]);
				area = 0.5*(a-c).cross(b-c).norm();
				i = 3*(el->second[0] - NODE_BASE_INDEX); // first dof of node a (NOTE THAT THIS ASSUMES ALL SURFACE ELEMENTS HAVE LOWER IDs THAN ALL FRACTURE ELEMENTS!)
				j = 3*(el->second[1] - NODE_BASE_INDEX); // first dof of node b
				k = 3*(el->second[2] - NODE_BASE_INDEX); // first dof of node c
				// the translation block contribution from each element to a node is 1/3*area*I
				T( 0, i ) += area/3.0; T( 1, 1+i ) += area/3.0; T( 2, 2+i ) += area/3.0;
				T( 0, j ) += area/3.0; T( 1, 1+j ) += area/3.0; T( 2, 2+j ) += area/3.0;
				T( 0, k ) += area/3.0; T( 1, 1+k ) += area/3.0; T( 2, 2+k ) += area/3.0;
				// the rotational block contribution from each element to a node is 1/12*area*M
				// where M is the matrix representing a cross product from the right with r
				// and   r is 2*a+b+c for node a; a+2*b+c for node b etc.
				// M = [ 0 rz -ry ; -rz 0 rx ; ry -rx 0 ]
				r = 2.0*a+b+c;
				T( 3, 1+i ) += area/12.0 * r[2]; T( 3, 2+i ) -= area/12.0 * r[1];
				T( 4,   i ) -= area/12.0 * r[2]; T( 4, 2+i ) += area/12.0 * r[0];
				T( 5,   i ) += area/12.0 * r[1]; T( 5, 1+i ) -= area/12.0 * r[0];
				r = a+2.0*b+c;
				T( 3, 1+j ) += area/12.0 * r[2]; T( 3, 2+j ) -= area/12.0 * r[1];
				T( 4,   j ) -= area/12.0 * r[2]; T( 4, 2+j ) += area/12.0 * r[0];
				T( 5,   j ) += area/12.0 * r[1]; T( 5, 1+j ) -= area/12.0 * r[0];
				r = a+b+2.0*c;
				T( 3, 1+k ) += area/12.0 * r[2]; T( 3, 2+k ) -= area/12.0 * r[1];
				T( 4,   k ) -= area/12.0 * r[2]; T( 4, 2+k ) += area/12.0 * r[0];
				T( 5,   k ) += area/12.0 * r[1]; T( 5, 1+k ) -= area/12.0 * r[0];
			}
			T /= totalArea;
			T *= (D.trace()/M); // scale for the regularizer is somewhat arbitrary, it should be about as "important" as D itself
			fN = D*fN; // .transpose() not coded since D is symmetric
			luS.compute(D*D + T.transpose()*T);
			//cout << "T =[" << endl << T << "];" << endl;
			//cout << "DD =[" << endl << D.transpose()*D << "];" << endl;
			//cout << "DfN=[" << endl << fN << "];" << endl;
		}
		return 0;
    }

    int HyENAWrapperImpl::addCrack(
		std::set<unsigned int>& crack_tip_nodes
	){
		// update matrices for next solve ...
        fsol_type fsol(youngsModulus,poissonsRatio);
        
		if(mesh!=NULL) delete mesh;
        mesh = new Mesh(nodes,elems);
		if(uh!=NULL) delete uh;
		uh = new diri_handler_type(*mesh);
        if(th!=NULL) delete th;
		th = new neum_handler_type(*mesh);

        CrackBCSetter<E,G,DA,P,DS> diri_set_bc(*uh, boundaryConditions, DIRICHLET, regions, crack_tip_nodes, cu_IDtoIDX, is2D);
        diri_set_bc(diri_data);

        CrackBCSetter<E,G,NA,P,NS> neum_set_bc(*th, boundaryConditions, NEUMANN, regions, crack_tip_nodes, cu_IDtoIDX);
        neum_set_bc(neum_data);

//        printf("%% cu_IDtoIDX\n");
//        for(id_map::iterator it=cu_IDtoIDX.begin(); it!=cu_IDtoIDX.end(); ++it)
//            printf("%% ID %d --> IDX %d\n",it->first, it->second);

        diri_traits_type::extractor_type diri_extractor(*uh);

		int n_old_cu_ldofs = L;
		dk_ldofs.clear(); du_ldofs.clear();
		cu_ldofs.clear();
		diri_extractor(dk_ldofs, KNOWN);
        diri_extractor(du_ldofs, UNKNOWN);
        diri_extractor(cu_ldofs, CUNKNOWN);

        neum_traits_type::extractor_type neum_extractor(*th);
		nk_ldofs.clear(); nu_ldofs.clear();
        neum_extractor(nk_ldofs, KNOWN);
        neum_extractor(nu_ldofs, UNKNOWN);

        GalerkinAssembler<dlp_kernel_type> assemble_dlp(fsol, iorder);
        GalerkinAssembler<hso_kernel_type> assemble_hso(fsol, iorder);

		M = du_ldofs.size();
		N = nu_ldofs.size();
		L = cu_ldofs.size();
		// split cu_ldofs in old and new ones
        std::vector<const LDof<E,G,DA,P,DS>*> new_cu_ldofs, old_cu_ldofs;
        // all cu_ldofs with IDX lower than n_old_cu_ldofs are old, others are new
        for(int i=0; i<cu_ldofs.size(); ++i){
            if(cu_ldofs[i]->getIDX() < n_old_cu_ldofs)
                old_cu_ldofs.insert(old_cu_ldofs.end(), cu_ldofs[i]);
            else
                new_cu_ldofs.insert(new_cu_ldofs.end(), cu_ldofs[i]);
        }

//        printf("NEW cu_ldofs\n");
//        for(int i=0; i<new_cu_ldofs.size(); ++i){
//            printf("ID %d --> IDX %d **\n", new_cu_ldofs[i]->getID(), new_cu_ldofs[i]->getIDX());
//        }
//        printf("cu_ldofs\n");
//        for(int i=0; i<cu_ldofs.size(); ++i){
//            printf("ID %d --> IDX %d **\n", cu_ldofs[i]->getID(), cu_ldofs[i]->getIDX());
//        }      
        

		Kc.conservativeResize(N,L);
		Dc.conservativeResize(M,L);
		Dcc.conservativeResize(L,L);

		Kc.block(0,n_old_cu_ldofs, N,L-n_old_cu_ldofs).setZero();
		Dc.block(0,n_old_cu_ldofs,M,L-n_old_cu_ldofs).setZero();
		Dcc.block(0,n_old_cu_ldofs, n_old_cu_ldofs,L-n_old_cu_ldofs).setZero();
		Dcc.block(n_old_cu_ldofs,0, L-n_old_cu_ldofs,n_old_cu_ldofs).setZero();
		Dcc.block(n_old_cu_ldofs,n_old_cu_ldofs,L-n_old_cu_ldofs,L-n_old_cu_ldofs).setZero();

		assemble_dlp(nu_ldofs, new_cu_ldofs, Kc,0,0,0,n_old_cu_ldofs);
		assemble_hso(du_ldofs, new_cu_ldofs, Dc,0,0,0,n_old_cu_ldofs);
		assemble_hso(old_cu_ldofs, new_cu_ldofs, Dcc,0,0,0,n_old_cu_ldofs,n_old_cu_ldofs);
        assemble_hso(new_cu_ldofs,new_cu_ldofs, Dcc,0,0,n_old_cu_ldofs,n_old_cu_ldofs);
        // assembly runs over old crack ldofs per rows and new crack ldofs per columns
        // results are written into matrix D, with indices shifted by (M,M)
        // only entries to the left and below of (M,M+n_old_cu_ldofs) (incl.) are written
        // only entries above M+n_old_cu_ldofs (excl.) are written
        // this ensures that only the block corresponding to Dcoc is modified in D.

		Dcc.block(n_old_cu_ldofs,0, L-n_old_cu_ldofs,n_old_cu_ldofs) =
			Dcc.block(0,n_old_cu_ldofs, n_old_cu_ldofs,L-n_old_cu_ldofs).transpose().eval();

		// cout << endl;
		// cout << "Kc =[" << endl << Kc  << "];" << endl;
		// cout << "Dc =[" << endl << Dc  << "];" << endl;
		// cout << "Dcc=[" << endl << Dcc << "];" << endl;
		// cout << endl << "%";

		return 0;
    }

	int HyENAWrapperImpl::updateRHS(vect3d_map& newBcData, id_set& crack_tip_nodes){
		for(vect3d_map::iterator it=newBcData.begin(); it!=newBcData.end(); ++it){
			if( boundaryConditions.count( it->first ) ){
				boundaryConditions[it->first].value.assign(3,0.0);
				boundaryConditions[it->first].value[0] = it->second[0];
				boundaryConditions[it->first].value[1] = it->second[1];
				boundaryConditions[it->first].value[2] = it->second[2];
			}
		}
		
		cu_IDtoIDX.clear();
        fsol_type fsol(youngsModulus,poissonsRatio);
//		if(mesh!=NULL) delete mesh;
//		mesh = new Mesh(nodes,elems);
//		if(uh!=NULL) delete uh;
//		uh = new diri_handler_type(*mesh);
//		if(th!=NULL) delete th;
//		th = new neum_handler_type(*mesh);

        CrackBCSetter<E,G,DA,P,DS> diri_set_bc(*uh, boundaryConditions, DIRICHLET, regions, crack_tip_nodes, cu_IDtoIDX, is2D);
        diri_set_bc(diri_data);

        CrackBCSetter<E,G,NA,P,NS> neum_set_bc(*th, boundaryConditions, NEUMANN, regions, crack_tip_nodes, cu_IDtoIDX);
        neum_set_bc(neum_data);
    
//		diri_traits_type::extractor_type diri_extractor(*uh);
//		dk_ldofs.clear(); du_ldofs.clear(); cu_ldofs.clear();
//		diri_extractor(dk_ldofs, KNOWN);
//		diri_extractor(du_ldofs, UNKNOWN);
//		diri_extractor(cu_ldofs, CUNKNOWN);
//
//		neum_traits_type::extractor_type neum_extractor(*th);
//		nk_ldofs.clear(); nu_ldofs.clear();
//		neum_extractor(nk_ldofs, KNOWN);
//		neum_extractor(nu_ldofs, UNKNOWN);
//
//		M = du_ldofs.size();
//		N = nu_ldofs.size();
//		L = cu_ldofs.size();

        GalerkinAssembler<slp_kernel_type> assemble_slp(fsol, iorder);
        GalerkinAssembler<dlp_kernel_type> assemble_dlp(fsol, iorder);
        GalerkinAssembler<adp_kernel_type> assemble_adp(fsol, iorder);
        GalerkinAssembler<hso_kernel_type> assemble_hso(fsol, iorder);
        MassAssembler<neum_traits_type, diri_traits_type> assemble_mass_dlp(iorder);
        MassAssembler<diri_traits_type, neum_traits_type> assemble_mass_adp(iorder);

        fD.resize(N); fD.setZero();
        fN.resize(M); fN.setZero();

		assemble_slp(nu_ldofs, nk_ldofs, neum_data, -1, fD);// fD += -V * gN
        assemble_dlp(nu_ldofs, dk_ldofs, diri_data, 1, fD);// fD += (1/2I + K) * gD
        assemble_mass_dlp(nu_ldofs, dk_ldofs, diri_data, 0.5, fD);

		assemble_adp(du_ldofs, nk_ldofs, neum_data, -1, fN);// fN += (1/2I - K') * gN
        assemble_mass_adp(du_ldofs, nk_ldofs, neum_data, 0.5, fN);
        assemble_hso(du_ldofs, dk_ldofs, diri_data, -1, fN);// fN += -D * gD

		if( N==0 ){ // we're solving a regularized Neumann problem (D'*D+T'*T)u=D'*fN
			fN = D*fN; // .transpose() not coded since D is symmetric
		}
		return 0;
    }

    int HyENAWrapperImpl::compute(){
		// preserve matrices while solving
		vector_type u(M+L), t(N);
		if( N>0 ){ // old version for mixed problems
			// // here we solve for displacement and COD in one system
			// // D_complete = [ D , Dc ; Dc' , Dcc ]
			// matrix_type Dcopy(M+L,M+L); //Dcopy.setZero();
			// matrix_type Kcopy(N,M+L);
			// vector_type fNcopy(M+L); //fNcopy.setZero();
			// // printf("\n\n old solve (N,M,L)=(%d,%d,%d)\n\n", N,M,L);
			// Dcopy.block(0,0,M,M) = D;
			// Dcopy.block(0,M,M,L) = Dc;
			// Dcopy.block(M,0,L,M) = Dc.transpose();
			// Dcopy.block(M,M,L,L) = Dcc;
			// Kcopy.block(0,0,N,M) = K;
			// Kcopy.block(0,M,N,L) = Kc;
			// fNcopy.block(0,0,M,1)= fN;
			// fNcopy.block(M,0,L,1).setZero();
			// // c = V^-1 * fD
			// vector_type c;
			// Eigen::PartialPivLU<matrix_type> lu(V);
			// c = lu.solve(fD);
			// // fN -= K' * c
			// fNcopy.noalias() -= Kcopy.transpose() * c;
			// // H = V^-1 * K;
			// matrix_type H;
			// H.noalias() = lu.inverse()*Kcopy;
			// // D += K' * H
			// Dcopy.noalias() += Kcopy.transpose() * H;
			// // u = D^-1 * fN
			// u = Dcopy.partialPivLu().solve(fNcopy);
			// // t = c + H * u
			// t.noalias() = c + H * u;
			
			// new solver: in init() we've already inverted the matrix blocks
			// corresponding to the non-crack DOFs and stored them in V, K, and D
			// such that A^1 = [V -K; K' D]
			matrix_type H1,H2,Dcc_copy=Dcc;
			vector_type cu(L),c1,c2;
			c1 = V*fD - K*fN;
			c2 = K.transpose()*fD + D*fN;
			if(L>0){
				H1 = -(V*Kc + K*Dc);
				H2 = D*Dc - K.transpose()*Kc;
				Dcc_copy -= (Kc.transpose()*H1 + Dc.transpose()*H2 );
				cu = Dcc_copy.ldlt().solve( Kc.transpose()*c1 + Dc.transpose()*c2 ); // this is actually -cu
				t  = c1 + H1*cu;
				u.block(0,0,M,1) = c2 + H2*cu;
				u.block(M,0,L,1) = -cu; // now we flip the sign back
			}else{
				t  = c1;
				u  = c2;
			}

		}else{ // solving the Neumann problem - needs to be regularized by init()
			if(L>0){ // there are fractures ...
				vector_type cu(L), uu = luS.solve(fN); //displacements in unfractured case
				matrix_type H;
				H.noalias() = luS.inverse()*D*Dc; // D.transpose not coded since D is symmetric
				cu = (Dc.transpose()*H - Dcc)
					.ldlt().solve( Dc.transpose()*uu ); // COD
				u.block(0,0,M,1) = uu - H*cu;
				u.block(M,0,L,1) = cu;
			}else{
				u = luS.solve(fN);
			}
		}

		//cout << endl;
		//cout << "uraw  =[" << endl << u   << "];" << endl;
		//cout << "qraw  =[" << endl << t   << "];" << endl;
		//cout << endl << "%";

        // post-processing (note that the ldof objects are managed by the DofHandler and will be deleted by the handler's destructor)

		displacements.resize(3*nodes.size());
		tractions.resize(3*elems.size());

		displacements.setZero();
		tractions.setZero();

		//printf("\ndk array size %d, u array size %d", diri_data.rows(), u.rows());
		MergeData merge_data;
		//printf("\nmerging dk");
		merge_data(dk_ldofs, diri_data, displacements);
		//printf("\nmerging du");
		merge_data(du_ldofs, u, displacements);
		//printf("\n%%merging cu\ncu_map=[ %% ID, node, IDX");
		merge_data(cu_ldofs, u, displacements, du_ldofs.size(), NODE_BASE_INDEX);
		//printf("]; %%");
		//printf("\nmerging neum");
		merge_data(nk_ldofs, neum_data, tractions);
		merge_data(nu_ldofs, t, tractions);

		return 0;
    }

	int HyENAWrapperImpl::buildBoundaryCnd(const vector<string>& specs){
		BCData bc;
		boundaryConditions.clear();
        cracks.clear();
		for(int i=0; i<specs.size(); ++i){ // parse and add ...
			int bnd_id;
			bnd_id=parseBoundaryCnd(specs[i].c_str(), bc);
			if(bnd_id>=0) boundaryConditions[bnd_id]=bc;
		}
		return boundaryConditions.size();
	}

    int HyENAWrapperImpl::parseBoundaryCnd(const char* spec, BCData& bc){
        //printf("parsing bnd-cnd %s\n", spec);
        bc.homogenous=false;
		bc.value.clear();
		bc.bc_type=NO_BC_TYPE;

        double fx,fy,fz;
        unsigned long bnd_id;
        int tokens; char check;
        tokens=sscanf(spec, "%lu(%lf,%lf,%lf)", &bnd_id, &fx, &fy, &fz);
        if(tokens==4){ // force bnd.cnd.
            //printf("%% ... boundary %d has external traction (%.3le, %.3le, %.3le)\n", bnd_id, fx,fy,fz);
            //...
            bc.bc_type=NEUMANN;
            bc.value.push_back(fx);
            bc.value.push_back(fy);
            bc.value.push_back(fz);
        }else{//not correct format for force bnd.cnd.
            tokens=sscanf(spec, "%lu[%lf,%lf,%lf]", &bnd_id, &fx, &fy, &fz);
            if(tokens==4){ // displacement bnd.cnd.
                //printf("%% ... boundary %d has external traction (%.3le, %.3le, %.3le)\n", bnd_id, fx,fy,fz);
                //...
                bc.bc_type=DIRICHLET;
                bc.value.push_back(fx);
                bc.value.push_back(fy);
                bc.value.push_back(fz);
            }else{
                check=-1;
                tokens=sscanf(spec, "%lu(fixed%c", &bnd_id, &check);
                //printf(" read %d tokens, check is %d\n",tokens, check);
                if(tokens==2 && check==')'){ //proper format of fixed boundary, add nonDOFs
                    //printf("%% ... boundary %d is fixed\n", bnd_id);
                    //...
                    bc.bc_type=DIRICHLET;
                    bc.homogenous=true;
                    //bc.value.push_back(0.0);
                    //bc.value.push_back(0.0);
                    //bc.value.push_back(0.0);
                }else{
                    check=-1;
                    tokens=sscanf(spec, "%lu(crack%c", &bnd_id, &check);
                    if(tokens==2 && check==')'){ //proper format of crack boundary, add nonDOFs
                        bc.bc_type=CRACK;
                        cracks.insert(bnd_id);
                        //printf("%% ... boundary %d is a crack\n", bnd_id);
                    }else{
                        fprintf(stderr, "invalid boundary condition format: %s -- ignored\n",spec);
                        return -1;
                    }
                }
            }
        }
        return bnd_id;
    }

    int HyENAWrapperImpl::addCrackBoundaryCnd(unsigned int bnd_id){
        if(boundaryConditions.count(bnd_id)!=0) return -1; // already have this boundary
        BCData bc;
        bc.homogenous=false; // does not apply to crack surfaces
		bc.value.clear(); // does not apply to crack surfaces (traction-free by default anyway)
        bc.bc_type=CRACK;
        cracks.insert(bnd_id);
        boundaryConditions[bnd_id]=bc;
        return 0;
    }
    
    void HyENAWrapperImpl::innerEval(vector_type& result, const vector_type& coords){
        vector<Point3> inner_pts;
		vector<const LDof<E,G,DA,P,DS>*> da_ldofs, dh_ldofs;
		unsigned int n_pts=coords.size()/3;

		fsol_type fsol(youngsModulus,poissonsRatio);
		InnerEvaluation<diri_traits_type,diri_traits_type,fsol_type,DLP_tag> inner_evaluation_dlp(fsol);
		vector_type tmp1, tmp2;

		inner_pts.assign(n_pts, Point3(0.0,0.0,0.0) );
        for(unsigned int i=0; i<n_pts; ++i){
            inner_pts[i][0]=coords[3*i  ];
            inner_pts[i][1]=coords[3*i+1];
            inner_pts[i][2]=coords[3*i+2];
		}

		diri_traits_type::extractor_type diri_extractor(*uh);
		diri_extractor(dh_ldofs, HOMOGEN);

		da_ldofs.insert(da_ldofs.end(),dk_ldofs.begin(),dk_ldofs.end());
		da_ldofs.insert(da_ldofs.end(),du_ldofs.begin(),du_ldofs.end());
		da_ldofs.insert(da_ldofs.end(),dh_ldofs.begin(),dh_ldofs.end());

        tmp1.resize( 3* n_pts ); tmp1.setZero();

        inner_evaluation_dlp(da_ldofs,displacements,inner_pts,-1.0,tmp1);

		result.resize(coords.size()); result.setZero();
        for(unsigned int i=0; i<n_pts; ++i){
            result[3*i  ]=tmp1[3*i  ];
            result[3*i+1]=tmp1[3*i+1];
            result[3*i+2]=tmp1[3*i+2];
		}
    }

	vector_type& HyENAWrapperImpl::computeCrackBaseDisplacements(){
		//have data:
		//vector_type displacements, tractions;
		//std::vector<const LDof<E,G,DA,P,DS>*> dk_ldofs, du_ldofs, cu_ldofs;
		//std::vector<const LDof<E,G,NA,P,NS>*> nk_ldofs, nu_ldofs;
		//write data: crackInnerEvals
		// - run an inner evaluation for each crack over all surfaces except the crack itself
		crackInnerEvals.resize(3*nodes.size());
		crackInnerEvals.setZero();

		map<unsigned int, vector<Point3> >       inner_pts;
		map<unsigned int, vector<unsigned int> > inner_nds;
		id_set done; // store which nodes are already done

//		map<unsigned int, vector<const LDof<E,G,DA,P,DS>*> > cu_ldofs_region;
		vector<const LDof<E,G,DA,P,DS>*> da_ldofs, dh_ldofs;
//		vector<const LDof<E,G,NA,P,NS>*> na_ldofs;
		unsigned int r,k;

		fsol_type fsol(youngsModulus,poissonsRatio);
		InnerEvaluation<diri_traits_type,diri_traits_type,fsol_type,DLP_tag> inner_evaluation_dlp(fsol);
//		InnerEvaluation<neum_traits_type,neum_traits_type,fsol_type,SLP_tag> inner_evaluation_slp(fsol);
		vector_type tmp1, tmp2;

		// build vectors of inner evaluations points
		for(elem_map::iterator it=elems.begin(); it!=elems.end(); ++it){ // for each element
			if( cracks.count( r=regions[it->first] ) !=0 ){ // which is a crack-element
				for(k=0; k < it->second.size(); ++k){ // check all nodes
					if( done.count( it->second[k] ) ==0 ){ // if a node is not done yet
						inner_pts[r].push_back( Point3( // insert it's coordinates into inner_pts
							nodes[it->second[k]][0], nodes[it->second[k]][1], nodes[it->second[k]][2]
						));
						inner_nds[r].push_back( it->second[k] ); // and it's ID into inner_nds
						done.insert( it->second[k] ); // now it's done
					}
				}
			}
		}

//		// build ldof vectors
//		for(k=0; k < cu_ldofs.size(); ++k){
//			r = regions[ cu_ldofs[k]->getGDof()->getSuperElement(0)->getElement()->getInputId() ];
//			cu_ldofs_region[r].push_back(cu_ldofs[k]);
//		}

		diri_traits_type::extractor_type diri_extractor(*uh);
		diri_extractor(dh_ldofs, HOMOGEN);
		da_ldofs.insert(da_ldofs.end(),dk_ldofs.begin(),dk_ldofs.end());
		da_ldofs.insert(da_ldofs.end(),du_ldofs.begin(),du_ldofs.end());
		da_ldofs.insert(da_ldofs.end(),dh_ldofs.begin(),dh_ldofs.end());
//		na_ldofs.insert(na_ldofs.end(),nk_ldofs.begin(),nk_ldofs.end());
//		na_ldofs.insert(na_ldofs.end(),nu_ldofs.begin(),nu_ldofs.end());

		for(id_set::iterator crack = cracks.begin(); crack!=cracks.end(); ++crack){
			tmp1.resize( 3* inner_pts[*crack].size() );
			tmp2.resize( 3* inner_pts[*crack].size() );
			tmp1.setZero(); tmp2.setZero();
		
			inner_evaluation_dlp(da_ldofs,displacements,inner_pts[*crack],-1.0,tmp1);
			//inner_evaluation_slp(na_ldofs,tractions    ,inner_pts[*crack], 1.0,tmp2);
			//tmp1+=tmp2;

//			for(id_set::iterator c_it = cracks.begin(); c_it!=cracks.end(); ++c_it){
//				if((*c_it)!=(*crack)){
//					tmp2.setZero();
//					inner_evaluation_dlp(cu_ldofs_region[*c_it],displacements,inner_pts[*crack],-0.5,tmp2);
//					tmp1+=tmp2;
//				}
//			}

			for(k=0; k < inner_nds[*crack].size(); ++k){
				unsigned int nd = inner_nds[*crack][k]-NODE_BASE_INDEX;
				crackInnerEvals[3*nd  ] = tmp1[3*k  ];
				crackInnerEvals[3*nd+1] = tmp1[3*k+1];
				crackInnerEvals[3*nd+2] = tmp1[3*k+2];
			}
			
		}

		return crackInnerEvals;
	}

	int HyENAWrapperImpl::updateRHSwithCOD(vect3d_map& newCOD, elem_map& support_elems, id_set& support_crack_tip_nodes){
		// given (estimated) crack opening displacements
		// at nodes that are part of support_elems
		// update the RHS (and the solution vectors) as
		// fN += Kc*u and fD -= Dc*u where u is the new COD
		
        fsol_type fsol(youngsModulus,poissonsRatio);
		//build a temporary mesh including only support_elems
		hyena::Mesh<E,G> tmpmesh(nodes,support_elems);
		diri_handler_type tmpuh(tmpmesh);
		neum_handler_type tmpth(tmpmesh);

		//set BCs on the temporary mesh
		id_map tmpcu_IDtoIDX; vector_type tmpdiri_data, tmpneum_data;
        CrackBCSetter<E,G,DA,P,DS> diri_set_bc(tmpuh, boundaryConditions, DIRICHLET, regions, support_crack_tip_nodes, tmpcu_IDtoIDX);
        diri_set_bc(tmpdiri_data);

        CrackBCSetter<E,G,NA,P,NS> neum_set_bc(tmpth, boundaryConditions, NEUMANN, regions, support_crack_tip_nodes, tmpcu_IDtoIDX);
        neum_set_bc(tmpneum_data);
    
		//only work with the cu (now known) dofs in the tmp mesh
		diri_traits_type::extractor_type diri_extractor(tmpuh);
		std::vector<const LDof<E,G,DA,P,DS>*> tmpcu_ldofs;
		diri_extractor(tmpcu_ldofs, CUNKNOWN);
		unsigned int tmpL = tmpcu_ldofs.size();

        GalerkinAssembler<dlp_kernel_type> assemble_dlp(fsol, iorder);
        GalerkinAssembler<hso_kernel_type> assemble_hso(fsol, iorder);
		vector_type fDup, fNup;
		fDup.resize(N); fDup.setZero();
        fNup.resize(M); fNup.setZero();
		tmpdiri_data.resize( 3*newCOD.size() );
		tmpdiri_data.setZero();
		unsigned int k=0;
		for(vect3d_map::iterator it=newCOD.begin(); it!=newCOD.end(); ++it){
			tmpdiri_data[k++] = it->second[0]; // post-inc!
			tmpdiri_data[k++] = it->second[1];
			tmpdiri_data[k++] = it->second[2];
		}

		//printf("\n%% *** RHS update:");
		//printf("\n%% *** nu %d, du %d tcu %d",nu_ldofs.size(),du_ldofs.size(),tmpcu_ldofs.size());
		//printf("\n%% *** data %d \n", tmpdiri_data.size());
		//printf("\n%% *** COD ID to IDX map:");
		//print_map(tmpcu_IDtoIDX);


        assemble_dlp(nu_ldofs, tmpcu_ldofs, tmpdiri_data,  1, fDup);// fD += K*gD
        assemble_hso(du_ldofs, tmpcu_ldofs, tmpdiri_data, -1, fNup);// fN -= D*gD

		//cout << endl << "fDup  =[" << endl << fDup << "];";
		//cout << endl << "fNup  =[" << endl << fNup << "];" << endl;

		//ToDo: figure out how to handle this case properly
		if( N==0 ){ // we're solving a regularized Neumann problem (D'*D+T'*T)u=D'*fN
			fNup = D*fNup; // .transpose() not coded since D is symmetric
		}
		fD += fDup; fN+=fNup;
		//test1=0.6*test1+fDup; test2=0.6*test2+fNup;
		//test1=0.4142*test1+fDup; test2=0.4142*test2+fNup; // const factor is 0.4142 = sqrt(2)-1
		//fD += test1; fN+=test2;
		return compute();
		//return 0;
    }

	
	HyENAWrapperImpl::~HyENAWrapperImpl(){
		if(mesh!=NULL) delete mesh;
		if(uh!=NULL) delete uh;
		if(th!=NULL) delete th;
	}




	
	int HyENAWrapperImpl::innerEvalTest(std::string fileName){
		/**
		//have data:
		//vector_type displacements, tractions;
		//std::vector<const LDof<E,G,DA,P,DS>*> dk_ldofs, du_ldofs, cu_ldofs;
		//std::vector<const LDof<E,G,NA,P,NS>*> nk_ldofs, nu_ldofs;
		typedef TriaTraits<E, G>::global_point_type global_point_type;
		std::vector<global_point_type> pnts;
		unsigned int N_pts=8;
		double l = 1.2;
		double delta = l/(N_pts-1);
		double x,y,z;
		//The Visualization Toolkit supports five different dataset formats: structured points, structured grid, rectilinear 
		//grid, unstructured grid, and polygonal data. Data with implicit topology (structured data such as vtkImageData and 
		//vtkStructuredGrid) are ordered with x increasing fastest, then y, then z. http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.html
		for(        int k=0; k<N_pts;++k) {
			for(    int j=0; j<N_pts;++j) {
				for(int i=0; i<N_pts;++i) {
					x=delta*i-l/2.;
					y=delta*j-l/2.;
					z=delta*k-l/2.;
					pnts.push_back(global_point_type(x,y,z));
				}
			}
		}
		static const unsigned int num_ldofs_per_gdof =
		   ProblemTraits<P>::num_ldofs_per_gdof;
		unsigned int M=N_pts*N_pts*N_pts*num_ldofs_per_gdof;
		vector_type tmp1(M), tmp2(M), u_inner(M);
		tmp1.setZero(); tmp2.setZero(); u_inner.setZero();
		
		std::vector<const LDof<E,G,DA,P,DS>*> da_ldofs;
		da_ldofs.insert(da_ldofs.end(),dk_ldofs.begin(),dk_ldofs.end());
		da_ldofs.insert(da_ldofs.end(),du_ldofs.begin(),du_ldofs.end());
		da_ldofs.insert(da_ldofs.end(),cu_ldofs.begin(),cu_ldofs.end());
		
		std::vector<const LDof<E,G,NA,P,NS>*> na_ldofs;
		na_ldofs.insert(na_ldofs.end(),nk_ldofs.begin(),nk_ldofs.end());
		na_ldofs.insert(na_ldofs.end(),nu_ldofs.begin(),nu_ldofs.end());

		fsol_type fsol(youngsModulus,poissonsRatio);
		
		InnerEvaluation<diri_traits_type,diri_traits_type,fsol_type,DLP_tag> inner_evaluation_dlp(fsol);
		inner_evaluation_dlp(da_ldofs,displacements,pnts,-1.0,tmp1);
		
		InnerEvaluation<neum_traits_type,neum_traits_type,fsol_type,SLP_tag> inner_evaluation_slp(fsol);
		inner_evaluation_slp(na_ldofs,tractions,pnts,1.0,tmp2);

		u_inner=tmp1+tmp2;

		ofstream out(fileName.c_str());
        
        if(!out.is_open()) return -1;
        out.precision(12);
        //out.setf(std::ios::scientific);
        //write header
        out << "# vtk DataFile Version 2.0" << endl;
        out << "VTK inner evaluation" << endl;
        out << "ASCII" << endl;
        out << "DATASET STRUCTURED_POINTS" << endl;
		out << "DIMENSIONS " << N_pts << " " << N_pts << " " << N_pts << endl;
		out << "ORIGIN " << -0.5*l << " " << -0.5*l << " " << -0.5*l << endl;
		out << "SPACING " << delta << " " << delta << " " << delta << endl;

		out << "POINT_DATA " << N_pts*N_pts*N_pts << endl;
		out << "VECTORS displacement double" << endl;
		for(int j=0;j<N_pts*N_pts*N_pts; ++j){
			out << u_inner[3*j  ] << " "
				<< u_inner[3*j+1] << " "
				<< u_inner[3*j+2] << endl;
			//out << pnts[j][0] << " "
			//	<< pnts[j][1] << " "
			//	<< pnts[j][2] << endl;
			//printf("\n (%+.3lf %+.3lf %+.3lf): %+.8le %+.8le %+.8le", pnts[j][0],pnts[j][1],pnts[j][2], u_inner[3*j  ],u_inner[3*j+1],u_inner[3*j+2]);
		}
		out.close(); /**/
		return 0; 
	} 
}
