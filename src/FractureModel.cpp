/* 
 * File:   FractureModel.cpp
 * Author: David
 * 
 * Created on 20. Jan 2014, 13:49
 */

#include "FractureModel.h"

#include "types.h"

#include <cmath>
#include <cstdio>
#include <float.h>

namespace FractureSim{
	unsigned int FractureModel::cpVersion=0; // use 0 for original version, 1 for improved and 2 for new (experimental) version
	void applyCompressiveFactor(Eigen::Vector3d& K, double cf); // file-local helper function

	bool DynamicVelocityFractureModel::fractureCriterion(
		Eigen::Vector3d& dn, const Eigen::Vector3d& K_in, const Eigen::Vector3d& pos,
		const Eigen::Vector3d& n1, const Eigen::Vector3d& n2, const Eigen::Vector3d& n3,
		double blend_w, double old_th, double old_th_min, double old_th_max
		// n1 is surface normal ("out-of-plane"), n2 is crack-tip normal ("in-plane"), n3 is crack-tip tangent direction
	) const { // SIFs (K_I, K_II, K_III) are stored in K[0], K[1], K[2]
        Eigen::Vector3d K=K_in; // copy of input param - I know this is annoying, I tried to be sneaky by passing a const-ref, but then the following trick wouldn't work, and I was too lazy to change the function declaration
		applyCompressiveFactor(K,matMdl.compressiveFactor(pos));
        
		double Kc_p = matMdl.fractureToughness(pos); // evaluate fracture toughness
		double Keff2= K[0]*K[0]+K[1]*K[1]+K[2]*K[2]/(1.0-nu); // squared effective stress intensity

		if( Keff2 >= (Kc_p*Kc_p) ){ // do fracture
			double v = cR*(1-Kc_p*Kc_p/Keff2); // propagation velocity
			v=std::max(minVelocity, v);

			// find direction of propagation
			double phi=0.5*atan(2*K[2]/K[0]/(1-2*nu)); // angle phi is twist around in-plane normal)  from http://perso.crans.org/~verbeke/Cours_MAGIS/Cours%20site%20MAGIS/Fili%C3%A8re%20Endommagement%20et%20rupture%20des%20mtx%20et%20des%20structures/M%C3%A9canique%20de%20la%20rupture/MLR-session-4-Course-LEFM-MixedMode.pdf
			
			K[0]= /*(K[0]<0?-1.0:1.0)*/ sqrt( K[0]*K[0] + K[2]*K[2]/(1.0-nu)); // group K1 and K3 since these drive the crack forward, K2 drives the crack sideways
			// initial theta (optimal for homogeneous material)
			double th=0.0, thMin=0.0, thMax=0.0;
			th = getHoopAngleAndInterval(K_in,pos,thMin,thMax);
			if( FractureModel::cpVersion==2){ // new and experimental version
				if( std::abs(old_th) <= M_PI)     th = old_th*blend_w + (1.0-blend_w)*th;
				if( std::abs(old_th_min) <= M_PI) thMin = old_th_min*blend_w + (1.0-blend_w)*thMin;
				if( std::abs(old_th_max) <= M_PI) thMax = old_th_max*blend_w + (1.0-blend_w)*thMax;
			}
			th = projectKcGradient(th, thMin,thMax, pos,n1,n2);

			dn[0]=v*dt; // build output
			dn[1]=th;
			dn[2]=phi; // twist angle not used in the propagation step atm. (variation of th along the crack tip handles this)
			return true;
		}
		return false;
	}


	double DynamicVelocityFractureModel::getHoopAngleAndInterval(
		const Eigen::Vector3d& K_in, const Eigen::Vector3d& pos, double& thMin, double& thMax
	) const {
		Eigen::Vector3d K=K_in; // copy of input param
		applyCompressiveFactor(K,matMdl.compressiveFactor(pos));

		K[0]= sqrt( K[0]*K[0] + K[2]*K[2]/(1.0-nu)); // group K1 and K3 since these drive the crack forward, K2 drives the crack sideways
		double th=0.0;
		if(std::abs(K[1]) > DBL_EPSILON) th = 2*atan( (K[0] - sqrt(K[0]*K[0] + 8*K[1]*K[1])) / (4*K[1]) ); // angle of crack propagation (theta is rotation around crack-tip tangent)
		thMin=th; thMax=th;

		// compute valid interval ...
		double Kth, dKth_dth, Kc_p, th_a,th_b,u, Ka,Kb;
        Kc_p = matMdl.fractureToughness(pos); // make sure we compute this only when stress is high enough

        Kth = K[0]*cos(th*0.5)*cos(th*0.5)*cos(th*0.5) - 3.0*K[1]*sin(th*0.5)*cos(th*0.5)*cos(th*0.5);
		if( Kth < Kc_p ) return th; //abort (not enough stress)

        double a,b,c,d, fa,fb,fc, s,fs, tmp, eps=FLT_EPSILON;
        // find lower root bracketing interval [a,th]
        bool found=false, mflag;
        for(d=M_PI/24.0; d<M_PI && !found; d+=M_PI/24.0){
            b = th-d;
            Kb= K[0]*cos(b*0.5)*cos(b*0.5)*cos(b*0.5) - 3.0*K[1]*sin(b*0.5)*cos(b*0.5)*cos(b*0.5);
            if(Kb < Kc_p) found=true;
        }
        if( found ){
			a=b;fa=Kb-Kc_p;
			// run on interval [a,th] - don't modify th or Kth
			b=th; fb=Kth-Kc_p;
			if(std::abs(fa)<std::abs(fb)){ //swap (a,b)
				tmp = a; a = b; b = tmp;
				tmp =fa; fa=fb; fb= tmp;
			}
			c=a; mflag=true;
			// Brent-Dekker iteration
			while( std::abs(fa) > eps && std::abs(b-a)>eps ) {
				fc= -Kc_p+ K[0]*cos(c*0.5)*cos(c*0.5)*cos(c*0.5) - 3.0*K[1]*sin(c*0.5)*cos(c*0.5)*cos(c*0.5);
				if( std::abs(fa-fc)>eps && std::abs(fb-fc)>eps ){ // can use inv quad interp
					s = a*fb*fc/((fa-fb)*(fa-fc)) + b*fa*fc/((fb-fa)*(fb-fc)) + c*fa*fb/((fc-fa)*(fc-fb));
				}else{ // use lin interp
					s = b - fb*(b-a)/(fb-fa);
				}
				if( (s < 0.25*(3.0*a+b) || s > b) ||
					(mflag==true  && std::abs(s-b) >= 0.5*std::abs(b-c)) ||
					(mflag==false && std::abs(s-b) >= 0.5*std::abs(c-d)) ||
					(mflag==true  && std::abs(b-c) < eps) ||
					(mflag==false && std::abs(c-d) < eps)
				){
					s = 0.5*(a+b); mflag=true;
				}else mflag=false;
				fs= -Kc_p+ K[0]*cos(s*0.5)*cos(s*0.5)*cos(s*0.5) - 3.0*K[1]*sin(s*0.5)*cos(s*0.5)*cos(s*0.5);
				d=c; c=b; fc=fb;
				if( fa*fs < 0.0){ b=s; fb=fs; }
				else { a=s; fa=fs; }
				if(std::abs(fa)<std::abs(fb)){ //swap (a,b)
					tmp = a; a = b; b = tmp;
					tmp =fa; fa=fb; fb= tmp;
				}
			}
			th_a = asin(sin(b));
		}else{
            //printf("\nth_a not found, last angle %.4lf with K %.1le, Kc %.3le (Kth %.3le)",b,Kb, Kc_p, Kth);
			th_a = asin(sin(th-M_PI)); // if we could not find a proper intersection it is most likely because the stress intensity in insanely high
        }

        // find upper root bracketing interval [th,b]
        found=false;
        for(d=M_PI/24.0; d<M_PI && !found; d+=M_PI/24.0){
            b = th+d;
            Kb= K[0]*cos(b*0.5)*cos(b*0.5)*cos(b*0.5) - 3.0*K[1]*sin(b*0.5)*cos(b*0.5)*cos(b*0.5);
            if(Kb < Kc_p) found=true;
        }
        if( found ){
			// run on interval [th,b] - don't modify th or Kth
			a=th; fa=Kth-Kc_p; fb=Kb-Kc_p;
			if(std::abs(fa)<std::abs(fb)){ //swap (a,b)
				tmp = a; a = b; b = tmp;
				tmp =fa; fa=fb; fb= tmp;
			}
			c=a; mflag=true;
			// Brent-Dekker iteration
			while( std::abs(fa) > eps && std::abs(b-a)>eps ) {
				fc= -Kc_p+ K[0]*cos(c*0.5)*cos(c*0.5)*cos(c*0.5) - 3.0*K[1]*sin(c*0.5)*cos(c*0.5)*cos(c*0.5);
				if( std::abs(fa-fc)>eps && std::abs(fb-fc)>eps ){ // can use inv quad interp
					s = a*fb*fc/((fa-fb)*(fa-fc)) + b*fa*fc/((fb-fa)*(fb-fc)) + c*fa*fb/((fc-fa)*(fc-fb));
				}else{ // use lin interp
					s = b - fb*(b-a)/(fb-fa);
				}
				if( (s < 0.25*(3.0*a+b) || s > b) ||
					(mflag==true  && std::abs(s-b) >= 0.5*std::abs(b-c)) ||
					(mflag==false && std::abs(s-b) >= 0.5*std::abs(c-d)) ||
					(mflag==true  && std::abs(b-c) < eps) ||
					(mflag==false && std::abs(c-d) < eps)
				){
					s = 0.5*(a+b); mflag=true;
				}else mflag=false;
				fs= -Kc_p+ K[0]*cos(s*0.5)*cos(s*0.5)*cos(s*0.5) - 3.0*K[1]*sin(s*0.5)*cos(s*0.5)*cos(s*0.5);
				d=c; c=b; fc=fb;
				if( fa*fs < 0.0){ b=s; fb=fs; }
				else { a=s; fa=fs; }
				if(std::abs(fa)<std::abs(fb)){ //swap (a,b)
					tmp = a; a = b; b = tmp;
					tmp =fa; fa=fb; fb= tmp;
				}
			}
			th_b = asin(sin(a));
		}else{
            //printf("\nth_b not found, last angle %.4lf with K %.1le, Kc %.3le (Kth %.3le)",b,Kb, Kc_p, Kth);
            th_b = asin(sin(th+M_PI)); // if we could not find a proper intersection it is most likely because the stress intensity in insanely high
        }
               
        if(th_a>th_b){
            u=th_a ; th_a=th_b ;th_b=u ;
        }
		thMin=th_a; thMax=th_b;
		return th;
	}

	bool DynamicVelocityFractureModel::fractureCriterion(
			Eigen::Vector3d& dn, const Eigen::Vector3d& K
	) const { printf("\nsimple fracture criterion\n");
		// mixed I,II mode max tensile stress criterion from http://www.win.tue.nl/analysis/reports/rana07-23.pdf
		//   K1-->K[0], K2-->K[1], K3-->K[2]
        double den  = K[0]*K[0] + 12*K[1]*K[1] - K[0]*sqrt(K[0]*K[0] + 8*K[1]*K[1]);
        double Kmax = 4*M_SQRT2*K[1]*K[1]*K[1]*( K[0] + 3*sqrt(K[0]*K[0]+8*K[1]*K[1]) ) / ( den*sqrt(den) );
        Kmax=std::abs(Kmax)* ((K[0]>0.0) - (K[0]<0.0)); // not sure what the sign of Kmax indicates, but its not tensile vs. compressive as for K1
        //double Kmax=K1;
		//printf("\nK1,2,3 (%.1le, %.1le, %.1le) max %.1le",K1,K2,K3,Kmax);
		if( Kmax >= Kc /*|| std::abs(K3)>=Kc*/ ){
			double th,phi; // angle of crack propagation
            th = 2*atan( (K[0] - sqrt(K[0]*K[0] + 8*K[1]*K[1])) / (4*K[1]) ); // from http://www.win.tue.nl/analysis/reports/rana07-23.pdf
			phi= 0.5*atan(2*K[2]/K[0]/(1-2*nu)); // from http://perso.crans.org/~verbeke/Cours_MAGIS/Cours%20site%20MAGIS/Fili%C3%A8re%20Endommagement%20et%20rupture%20des%20mtx%20et%20des%20structures/M%C3%A9canique%20de%20la%20rupture/MLR-session-4-Course-LEFM-MixedMode.pdf
			dn[1]=th;
			dn[2]=phi;
			// dynamic crack velocity
			double v = cR*(1-Kc*Kc/(Kmax*Kmax));
			//clamp v to a minimal value to avoid very short propagation steps
			v=std::max(minVelocity, v);
            //th=std::max(-M_PI_4,std::min(M_PI_4,th));
			dn[0]=v*dt;
			//printf("%% dynamic crack propagation at %.3lf cR\n",(1-Kc*Kc/(Kmax*Kmax)));

			return true;
		} // How to consider compressive fracture or pure mode III propagation (tearing)?
		return false;
	}
    
	double DynamicVelocityFractureModel::projectKcGradient(
			double th, double th_a, double th_b, 
			const Eigen::Vector3d& pos, const Eigen::Vector3d& n1, const Eigen::Vector3d& n2
	) const{
		Eigen::Vector3d dKc;
		double Kc_p = matMdl.fractureToughness(pos), eps=FLT_EPSILON;
		matMdl.fractureToughnessGradient(dKc, pos);

		if( dKc.norm() < eps) return th; //abort (locally homogeneous material)

		// represent dKc in (n1,n2) basis
        double u = dKc.dot(n2);
        double v = dKc.dot(n1);
		double w = 1.0/(1.0 + sqrt(u*u+v*v) / Kc_p *maxStep ); // compare magnitude of (projected) Kc-gradient to (Kc/length) with some characteristic length << make this the bem resolution
		u=atan2(-v,-u);
        if(u<th_a)      u=th_a;
		else if(u>th_b) u=th_b;
		th=(1.0-w)*u + w*th;
		//printf("\n%% weight %.3lf (mxs %.3le) grad-norm %.1lg, limits (%.3lf , %.3lf)", w, maxStep, sqrt(u*u+v*v),th_a,th_b);
		//printf("\n%% new theta %.4lf\t candidate was %.4lf",th,u);
        return th;
	}


    double DynamicVelocityFractureModel::findValidKcGradient(
            double th, const Eigen::Vector3d& K,
            const Eigen::Vector3d& pos, const Eigen::Vector3d& n1, const Eigen::Vector3d& n2
    ) const{
		printf("\nOLD VERSION, DON'T USE - DynamicVelocityFractureModel::findValidKcGradient\n");
        // idea: based on K, find the range of th such that Kth >= Kc(pos)
        // then find which th in this range provides the most decrease of Kc
        // according to the gradient of Kc at pos
        Eigen::Vector3d dKc;
        double Kth, dKth_dth, Kc_p, th_a,th_b,u, eps=FLT_EPSILON; //1e-3;
        double Ka,Kb; //for testing
        Kc_p = matMdl.fractureToughness(pos); // make sure we compute only when stress is high enough
		matMdl.fractureToughnessGradient(dKc, pos);

        Kth = K[0]*cos(th*0.5)*cos(th*0.5)*cos(th*0.5) - 3.0*K[1]*sin(th*0.5)*cos(th*0.5)*cos(th*0.5);
		if( Kth < (Kc_p-10.0*eps) || dKc.norm() < /**0.1*Kc_p/*/eps/**/) return th; //abort (not enough stress or locally homogeneous material)

        double a,b,c,d, fa,fb,fc, s,fs, tmp;
        // find lower root bracketing interval [a,th]
        bool found=false, mflag;
        for(d=M_PI/24.0; d<M_PI && !found; d+=M_PI/24.0){
            b = th-d;
            Kb= K[0]*cos(b*0.5)*cos(b*0.5)*cos(b*0.5) - 3.0*K[1]*sin(b*0.5)*cos(b*0.5)*cos(b*0.5);
            if(Kb < Kc_p) found=true;
        }
        if(!found){
            //printf("\nth_a not found, last angle %.4lf with K %.1le, Kc %.3le (Kth %.3le)",b,Kb, Kc_p, Kth);
            return th; //abort
        }
        a=b;fa=Kb-Kc_p;
        // run on interval [a,th] - don't modify th or Kth
        b=th; fb=Kth-Kc_p;
        if(std::abs(fa)<std::abs(fb)){ //swap (a,b)
            tmp = a; a = b; b = tmp;
            tmp =fa; fa=fb; fb= tmp;
        }
        c=a; mflag=true;
        // Brent-Dekker iteration
        while( std::abs(fa) > eps && std::abs(b-a)>eps ) {
            fc= -Kc_p+ K[0]*cos(c*0.5)*cos(c*0.5)*cos(c*0.5) - 3.0*K[1]*sin(c*0.5)*cos(c*0.5)*cos(c*0.5);
            if( std::abs(fa-fc)>eps && std::abs(fb-fc)>eps ){ // can use inv quad interp
                s = a*fb*fc/((fa-fb)*(fa-fc)) + b*fa*fc/((fb-fa)*(fb-fc)) + c*fa*fb/((fc-fa)*(fc-fb));
            }else{ // use lin interp
                s = b - fb*(b-a)/(fb-fa);
            }
            if( (s < 0.25*(3.0*a+b) || s > b) ||
                (mflag==true  && std::abs(s-b) >= 0.5*std::abs(b-c)) ||
                (mflag==false && std::abs(s-b) >= 0.5*std::abs(c-d)) ||
                (mflag==true  && std::abs(b-c) < eps) ||
                (mflag==false && std::abs(c-d) < eps)
            ){
                s = 0.5*(a+b); mflag=true;
            }else mflag=false;
            fs= -Kc_p+ K[0]*cos(s*0.5)*cos(s*0.5)*cos(s*0.5) - 3.0*K[1]*sin(s*0.5)*cos(s*0.5)*cos(s*0.5);
            d=c; c=b; fc=fb;
            if( fa*fs < 0.0){ b=s; fb=fs; }
            else { a=s; fa=fs; }
            if(std::abs(fa)<std::abs(fb)){ //swap (a,b)
                tmp = a; a = b; b = tmp;
                tmp =fa; fa=fb; fb= tmp;
            }
        }
        th_a = asin(sin(b));
        
        // find upper root bracketing interval [th,b]
        found=false;
        for(d=M_PI/24.0; d<M_PI && !found; d+=M_PI/24.0){
            b = th+d;
            Kb= K[0]*cos(b*0.5)*cos(b*0.5)*cos(b*0.5) - 3.0*K[1]*sin(b*0.5)*cos(b*0.5)*cos(b*0.5);
            if(Kb < Kc_p) found=true;
        }
        if(!found){
            //printf("\nth_b not found");
            return th; //abort
        }
        // run on interval [th,b] - don't modify th or Kth
        a=th; fa=Kth-Kc_p; fb=Kb-Kc_p;
        if(std::abs(fa)<std::abs(fb)){ //swap (a,b)
            tmp = a; a = b; b = tmp;
            tmp =fa; fa=fb; fb= tmp;
        }
        c=a; mflag=true;
        // Brent-Dekker iteration
        while( std::abs(fa) > eps && std::abs(b-a)>eps ) {
            fc= -Kc_p+ K[0]*cos(c*0.5)*cos(c*0.5)*cos(c*0.5) - 3.0*K[1]*sin(c*0.5)*cos(c*0.5)*cos(c*0.5);
            if( std::abs(fa-fc)>eps && std::abs(fb-fc)>eps ){ // can use inv quad interp
                s = a*fb*fc/((fa-fb)*(fa-fc)) + b*fa*fc/((fb-fa)*(fb-fc)) + c*fa*fb/((fc-fa)*(fc-fb));
            }else{ // use lin interp
                s = b - fb*(b-a)/(fb-fa);
            }
            if( (s < 0.25*(3.0*a+b) || s > b) ||
                (mflag==true  && std::abs(s-b) >= 0.5*std::abs(b-c)) ||
                (mflag==false && std::abs(s-b) >= 0.5*std::abs(c-d)) ||
                (mflag==true  && std::abs(b-c) < eps) ||
                (mflag==false && std::abs(c-d) < eps)
            ){
                s = 0.5*(a+b); mflag=true;
            }else mflag=false;
            fs= -Kc_p+ K[0]*cos(s*0.5)*cos(s*0.5)*cos(s*0.5) - 3.0*K[1]*sin(s*0.5)*cos(s*0.5)*cos(s*0.5);
            d=c; c=b; fc=fb;
            if( fa*fs < 0.0){ b=s; fb=fs; }
            else { a=s; fa=fs; }
            if(std::abs(fa)<std::abs(fb)){ //swap (a,b)
                tmp = a; a = b; b = tmp;
                tmp =fa; fa=fb; fb= tmp;
            }
        }
        th_b = asin(sin(a));
               
        if(th_a>th_b){
            u=th_a ; th_a=th_b ;th_b=u ;
        }
        // represent dKc in (n1,n2) basis
        u = dKc.dot(n2);
        double v = dKc.dot(n1);
		double w = 1.0/(1.0 + sqrt(u*u+v*v) / Kc_p *maxStep ); // compare magnitude of (projected) Kc-gradient to (Kc/length) with some characteristic length << make this the bem resolution
		u=atan2(-v,-u);
        if(u<th_a)      u=th_a;
		else if(u>th_b) u=th_b;
		th=(1.0-w)*u + w*th; //printf("\n %.3lf (mxs %.3le)", w, maxStep);
//        printf("\nnew theta %.4lf\t candidate was %.4lf",th,u);
        return th;
    }
    
    double DynamicVelocityFractureModel::findMaxKthMinusKc(
        double th, const Eigen::Vector3d& K,
        const Eigen::Vector3d& pos, const Eigen::Vector3d& n1, const Eigen::Vector3d& n2
    ) const{
		printf("\nOLD VERSION, DON'T USE - DynamicVelocityFractureModel::findMaxKthMinusKc\n");
		// instead of computing the direction of max. hoop stress (mode I/II) as earlier versions,
		// search for the direction that maximizes the difference of hoop stress to fracture toughness
		// i.e. we want theta such as to find max(Kmax(theta)-Kc(pos+0.5*v(theta)))
		// (at about 1/2 a timestep (or a complete one?) away from the current position)

		// Based on eq. (33), (34) in http://www.win.tue.nl/analysis/reports/rana07-23.pdf
		// first compute the angle theta which maximizes mixed I/II SIF as
		// th = 2 atan( ( K1 - sqrt(K1^2 + 8 K2^2) ) / (4 K2) )
		// for any given theta, the mixed-mode SIF is
		// K_th = K1 cos(th/2)^3 - 3 K2 sin(th/2) cos(th/2)^2
		// then check the gradient of (toughness - K_th) to maximize the difference
        // run a Brent-Dekker root finding on d/d_theta( toughness - K_th )
        Eigen::Vector3d dKc;
        double Kth, dKth_dth, Kc_p, dKc_dn;
        double a,b,c,d, fa,fb,fc, s,fs, tmp, eps = FLT_EPSILON;
        // set up Brent-Dekker algorithm --> find root bracketing interval
        a = th;
        fa= evalKthMinusKcGradient(
            Kth, dKth_dth, Kc_p, dKc_dn, dKc,
            a, K, pos, n1, n2
        );
        bool found=false, mflag;
        for(d=M_PI/24.0; d<M_PI_2 && !found; d+=M_PI/24.0){
            b = th+d;
            fb= evalKthMinusKcGradient(
                Kth, dKth_dth, Kc_p, dKc_dn, dKc,
                b, K, pos, n1, n2
            );
            if(fa*fb<0.0) found=true;
            else{
                b = th-d;
                fb= evalKthMinusKcGradient(
                    Kth, dKth_dth, Kc_p, dKc_dn, dKc,
                    b, K, pos, n1, n2
                );
                if(fa*fb<0.0) found=true;
            }
        }
        //if(!found)printf("\nsearch interval [%.4lf, %.4lf] found=%d",a,b,found);
        if(std::abs(fa)<std::abs(fb)){ //swap (a,b)
            tmp = a; a = b; b = tmp;
            tmp =fa; fa=fb; fb= tmp;
        }
        c=a; mflag=true;
        // Brent-Dekker iteration
        while( found && std::abs(fa) > eps && std::abs(b-a)>eps ) {
            fc= evalKthMinusKcGradient(
                Kth, dKth_dth, Kc_p, dKc_dn, dKc,
                c, K, pos, n1, n2
            );
            if( std::abs(fa-fc)>eps && std::abs(fb-fc)>eps ){ // can use inv quad interp
                s = a*fb*fc/((fa-fb)*(fa-fc)) + b*fa*fc/((fb-fa)*(fb-fc)) + c*fa*fb/((fc-fa)*(fc-fb));
            }else{ // use lin interp
                s = b - fb*(b-a)/(fb-fa);
            }
            if( (s < 0.25*(3.0*a+b) || s > b) ||
                (mflag==true  && std::abs(s-b) >= 0.5*std::abs(b-c)) ||
                (mflag==false && std::abs(s-b) >= 0.5*std::abs(c-d)) ||
                (mflag==true  && std::abs(b-c) < eps) ||
                (mflag==false && std::abs(c-d) < eps)
            ){
                s = 0.5*(a+b); mflag=true;
            }else mflag=false;
            fs= evalKthMinusKcGradient(
                Kth, dKth_dth, Kc_p, dKc_dn, dKc,
                s, K, pos, n1, n2
            );
            d=c; c=b; fc=fb;
            if( fa*fs < 0.0){ b=s; fb=fs; }
            else { a=s; fa=fs; }
            if(std::abs(fa)<std::abs(fb)){ //swap (a,b)
                tmp = a; a = b; b = tmp;
                tmp =fa; fa=fb; fb= tmp;
            }
        }
        if(found){
            // theta must always be in [-pi/2,+pi/2]
            th = asin(sin(a));
        }
        return th;
    }
    
    
    double DynamicVelocityFractureModel::evalKthMinusKcGradient(
        double& Kth, double& dKth_dth, double& Kc_p, double& dKc_dn, Eigen::Vector3d& dKc,
        double th, const Eigen::Vector3d& K,
        const Eigen::Vector3d& pos, const Eigen::Vector3d& n1, const Eigen::Vector3d& n2
    ) const{
            // function and angular derivative evaluation **************************
        Eigen::Vector3d r  = n1*sin(th) + n2*cos(th); // direction given by theta in (n1,n2) plane
		Eigen::Vector3d n  = n1*cos(th) - n2*sin(th); // normal to r in (n1,n2) plane
        Eigen::Vector3d x  = pos + 0.5*cR*dt*r;
        Kc_p = matMdl.fractureToughness(x); // evaluate fracture toughness
		matMdl.fractureToughnessGradient(dKc, x); // evaluate gradient of fracture toughness
        dKc_dn = dKc.dot(n); // directional derivative of toughness along the normal to theta ~> locally angular derivative around theta
        Kth = K[0]*cos(th*0.5)*cos(th*0.5)*cos(th*0.5) - 3.0*K[1]*sin(th*0.5)*cos(th*0.5)*cos(th*0.5);
        dKth_dth = (3.0*K[1]*sin(th)*sin(th*0.5)*0.5 - cos(th*0.5)*(3.0*K[1]*cos(th)*0.25 + 3.0*K[0]*sin(th)*0.25 + 3.0*K[1]*0.25));
        // *********************************************************************
        return dKc_dn - dKth_dth;
    }

	bool ConstantVelocityFractureModel::fractureCriterion(
			Eigen::Vector3d& dn, const Eigen::Vector3d& K
	) const {
		// mixed I,II mode max tensile stress criterion from http://www.win.tue.nl/analysis/reports/rana07-23.pdf
        double den  = K[0]*K[0] + 12*K[1]*K[1] - K[0]*sqrt(K[0]*K[0] + 8*K[1]*K[1]);
        double Kmax = 4*M_SQRT2*K[1]*K[1]*K[1]*( K[0] + 3*sqrt(K[0]*K[0]+8*K[1]*K[1]) ) / ( den*sqrt(den) );
        Kmax=std::abs(Kmax)* ((K[0]>0.0) - (K[0]<0.0));
        
		if( Kmax >= Kc ){
			double th; // angle of crack propagation
            th = 2*atan( (K[0] - sqrt(K[0]*K[0] + 8*K[1]*K[1])) / (4*K[1]) ); // from http://www.win.tue.nl/analysis/reports/rana07-23.pdf
			//simple implementation (constant crack velocity)
			dn[0]=v*dt;
			dn[1]=th;
			dn[2]=0.0; // no twist angle
			return true;
		}else{
			return false;
		}
	}

	void applyCompressiveFactor(Eigen::Vector3d& K, double cf){
		double th_mx,th_mi,K_mx,K_mi;
		if( (K[0]<0 || cf<1.0)  && cf>=0.0){ // cf must be positive to make sense, otherwise don't change anything

			if(std::abs(K[1]) > DBL_EPSILON){
				th_mx= 2*atan( (K[0] - sqrt(K[0]*K[0] + 8*K[1]*K[1])) / (4*K[1]) );
				th_mi= 2*atan( (K[0] + sqrt(K[0]*K[0] + 8*K[1]*K[1])) / (4*K[1]) );
				K_mx = K[0]*cos(th_mx*0.5)*cos(th_mx*0.5)*cos(th_mx*0.5) - 3.0*K[1]*sin(th_mx*0.5)*cos(th_mx*0.5)*cos(th_mx*0.5);
				K_mi = K[0]*cos(th_mi*0.5)*cos(th_mi*0.5)*cos(th_mi*0.5) - 3.0*K[1]*sin(th_mi*0.5)*cos(th_mi*0.5)*cos(th_mi*0.5);
			}else{
				K_mx=K[0]; K_mi=K[0]; // reduces (K_mx < (-K_mi/cf)) to (K[0]<0)
			}

			if( K_mx < (-K_mi/cf) ){ // flip sign if we want to go for the minimum instead
				// cf means that compressive toughness is cf times tensile toughness
				K[0]*=-1.0/cf;
				K[1]*=-1.0/cf; // changed my mind on these - we do need to reduce all SIFs if we commit to going for compressive propagation
				K[2]*=-1.0/cf;
			}
        }
	}
}
