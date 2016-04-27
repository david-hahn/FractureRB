/* 
 * File:   FractureModel.h
 * Author: David
 *
 * Created on 20. Jänner 2014, 13:49
 */

#ifndef FRACTUREMODEL_H
#define	FRACTUREMODEL_H

#include <cstdio>
#include <Eigen/Dense>
#include "MaterialModel.h"

namespace FractureSim{

	//abstract class
    class FractureModel {
    public:
		// this affects the way we compute the crack propagation direction
		// in the FractureModel and SubsampledCrackTip classes
		// use 0 for the original version (linear blending to the previous time-step)
		//     1 for the improved version (quadratic blending)
		// and 2 for the experimental version (blending only affects the homogeneous material case,
		//       the toughness gradient is applied later, and the valid direction interval is smoothed)
		static unsigned int cpVersion;

		FractureModel(){}
		virtual ~FractureModel(){}

		/* Evaluate the fracture criterion for a crack-tip segment
		 * given stress intensity factors sifs:=(K1,K2,K3)
		 * in the local coordinate system such that
		 * K1 -> opening, K2 -> shearing, K3 -> tearing
		 * Returns true if the crack segment should propagate or
		 * false if the segment should arrest.
		 * In dn, the crack advance will be stored, if true is returned.
		 * Components of dn are (dr/dt, theta, phi)
		 */
		virtual bool fractureCriterion(
			Eigen::Vector3d& dn, const Eigen::Vector3d& sifs
		) const = 0; //abstract

		/* Newer version for spatially variable fracture toughness.
		 * Additional input parameters:
		 * Coordinates where the fracture criterion is evaluated pos:=(x,y,z)
		 * and the local coordinate system (n1,n2,n3) in which SIFs are given,
		 * such that K1 is oriented in n1 direction, etc.
		 */
		virtual bool fractureCriterion(
			Eigen::Vector3d& dn, const Eigen::Vector3d& sifs, const Eigen::Vector3d& pos,
			const Eigen::Vector3d& n1, const Eigen::Vector3d& n2, const Eigen::Vector3d& n3,
			double blend_w=0.0, double old_th=0.0, double old_th_min=0.0, double old_th_max=0.0 // params for cpVersion==2
		) const = 0; //abstract

		/* Computes the homogeneous out-of-plane propagation direction and valid interval
		 * i.e. the direction a crack-front marker would propagate under the given SIFs
		 * in the absence of a toughness gradient (this does NOT evaluate a fracture criterion)
		 */
		virtual double getHoopAngleAndInterval(const Eigen::Vector3d& K_in, const Eigen::Vector3d& pos, double& thMin, double& thMax) const = 0;

		/* NEW FOR FractureRB:
		 * return the current timestep
		 */
		inline double getTimeStep(){ return dt; }

	protected:
		double dt; // NEW FOR FractureRB: moved to base class
	};

	class ConstantVelocityFractureModel : public FractureModel {
	public:
		/* The crack will propagate with velocity v_ in the direction
		 * of max. hoop stress, if the corresponding SIF Kmax >= Kc_.
		 * Only mode I and II fracture supported! Compressive fracture not supported!
		 */
		ConstantVelocityFractureModel(double Kc_=1.0, double v_=1.0, double dt_=1.0){
			Kc=Kc_; v=v_; dt=dt_;
		}
		
		virtual bool fractureCriterion(
			Eigen::Vector3d& dn, const Eigen::Vector3d& sifs
		) const;

		virtual bool fractureCriterion(
			Eigen::Vector3d& dn, const Eigen::Vector3d& sifs, const Eigen::Vector3d&,
			const Eigen::Vector3d& n1, const Eigen::Vector3d& n2, const Eigen::Vector3d& n3,
			double blend_w=0.0, double old_th=0.0, double old_th_min=0.0, double old_th_max=0.0
		) const {
			printf("\nNOT IMPLEMENTED (ConstantVelocityFractureModel::fractureCriterion)\n");
			return false;
		}
		virtual double getHoopAngleAndInterval(const Eigen::Vector3d& K_in, const Eigen::Vector3d& pos, double& thMin, double& thMax) const{
			printf("\nNOT IMPLEMENTED (ConstantVelocityFractureModel::getHoopAngleAndInterval)\n");
			return 0.0;
		}

		virtual ~ConstantVelocityFractureModel(){}
	protected:
		double Kc, v;

	};

	/* The crack will propagate with velocity
	 * v = cR*(1-Kc*Kc/(Kmax*Kmax))
	 * where cR is the Rayleigh wave speed
	 * if Kmax > Kc, Kmax is max. hoop stress in-plane (K1,K2)
	 */
	class DynamicVelocityFractureModel : public FractureModel {
	public:
		DynamicVelocityFractureModel(MaterialModel& matModel, double maxStepSize)
			: matMdl(matModel)
		{
			Kc=-1.0; //no longer used
			maxStep = maxStepSize;
			E = matMdl.getE();
			nu= matMdl.getNu();
			cR=0.9*sqrt(E/(2.0*(1.0+nu)*matMdl.getDensity())); // approx., see Gross & Seelig, Fracture Mechanics, Springer 2011
			dt=maxStep / cR;
			minVelocity = 0.1*cR;
            printf("\n%% ... timestep is %.3lg, Rayleigh wave speed is %.3lg",dt, cR);
		}

		virtual bool fractureCriterion(
			Eigen::Vector3d& dn, const Eigen::Vector3d& sifs
		) const;

		virtual bool fractureCriterion(
			Eigen::Vector3d& dn, const Eigen::Vector3d& sifs, const Eigen::Vector3d& pos,
			const Eigen::Vector3d& n1, const Eigen::Vector3d& n2, const Eigen::Vector3d& n3,
			double blend_w=0.0, double old_th=0.0, double old_th_min=0.0, double old_th_max=0.0 // params for cpVersion==2
		) const;

		//for cpVersion==2
		virtual double getHoopAngleAndInterval(const Eigen::Vector3d& K_in, const Eigen::Vector3d& pos, double& thMin, double& thMax) const;

		virtual ~DynamicVelocityFractureModel(){}
	protected:
		double Kc, cR, E, nu, minVelocity, maxStep;
		const MaterialModel& matMdl;
        
		double projectKcGradient(
			double th, double thMin, double thMax, 
			const Eigen::Vector3d& pos, const Eigen::Vector3d& n1, const Eigen::Vector3d& n2
		) const;
        double findValidKcGradient(
            double th, const Eigen::Vector3d& K,
            const Eigen::Vector3d& pos, const Eigen::Vector3d& n1, const Eigen::Vector3d& n2
        ) const;
        double findMaxKthMinusKc(
            double th, const Eigen::Vector3d& K,
            const Eigen::Vector3d& pos, const Eigen::Vector3d& n1, const Eigen::Vector3d& n2
        ) const;
        double evalKthMinusKcGradient(
            double& Kth, double& dKth_dth, double& Kc_p, double& dKc_dn, Eigen::Vector3d& dKc,
            double th, const Eigen::Vector3d& K,
            const Eigen::Vector3d& pos, const Eigen::Vector3d& n1, const Eigen::Vector3d& n2
        ) const;
	};
}

#endif	/* FRACTUREMODEL_H */
