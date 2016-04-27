/* 
 * File:   MaterialModel.h
 * Author: David
 *
 * Created on 24. März 2014
 */

#ifndef MATERIALMODEL_H
#define	MATERIALMODEL_H

#include <Eigen/Dense>
#include <string>

namespace FractureSim{
	// factory method
	class MaterialModel;
	MaterialModel* createMaterialModel(std::string spec, double youngsMod, double poissonsRatio, double density, double strength, double toughness, double compressiveFactor);

	//abstract class, use MaterialModelImpl.h for useable subclasses
	class MaterialModel {
    public:
		MaterialModel(double E=1.0, double nu=0.3, double rho=1e3)
			: E(E), nu(nu), rho(rho) {}
		MaterialModel(const MaterialModel& ori)
			: E(ori.E), nu(ori.nu), rho(ori.rho) {}
		virtual MaterialModel* clone() const =0;
		virtual ~MaterialModel(){}

		virtual double tensileStrength(const Eigen::Vector3d& x) const =0; // returns the critical tensile stress at which a fracture will start 
		virtual double compressiveFactor(const Eigen::Vector3d& x) const =0; // returns a factor specifying how much stronger & tougher the material is under compression as opposed to tension

		virtual double fractureToughness(const Eigen::Vector3d& x) const =0; // returns the critical stress intensity at which cracks will propagate
		virtual void   fractureToughnessGradient(Eigen::Vector3d& dKc, const Eigen::Vector3d& x) const =0;

		virtual double getE() const {return E;}
		virtual double getNu() const {return nu;}
		virtual double getDensity() const {return rho;}
		virtual void setE(double value){E=value;}
		virtual void setNu(double value){nu=value;}
		virtual void setDensity(double value){rho=value;}
	protected:
		double E,nu,rho; // store Young's modulus, Poisson's ratio and density
	};
}

#endif	/* MATERIALMODEL_H */
