/* 
 * File:   MaterialModel.h
 * Author: David
 *
 * Created on 27 June 2014
 */

#ifndef MATERIALMODELIMPL_H
#define	MATERIALMODELIMPL_H

#include "MaterialModel.h"

namespace FractureSim{

	class HomogeneousMaterialModel : public MaterialModel {
	public:
		HomogeneousMaterialModel(
			double E=1.0, double nu=0.3, double rho=1e3,
			double strength=1.0, double toughness=1.0, double compressiveFactor=3.0
		) : MaterialModel(E,nu,rho), Sc(strength), Kc(toughness), cf(compressiveFactor) {}
		HomogeneousMaterialModel(const HomogeneousMaterialModel& ori)
			: MaterialModel(ori.E,ori.nu,ori.rho), Sc(ori.Sc), Kc(ori.Kc), cf(ori.cf) {}
		virtual HomogeneousMaterialModel* clone() const {return(new HomogeneousMaterialModel(*this)); }
		virtual ~HomogeneousMaterialModel(){}

		inline virtual double tensileStrength(const Eigen::Vector3d& x) const {return Sc;}
		inline virtual double compressiveFactor(const Eigen::Vector3d& x) const {return cf;}

		inline virtual double fractureToughness(const Eigen::Vector3d& x) const { return Kc;}
		inline virtual void   fractureToughnessGradient(Eigen::Vector3d& dKc, const Eigen::Vector3d& x) const { dKc.setZero();}
	protected:
		double Kc,Sc,cf;
	};

	class DummyMaterialModel : public HomogeneousMaterialModel {
	public:
		DummyMaterialModel(
			double E=1.0, double nu=0.3, double rho=1e3,
			double strength=1.0, double toughness=1.0, double compressiveFactor=3.0,
			double toughnessFactor=1.0
		) : HomogeneousMaterialModel(E,nu,rho,strength,toughness,compressiveFactor), KcDelta(toughnessFactor) {}
		DummyMaterialModel(const DummyMaterialModel& ori)
			: HomogeneousMaterialModel(ori.E,ori.nu,ori.rho,ori.Sc,ori.Kc,ori.cf), KcDelta(ori.KcDelta) {}
		virtual DummyMaterialModel* clone() const {return(new DummyMaterialModel(*this)); }
		virtual ~DummyMaterialModel(){}

		virtual double fractureToughness(const Eigen::Vector3d& x) const;
		virtual void   fractureToughnessGradient(Eigen::Vector3d& dKc, const Eigen::Vector3d& x) const;
	protected:
		double KcDelta;
	};

	class LayeredMaterialModel : public HomogeneousMaterialModel {
	public:
		LayeredMaterialModel(
			double E=1.0, double nu=0.3, double rho=1e3,
			double strength=1.0, double compressiveFactor=3.0,
			double toughnessMin=1.0, double toughnessMax=1.0,
			double frequency=1.0, double phase=0.0,
            Eigen::Vector3d normal = Eigen::Vector3d::UnitX(),
            double strengthMax=1.0,	double strengthFreq=0.0, double strengthPhase=M_PI,
            Eigen::Vector3d strengthNormal = Eigen::Vector3d::UnitX()
		) : HomogeneousMaterialModel(E,nu,rho,strength,toughnessMin,compressiveFactor),
			KcHard(toughnessMax), f(frequency), p(phase), n(normal),
            ScHard(strengthMax), sf(strengthFreq), sp(strengthPhase), sn(strengthNormal) {}
		LayeredMaterialModel(const LayeredMaterialModel& ori)
			: HomogeneousMaterialModel(ori.E,ori.nu,ori.rho,ori.Sc,ori.Kc,ori.cf),
			KcHard(ori.KcHard), f(ori.f) , p(ori.p) , n(ori.n),
            ScHard(ori.ScHard),sf(ori.sf),sp(ori.sp),sn(ori.sn){}
		virtual LayeredMaterialModel* clone() const {return(new LayeredMaterialModel(*this)); }
		virtual ~LayeredMaterialModel(){}

        virtual double tensileStrength(const Eigen::Vector3d& x) const;
		virtual double fractureToughness(const Eigen::Vector3d& x) const;
		virtual void   fractureToughnessGradient(Eigen::Vector3d& dKc, const Eigen::Vector3d& x) const;
	protected:
		double KcHard, f ,p ;
		double ScHard, sf,sp;
		Eigen::Vector3d n,sn;
	};

	class TubedMaterialModel : public HomogeneousMaterialModel {
	public:
		TubedMaterialModel(
			double E=1.0, double nu=0.3, double rho=1e3,
			double strength=1.0, double compressiveFactor=3.0,
			double toughnessMin=1.0, double toughnessMax=1.0,
			double frequency=1.0, double phase=0.0, Eigen::Vector3d normal = Eigen::Vector3d::UnitX()
		) : HomogeneousMaterialModel(E,nu,rho,strength,toughnessMin,compressiveFactor),
			KcHard(toughnessMax) {}
		TubedMaterialModel(const TubedMaterialModel& ori)
			: HomogeneousMaterialModel(ori.E,ori.nu,ori.rho,ori.Sc,ori.Kc,ori.cf),
			KcHard(ori.KcHard) {}
		virtual TubedMaterialModel* clone() const {return(new TubedMaterialModel(*this)); }
		virtual ~TubedMaterialModel(){}

		virtual double fractureToughness(const Eigen::Vector3d& x) const;
		virtual void   fractureToughnessGradient(Eigen::Vector3d& dKc, const Eigen::Vector3d& x) const;
	protected:
		double KcHard;
		//Eigen::Vector3d n;
	};

	// use an OpenVDB sparse grid data file as the material's toughness and strength maps
	class VDBMaterialModel : public HomogeneousMaterialModel {
	public:
		VDBMaterialModel(
			double E=1.0, double nu=0.3, double rho=1e3,
			double strength=1.0, double toughness=1.0, double compressiveFactor=3.0,
            std::string file="material.vdb", std::string toughnessGrid="toughness", std::string strengthGrid="strength",
            bool udf=false, double scale=1.0, double voxelSize=1.0
        )   : HomogeneousMaterialModel(E,nu,rho,strength,toughness,compressiveFactor),
            data(NULL), vdbFile(file), kName(toughnessGrid), sName(strengthGrid),
            useUDF(udf), scale(scale), voxelSize(voxelSize) {}
		VDBMaterialModel(const VDBMaterialModel& ori)
            : HomogeneousMaterialModel(ori.E,ori.nu,ori.rho,ori.Sc,ori.Kc,ori.cf),
            data(NULL), vdbFile(ori.vdbFile), kName(ori.kName), sName(ori.sName),
            useUDF(ori.useUDF), scale(ori.scale), voxelSize(ori.voxelSize) {}
		virtual VDBMaterialModel* clone() const {return(new VDBMaterialModel(*this)); }
		virtual ~VDBMaterialModel(){ if(data!=NULL) delete data; }

        virtual double tensileStrength(const Eigen::Vector3d& x) const;
        virtual double fractureToughness(const Eigen::Vector3d& x) const;
		virtual void   fractureToughnessGradient(Eigen::Vector3d& dKc, const Eigen::Vector3d& x) const;
	protected:
        std::string vdbFile, kName, sName;
        double scale, voxelSize;
        bool useUDF;
    private:
        void init() const;
        class vdbData; //fwd. decl
        mutable vdbData* data;
		mutable bool haveToughnessMap, haveStrengthMap;
		template <class Vec3> Vec3 toPeriodicIndexSpace(const Eigen::Vector3d& x) const;
	};
}
#endif	/* MATERIALMODELIMPL_H */
