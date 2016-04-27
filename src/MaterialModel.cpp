/* 
 * File:   MaterialModel.cpp
 * Author: David
 * 
 * Created on 24. March 2014
 */

#include "MaterialModelImpl.h"

#include <cmath>
#include <cstdio>
using namespace std;

#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
namespace vdb = openvdb::v2_2_0;

namespace FractureSim{
	MaterialModel* createMaterialModel(string spec, double youngsMod, double poissonsRatio, double density, double strength, double toughness, double compress){
		if(spec.compare("?")==0){ // print options
			printf(
				"\nthe following material model specifications are implemented:\n"
				"default ... homogeneous material\n"
				"'lay(nx,ny,nz,f,p,K)' ... layered material with normal (nx,ny,nz), frequency f, phase p and secondary toughness K\n"
				"'test' ... dummy material model (used for testing)\n" );
		}
		if(spec.compare(0,3,"lay")==0){ // format: "lay(nx,ny,nz,f,p,Kc_hard)" or "lay(nx,ny,nz,f,p,K_hard/snx,sny,snz,sf,sp,Sc_hard)"
			double nx,ny,nz,f,p,k;
			char c;
			int done = sscanf(spec.substr(3).c_str(), "(%lg,%lg,%lg,%lg,%lg,%lg%c", &nx,&ny,&nz,&f,&p,&k,&c);
			if(done == 7 && c == ')'){
				Eigen::Vector3d n(nx,ny,nz); n.normalize();
				return new LayeredMaterialModel( youngsMod, poissonsRatio,density, strength, compress, toughness, k, 2.0*M_PI*f, p, n);
			}else if(done == 7 && c == '/'){
                // also use varying strength
                double snx,sny,snz,sf,sp,s;
                int done = sscanf(spec.substr( spec.find_first_of('/') ).c_str(), 
                    "/%lg,%lg,%lg,%lg,%lg,%lg%c", &snx,&sny,&snz,&sf,&sp,&s,&c);
                if(done == 7 && c == ')'){
    				Eigen::Vector3d  n(nx,ny,nz);     n.normalize();
                    Eigen::Vector3d sn(snx,sny,snz); sn.normalize();
                    return new LayeredMaterialModel(
                        youngsMod, poissonsRatio,density, strength, compress, toughness,
                        k, 2.0*M_PI*f, p, n,    s, 2.0*M_PI*sf, sp, sn );
                }
            }
		}
		if(spec.compare("tube")==0){
			return new TubedMaterialModel(youngsMod, poissonsRatio,density, strength, compress, toughness, 2.0*toughness);
		}
		if(spec.compare("test")==0){
			return new DummyMaterialModel(youngsMod, poissonsRatio,density, strength,toughness, compress, 3.0);
		}
        if(spec.compare(0,3,"vdb")==0){ // format: "vdb(filename/toughnessGridName/strengthGridName,scaleFactor,voxelSize,useUDF)"
			std::string file, toughname="toughness", strengthname="strength";
			file=spec.substr(4, spec.find_first_of(',')-4);                     //printf("\n1)     file=\"%s\"",file.c_str());
			if(file.find_first_of('/')!=file.npos){
				toughname=file.substr(file.find_first_of('/')+1);               //printf("\n2)    tough=\"%s\"", toughname.c_str());
				file=file.substr(0,file.find_first_of('/'));                    //printf("\n2)     file=\"%s\"",file.c_str());
				strengthname=toughname;
			}
			if(toughname.find_first_of('/')!=toughname.npos){
				strengthname=toughname.substr(toughname.find_first_of('/')+1);  //printf("\n3) strength=\"%s\"", strengthname.c_str());
				toughname=toughname.substr(0,toughname.find_first_of('/'));     //printf("\n3)    tough=\"%s\"",toughname.c_str());
			}                                                                   //printf("\n4) remainder = \"%s\"", spec.substr(spec.find_first_of(',')+1).c_str());
			double scale, voxelSize; int readUDF; char check=0;
			sscanf( spec.substr(spec.find_first_of(',')+1).c_str(), "%lf,%lf,%d%c", &scale, &voxelSize, &readUDF, &check);
			if(check!=')'){
				printf("\nillegal definition of VDB material: %s, use "
					   "\"vdb(filename/toughnessGridName/strengthGridName,scaleFactor,voxelSize,useUDF)\"", spec.c_str());
			}else{
				printf("\n%% ... VDB material model: %s / %s / %s scale %.3lg, voxel size %.3lg, udf %d",
					file.c_str(), toughname.c_str(), strengthname.c_str(), scale, voxelSize, (readUDF!=0));
				return new VDBMaterialModel(
					youngsMod, poissonsRatio,density, strength,toughness, compress,
					file,toughname,strengthname,readUDF!=0, scale ,voxelSize);
			}
        }
        printf("\n%% ... using default material model");
		return new HomogeneousMaterialModel( youngsMod, poissonsRatio,density, strength,toughness, compress);
	}

	// DummyMaterialModel ==================================================================================
	//// "y-split" - doesn't have a proper gradient, testing how badly this affects results
	//double DummyMaterialModel::fractureToughness(const Eigen::Vector3d& x) const{
	//	return Kc* (1.0 + KcDelta* (x[1]>0.0?1.0:0.0));
	//}
	//void DummyMaterialModel::fractureToughnessGradient(Eigen::Vector3d& dKc, const Eigen::Vector3d& x) const{
	//	dKc.setZero(); double h=1e-2;
	//	dKc[1]=2.0*h*(fractureToughness(x+h*Eigen::Vector3d::UnitY())-fractureToughness(x-h*Eigen::Vector3d::UnitY()));
	//}

	// cos-"bubbles"
	double DummyMaterialModel::fractureToughness(const Eigen::Vector3d& x) const{
        const double w=15*M_PI;
		return Kc* (1.0 + KcDelta* (0.5 + 0.5*cos(w*x[0])*cos(w*x[1])*cos(w*x[2])));
	}
	void DummyMaterialModel::fractureToughnessGradient(Eigen::Vector3d& dKc, const Eigen::Vector3d& x) const{
    	const double w=15*M_PI;
		dKc[0]= -Kc*KcDelta* 0.5*w*sin(w*x[0])*cos(w*x[1])*cos(w*x[2]);
		dKc[1]= -Kc*KcDelta* 0.5*w*cos(w*x[0])*sin(w*x[1])*cos(w*x[2]);
		dKc[2]= -Kc*KcDelta* 0.5*w*cos(w*x[0])*cos(w*x[1])*sin(w*x[2]);
	}
	// =====================================================================================================

	// LayeredMaterialModel ==================================================================================
	double LayeredMaterialModel::tensileStrength(const Eigen::Vector3d& x) const {
		double w = 0.5 + 0.5*cos(sf * x.dot(sn) + sp); // w is in [0, 1]
//        printf("\n strength eval at (%+.2lf, %+.2lf, %+.2lf) is %.3lg (w=%.3lf)",
//            x[0],x[1],x[2] , (1.0-w)*Sc + w*ScHard ,w);
		return (1.0-w)*Sc + w*ScHard;
    }
    double LayeredMaterialModel::fractureToughness(const Eigen::Vector3d& x) const{
		double w = 0.5 + 0.5*cos(f * x.dot(n) + p); // w is in [0, 1]
		return (1.0-w)*Kc + w*KcHard;
	}
	void LayeredMaterialModel::fractureToughnessGradient(Eigen::Vector3d& dKc, const Eigen::Vector3d& x) const{
		dKc= (KcHard-Kc)* (-0.5*f)*n*sin( f * x.dot(n) + p );
	}
	// =====================================================================================================

	// TubedMaterialModel ==================================================================================
	double TubedMaterialModel::fractureToughness(const Eigen::Vector3d& x) const{
		double f1=15*M_PI, f2=30*M_PI, p1=0.0, p2=0.0;
		Eigen::Vector3d n1(1.0, 0.0, 0.0), n2(0.0, 0.0, 1.0);
		n1.normalize(); n2.normalize();

		double w = 0.5+0.5*cos(f1*x.dot(n1)+p1)*cos(f2*x.dot(n2)+p2); // w is in [0, 1]
		//printf("\n %.3le ", (1-w)*Kc + w*KcHard);
		return (1-w)*Kc + w*KcHard;
	}
	void TubedMaterialModel::fractureToughnessGradient(Eigen::Vector3d& dKc, const Eigen::Vector3d& x) const{
		double f1=15*M_PI, f2=30*M_PI, p1=0.0, p2=0.0;
		Eigen::Vector3d n1(1.0, 0.0, 0.0), n2(0.0, 0.0, 1.0);
		n1.normalize(); n2.normalize();

		dKc= (KcHard-Kc)*(
			-0.5*f1*n1*sin(f1*x.dot(n1)+p1)*cos(f2*x.dot(n2)+p2)
			-0.5*f2*n2*cos(f1*x.dot(n1)+p1)*sin(f2*x.dot(n2)+p2)
		);
		//printf(" grad( %.3le %.3le %.3le )",dKc[0],dKc[1],dKc[2]);
	}
	// =====================================================================================================

	// VDBMaterialModel ==================================================================================
	class VDBMaterialModel::vdbData{
    public:
        vdbData() : ux(1.0,0.0,0.0), uy(0.0,1.0,0.0), uz(0.0,0.0,1.0) {};
        vdb::FloatGrid::Ptr toughnessGrid, strengthGrid;
        vdb::math::Transform::Ptr xform;
        vdb::Vec3d ux,uy,uz;
		vdb::math::CoordBBox bbox;
    };
    void VDBMaterialModel::init() const{
        haveToughnessMap=false; haveStrengthMap=false;
        if(data==NULL) data = new vdbData();
        vdb::initialize();
        vdb::io::File inFile(vdbFile);
        if(!inFile.open()){
			printf("\n%% Failed to open file %s", vdbFile.c_str());
		}else{
			vdb::GridBase::Ptr baseGrid;
			if( inFile.hasGrid(kName) ){
				baseGrid = inFile.readGrid(kName);
				data->toughnessGrid=vdb::gridPtrCast<vdb::FloatGrid>(baseGrid);
				data->bbox = data->toughnessGrid->evalActiveVoxelBoundingBox();
				haveToughnessMap=true;
			}
			if( kName.compare(sName)) data->strengthGrid = data->toughnessGrid; // share the grid ptr when using the same data for both strength and toughness
			else if( inFile.hasGrid(sName) ){
				baseGrid = inFile.readGrid(sName);
				data->strengthGrid=vdb::gridPtrCast<vdb::FloatGrid>(baseGrid);
				data->bbox.expand( data->strengthGrid->evalActiveVoxelBoundingBox() );
				haveStrengthMap=true;
			}
			inFile.close();
			data->xform = vdb::math::Transform::createLinearTransform(voxelSize);
			//printf("\n%% VDB material model initialized%s%s\n",
			//	haveToughnessMap?", toughness map loaded":"",
			//	haveStrengthMap ? ", strength map loaded":""
			//);
		}
    }
    double VDBMaterialModel::tensileStrength(const Eigen::Vector3d& x) const{
        if(data==NULL) init();
        if(!haveStrengthMap) return Sc;
        vdb::tools::GridSampler<vdb::FloatGrid::ConstAccessor, vdb::tools::BoxSampler> sampler(
            data->strengthGrid->getConstAccessor(), *(data->xform)
        );
		vdb::Vec3d xi = toPeriodicIndexSpace<vdb::Vec3d>(x);
        float gridValue = sampler.isSample(xi) / data->strengthGrid->background();
        if(useUDF) gridValue=std::abs(gridValue);
        return Sc + scale*gridValue;

    }
	double VDBMaterialModel::fractureToughness(const Eigen::Vector3d& x) const{
        if(data==NULL) init();
        if(!haveToughnessMap) return Kc;
        vdb::tools::GridSampler<vdb::FloatGrid::ConstAccessor, vdb::tools::BoxSampler> sampler(
            data->toughnessGrid->getConstAccessor(), *(data->xform)
        );
		vdb::Vec3d xi = toPeriodicIndexSpace<vdb::Vec3d>(x);
        float gridValue = sampler.isSample(xi) / data->toughnessGrid->background();
        if(useUDF) gridValue=std::abs(gridValue);
		//printf("\n eval Kc at (%.3lf,%.3lf,%.3lf) df %.3lf Kc %.3le",x[0],x[1],x[2],gridValue,Kc + scale*gridValue);
        return Kc + scale*gridValue;
    }
	void VDBMaterialModel::fractureToughnessGradient(Eigen::Vector3d& dKc, const Eigen::Vector3d& x) const{
		if(data==NULL) init();
        if(!haveToughnessMap){
            dKc.setZero();
            return;
        }
        vdb::tools::GridSampler<vdb::FloatGrid::ConstAccessor, vdb::tools::BoxSampler> sampler(
            data->toughnessGrid->getConstAccessor(), *(data->xform)
        );
        float x0,x1,y0,y1,z0,z1, bg=data->toughnessGrid->background();
		vdb::Vec3d xi = toPeriodicIndexSpace<vdb::Vec3d>(x);
        x0=sampler.isSample(xi-0.25*data->ux) / bg;
        x1=sampler.isSample(xi+0.25*data->ux) / bg;
        y0=sampler.isSample(xi-0.25*data->uy) / bg;
        y1=sampler.isSample(xi+0.25*data->uy) / bg;
        z0=sampler.isSample(xi-0.25*data->uz) / bg;
        z1=sampler.isSample(xi+0.25*data->uz) / bg;
        if(useUDF){
            x0=std::abs(x0); x1=std::abs(x1);
            y0=std::abs(y0); y1=std::abs(y1);
            z0=std::abs(z0); z1=std::abs(z1);
		}
        dKc[0] = scale*(x1-x0)*2.0/voxelSize;
        dKc[1] = scale*(y1-y0)*2.0/voxelSize;
        dKc[2] = scale*(z1-z0)*2.0/voxelSize;
		//if(dKc.norm() > FLT_EPSILON){
		//	printf("\n x(%.3lf, %.3lf) grad_x %.3le",x0,x1,dKc[0]);
		//	printf("\n y(%.3lf, %.3lf) grad_y %.3le",y0,y1,dKc[1]);
		//	printf("\n z(%.3lf, %.3lf) grad_z %.3le",z0,z1,dKc[2]);
		//}
    }
	template <class Vec3>
	Vec3 VDBMaterialModel::toPeriodicIndexSpace(const Eigen::Vector3d& x) const{
		Vec3 xi = data->xform->worldToIndex(Vec3(x[0],x[1],x[2]));
		//printf("\n pb_map w(%.3lf,%.3lf,%.3lf) -> i(%.1lf,%.1lf,%.1lf)", x[0],x[1],x[2], xi[0],xi[1],xi[2]);
		vdb::math::BBox<Vec3> box ( data->bbox.min().asVec3d(), data->bbox.max().asVec3d() );
		if( box.isInside(xi) ) return xi;
		xi-=box.min(); // shift min to (0,0,0)
		//printf("\n--> s(%.1lf,%.1lf,%.1lf)", xi[0],xi[1],xi[2]);
		xi[0] = std::fmod( xi[0], box.extents()[0] ); if(xi[0]<0.0) xi[0]+=box.extents()[0];
		xi[1] = std::fmod( xi[1], box.extents()[1] ); if(xi[1]<0.0) xi[1]+=box.extents()[1];
		xi[2] = std::fmod( xi[2], box.extents()[2] ); if(xi[2]<0.0) xi[2]+=box.extents()[2];
		//printf("\n--> M(%.1lf,%.1lf,%.1lf), extents (%.1lf,%.1lf,%.1lf)", xi[0],xi[1],xi[2], box.extents()[0],box.extents()[1],box.extents()[2]);
		xi+=box.min(); // undo shift
		//printf("\n--> b(%.1lf,%.1lf,%.1lf) in box (%.0lf,%.0lf,%.0lf)-(%.0lf,%.0lf,%.0lf)", xi[0],xi[1],xi[2], box.min()[0],box.min()[1],box.min()[2], box.max()[0],box.max()[1],box.max()[2]);
		return xi;
	}
	// =====================================================================================================
}
