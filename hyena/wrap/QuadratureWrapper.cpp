#include "QuadratureWrapper.h"
#include "hyena/core/quadrature/gausstria.hpp"

namespace FractureSim{
	using namespace hyena;
	unsigned int getTriangleGaussNumPoints(unsigned int o){
		return gTriaPointsPerOrder[o-1];
	}
	double getTriangleGaussWeight(unsigned int o, unsigned int n){
		return gwTria[gTriaAddress[gTriaPointsPerOrder[o-1]-1]+n];
	}
	Eigen::Vector3d& getTriangleGaussPoint(unsigned int o, unsigned int n, Eigen::Vector3d& p){
		unsigned int pos( (gTriaAddress[gTriaPointsPerOrder[o-1]-1] + n)*3 );
		p[0]=gpTria[pos+0];
		p[1]=gpTria[pos+1];
		p[2]=gpTria[pos+2];
		return p;
	}
}