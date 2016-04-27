#ifndef QUADRATUREWRAPPER_H
#define	QUADRATUREWRAPPER_H

#include "types.h"

#ifdef _MSC_VER
	#ifndef HYENA_DLL
		#define HYENA_DLL __declspec(dllimport)
	#endif
#else
	// leave empty as GCC exports all symbols by default
	#define HYENA_DLL
#endif

namespace FractureSim{
	/* Returns the number of quadrature points
	 * for the Gaussian quadrature rule of order o on triangles.
	 */
	HYENA_DLL unsigned int getTriangleGaussNumPoints(unsigned int o);
	/* Returns the weight of the n-th quadrature point
	 * for the Gaussian quadrature rule of order o on triangles.
	 * Note sum of weights per rule is 1.
	 */
	HYENA_DLL double getTriangleGaussWeight(unsigned int o, unsigned int n);
	/* Copies the barycentric coordinates of the n-th quadrature point
	 * for the Gaussian quadrature rule of order o on triangles into p.
	 * Returns the input reference for convenience.
	 */
	HYENA_DLL Eigen::Vector3d& getTriangleGaussPoint(unsigned int o, unsigned int n, Eigen::Vector3d& p);
}

#endif