#ifndef _DIST_PT_TRI_INCLUDED
#define _DIST_PT_TRI_INCLUDED
// note: somewhere in the OpenVDB library there might be a similar function, maybe we can call that instead?

inline double clamp(double x, double a, double b){ return x < a ? a : (x > b ? b : x); }
// compute the squared closest distance in 3D of a point (p) to a triangle with vertices (t_a,t_b,t_c)
// similiar to http://www.gamedev.net/topic/552906-closest-point-on-triangle/
template< typename point_type >
double sqDistPtTri(
	const point_type& p,
	const point_type& t_a, const point_type& t_b, const point_type& t_c,
	double& s, double& t, bool& outsideFlag
){
	Eigen::Vector3d
		e0(t_b[0]-t_a[0],t_b[1]-t_a[1],t_b[2]-t_a[2]),
		e1(t_c[0]-t_a[0],t_c[1]-t_a[1],t_c[2]-t_a[2]),
		v0(t_a[0]-  p[0],t_a[1]-  p[1],t_a[2]-  p[2]);
	double
		a = e0.dot(e0),
		b = e0.dot(e1),
		c = e1.dot(e1),
		d = e0.dot(v0),
		e = e1.dot(v0);
	double det = a*c - b*b;
	s = b*e - c*d;
	t = b*d - a*e;
	
	if ( s + t < det ){
		if ( s < 0.0 ){
			if ( t < 0.0 ){
				if ( d < 0.0 ){
					s = clamp( -d/a, 0.0, 1.0 );
					t = 0.0;
				}else{
					s = 0.0;
					t = clamp( -e/c, 0.0, 1.0 );
				}
			}else{
				s = 0.0;
				t = clamp( -e/c, 0.0, 1.0 );
			}
		}else if ( t < 0.0 ){
			s = clamp( -d/a, 0.0, 1.0 );
			t = 0.0;
		}else{
			s /= det;
			t /= det;
		}
	}else{
		if ( s < 0.0 ){
			if ( (c+e) > (b+d) ){
				s = clamp( ((c+e) - (b+d)) / (a-2*b+c), 0.0, 1.0 );
				t = 1-s;
			}else{
				t = clamp( -e/c, 0.0, 1.0 );
				s = 0.0;
			}
		}else if ( t < 0.0 ){
			if ( a+d > b+e ){
				s = clamp( (c+e-b-d) / (a-2*b+c), 0.0, 1.0 );
				t = 1-s;
			}else{
				s = clamp( -e/c, 0.0, 1.0 );
				t = 0.0;
			}
		}else{
			s = clamp( (c+e-b-d) / (a-2*b+c), 0.0, 1.0 );
			t = 1.0 - s;
		}
	}
	outsideFlag = (e0.cross(e1).dot(v0) < 0); 
	return ((s*e0 + t*e1)+v0).squaredNorm();
}

#endif