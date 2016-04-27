// fragments that have less volume than FRAGMENT_SMALL_THR times the original volume are treated as small fragments
// small fragments will get a standard rigid body and will not be breakable anymore
#ifndef FRAGMENT_SMALL_THR
#define FRAGMENT_SMALL_THR (0.02)
#endif
// small fragments get convex hull collision shapes if smaller than FRAGMENT_SMALL_CONV, otherwise meshes are used
// meshes are not working well for very small objects
#ifndef FRAGMENT_SMALL_CONV
#define FRAGMENT_SMALL_CONV (0.005)
#endif
// if a fragment has more volume than FRAGMENT_ORI_THR times the original volume we'll replace
// the parent's implicit surface with the fragment's
// while keeping the BEM mesh etc. intact
#ifndef FRAGMENT_ORI_THR
#define FRAGMENT_ORI_THR (0.95)
#endif

// ignore fragments that have less than FRAGMENT_IGNORE_THR times the original volume
#ifndef FRAGMENT_IGNORE_THR
#define FRAGMENT_IGNORE_THR (1.0e-6)
#endif

// threshold for level-set segmentation in voxel units
#ifndef SEG_HDL_THR
#define SEG_HDL_THR (1.0)
#endif

// control the resolution of collision meshes generated for fragments
// we'll use at most COLL_RES*m triangles where m is the number of surface elements in the BEM mesh
#ifndef COLL_RES
#define COLL_RES (5)
#endif
