/*
mySegment.h
defines the function "mySegment":
std::vector<typename GridType::Ptr> mySegment(GridType& grid, double handleThreshold, bool useTiles)

This is used to split a zero-level-set surface stored in grid into separate connected components
running a breath-first-search.

The parameter handleThreshold can be used to remove small handles during segmentation,
it should not exeed the voxel-size of the input grid, otherwise it will cause artifacts.

If useTiles is true, the algorithm runs on interior (grid-value < 0) tiles as well as on voxels
this is required if the object has holes far away from the narrow-band, which would otherwise be detected as
separate components, rather than being kept in the containing segment
otherwise (useTiles=false) the algorithm runs only on the narrow-band, which is obviously faster.

Note that this function is not fully generic
ie. it assumes that the tree structure of the input grid has 4 levels (root, internal1, internal2, leaf)
and that the voxel-size in the input grid is the same in all directions (voxels are cubes).

*/


#ifndef _MY_SEGMENT_INCL_
#define _MY_SEGMENT_INCL_

#include <openvdb/openvdb.h>
#include <openvdb/util/Util.h>
#include <openvdb/tools/Morphology.h>
#include <vector>
#include <deque>

// use #define _DEBUG_SEGMENT_
// modified version from openvdb::tools::internal::segment, can operate on tiles and voxels
// not fully generic, assumes that the tree structure is
// using 4 levels (root, internal1, internal2, leaf)
template<typename GridType>
std::vector<typename GridType::Ptr> mySegment(GridType& grid, double handleThreshold, bool useTiles){
	using namespace openvdb;
    typedef typename GridType::Ptr GridPtr;
    typedef typename GridType::TreeType TreeType;
    typedef typename TreeType::Ptr TreePtr;
    typedef typename TreeType::ValueType ValueType;
	typedef typename TreeType::RootNodeType  RootType;  // level 3 RootNode 
#ifdef _DEBUG_SEGMENT_
	assert(RootType::LEVEL == 3);
#endif
	typedef typename RootType::ChildNodeType Int1Type;  // level 2 InternalNode
	typedef typename Int1Type::ChildNodeType Int2Type;  // level 1 InternalNode
	typedef typename TreeType::LeafNodeType  LeafType;  // level 0 LeafNode
    Coord ijk, n_ijk;
	CoordBBox bndBox;
    ValueType value;
	openvdb::Index depth, n_depth;
	openvdb::Index dim[4]; // size of tile is (dim x dim x dim)
	dim[0]=Int1Type::dim();
	dim[1]=Int2Type::dim();
	dim[2]=LeafType::dim();
	dim[3]=1; //dim[depth==3] is a voxel

#ifdef _DEBUG_SEGMENT_
	unsigned long tvals=0, vvals=0; // count how many tiles we set active, how many ops we do etc.
	unsigned long push_back_counter=0, max_queue_size=0, voxels_done=0, tiles_done=0;
#endif
	for(typename GridType::ValueAllIter it= grid.beginValueAll(); it.test(); ++it){
		if(useTiles || it.isVoxelValue()){ // useTiles==true works also for interior fractures - they are kept as holes in the segment - but this is slower than just using voxels here
			it.setValue(it.getValue()+handleThreshold);
			it.setActiveState( it.getValue() <= 0.0 );
		}
#ifdef _DEBUG_SEGMENT_
		if(it.isTileValue()  && it.isValueOn()) ++tvals;
		if(it.isVoxelValue() && it.isValueOn()) ++vvals;
#endif
	}
#ifdef _DEBUG_SEGMENT_
		printf("\nstarting segmentation with %lu on-voxels and %lu on-tiles\n", vvals, tvals);
#endif

	typename tree::ValueAccessor<TreeType> sourceAcc(grid.tree());
    std::vector<GridPtr> segments;
    std::deque<Coord> coordList;
	try{while (grid.activeVoxelCount() > 0) {
        // Deep copy the grid's metadata (tree and transform are shared)
        GridPtr segment(new GridType(grid, ShallowCopy()));
        // Make the transform unique and insert an empty tree
        segment->setTransform(grid.transform().copy());
        TreePtr tree(new TreeType(grid.background()));
        segment->setTree(tree);

		coordList.clear();
        ijk=grid.beginValueOn().getCoord();
		sourceAcc.setValueOff(ijk);
		coordList.push_back(ijk);
        typename tree::ValueAccessor<TreeType> targetAcc(segment->tree());

        while (!coordList.empty()) {
#ifdef _DEBUG_SEGMENT_
			if(max_queue_size < coordList.size()) max_queue_size = coordList.size();
#endif

            ijk = coordList.back();
            coordList.pop_back();
            if( targetAcc.isValueOn(ijk)) continue;
			value = sourceAcc.getValue(ijk);

			if(sourceAcc.isVoxel(ijk)){
				targetAcc.setValue(ijk, value);
#ifdef _DEBUG_SEGMENT_
				++voxels_done;
#endif

				for (int n = 0; n < 26; n++) {
					n_ijk = ijk + util::COORD_OFFSETS[n];
					if(sourceAcc.probeValue(n_ijk, value)){
						if(sourceAcc.isVoxel(n_ijk)){
							sourceAcc.setValueOff(n_ijk); // switch off as this voxel has been visited
						}else{
							n_depth=sourceAcc.getValueDepth(n_ijk);
							if(n_depth!=-1) //not background value
								sourceAcc.addTile(RootType::LEVEL-n_depth, n_ijk,value,false);

						}
						if(!targetAcc.isValueOn(n_ijk)) {
							coordList.push_back(n_ijk);
#ifdef _DEBUG_SEGMENT_
							++push_back_counter;
#endif
						}
					}else{

						if(value >= 0.0 && value!=grid.background() && !targetAcc.isValueOn(n_ijk)) // simply copy outside voxels
							targetAcc.setValue(n_ijk, value); //++voxels_done; //don't count setting of outside values, we only want to check all active voxels are handled once
					}
				}
			}else{ //printf("\n tile! ");
				// "tiles" are values stored in a node of the tree which is not a leaf
				// to get the correct bounding box of the region represented by a tile
				// touch the leaf containing ijk (generates nodes on all tree-levels including the leaf)
				// get the containing node at the level of the tile
				// get the bounding box of that node
				// replace everything by a tile again (deletes the child branches)
				// note that addTile takes the level NOT the depth of the tile!
				depth = sourceAcc.getValueDepth(ijk); // first get the depth of the original value
				targetAcc.addTile(RootType::LEVEL-depth, ijk,value,true);
#ifdef _DEBUG_SEGMENT_
				++tiles_done;
#endif
				bndBox = CoordBBox::createCube(ijk & ~(dim[depth]-1), dim[depth]);
				//printf("at (%d,%d,%d) has depth %d, lvl %d, bbox size %d, on (%d,%d)",
				//	ijk[0],ijk[1],ijk[2], depth, RootType::LEVEL-depth, bndBox.dim()[bndBox.maxExtent()], sourceAcc.isValueOn(ijk), targetAcc.isValueOn(ijk));

				if(!bndBox.empty()){
                    Coord n_ijks[6];
                    n_ijks[0]=bndBox.min() + Coord(-1, 0, 0);
                    n_ijks[1]=bndBox.min() + Coord( 0,-1, 0);
                    n_ijks[2]=bndBox.min() + Coord( 0, 0,-1);
                    n_ijks[3]=bndBox.max() + Coord( 1, 0, 0);
                    n_ijks[4]=bndBox.max() + Coord( 0, 1, 0);
                    n_ijks[5]=bndBox.max() + Coord( 0, 0, 1);
                    Coord span[12];
                    span[ 0]= Coord(0,1,0); span[ 1]= Coord(0,0,1); // for n_ijks[0]
                    span[ 2]= Coord(1,0,0); span[ 3]= Coord(0,0,1); // for n_ijks[1]
                    span[ 4]= Coord(1,0,0); span[ 5]= Coord(0,1,0); // for n_ijks[2]
                    span[ 6]=Coord(0,-1,0); span[ 7]=Coord(0,0,-1); // for n_ijks[3]
                    span[ 8]=Coord(-1,0,0); span[ 9]=Coord(0,0,-1); // for n_ijks[4]
                    span[10]=Coord(-1,0,0); span[11]=Coord(0,-1,0); // for n_ijks[5]
					for (int n = 0; n < 6; n++) {
						/**/ // simple version, just check the 3 corner-neighbors of the max and min corner of the bounding box (could miss something in special cases but is probably good enough)
						// if we assume that the on-states of the tiles and voxels represents the inside of a level-set and the narrow-band is properly resolved, this should always work.
						if(sourceAcc.probeValue(n_ijks[n],value)){
							n_depth=sourceAcc.getValueDepth(n_ijks[n]);
							if(sourceAcc.isVoxel(n_ijks[n])){
								sourceAcc.setValueOff(n_ijks[n]);
								if( !targetAcc.isValueOn(n_ijks[n]) ){
									coordList.push_back(n_ijks[n]);
#ifdef _DEBUG_SEGMENT_
									++push_back_counter;
#endif
								}
							}else if(n_depth!=-1){ //not background value --> it' a tile
								sourceAcc.addTile(RootType::LEVEL-n_depth, n_ijks[n],value,false);
								if(!targetAcc.isValueOn(n_ijks[n]) ) { // remove this if target tiles are be kept off
									coordList.push_back(n_ijks[n]);
#ifdef _DEBUG_SEGMENT_
									++push_back_counter; //printf("added single coord\n");
#endif
								}
							}
						}
						/*/
						n_depth=sourceAcc.getValueDepth(n_ijks[n]);
						if(n_depth==-1){ //background, ignore
						}else if(n_depth <= depth){ //neighbor is a tile of the same size or larger
							sourceAcc.addTile(RootType::LEVEL-n_depth, n_ijks[n],value,false);
							if(!targetAcc.isValueOn(n_ijks[n]) ) { // remove this if target tiles are be kept off
								coordList.push_back(n_ijks[n]);
#ifdef _DEBUG_SEGMENT_
								++push_back_counter; //printf("added single coord\n");
#endif
								}
						}else{ // neighbors are a bunch of smaller tiles (or voxels)
							//printf("fill face %d size %d (my depth=%d n_depth=%d\n",n,dim[depth],depth,n_depth);
							// stride along directions span[2*n],span[2*n+1]
							// add ALL voxels touching the face of the tile
							for(int n1=0; n1<dim[depth]; ++n1){ //dim[depth]/dim[n_depth]
								for(int n2=0; n2<dim[depth]; ++n2){
									n_ijk = n_ijks[n];
									n_ijk[0]+=n1*span[2*n][0]+n2*span[2*n+1][0];
									n_ijk[1]+=n1*span[2*n][1]+n2*span[2*n+1][1];
									n_ijk[2]+=n1*span[2*n][2]+n2*span[2*n+1][2];
									if(sourceAcc.probeValue(n_ijk,value)){
										n_depth=sourceAcc.getValueDepth(n_ijk);
										if(sourceAcc.isVoxel(n_ijk)){
											sourceAcc.setValueOff(n_ijk);
											if( !targetAcc.isValueOn(n_ijk) ){
												coordList.push_back(n_ijk);
			#ifdef _DEBUG_SEGMENT_
												++push_back_counter;
			#endif
											}
										}else if(n_depth!=-1){ //not background value --> it' a tile
											sourceAcc.addTile(RootType::LEVEL-n_depth, n_ijk,value,false);
											if(!targetAcc.isValueOn(n_ijk) ) { // remove this if target tiles are be kept off
												coordList.push_back(n_ijk);
			#ifdef _DEBUG_SEGMENT_
												++push_back_counter; //printf("added single coord\n");
			#endif
											}
										}
									}
								}
							}
						}/**/
                    }
				}
            }
        }

        //grid.tree().pruneInactive(); // don't do this - it messes up all but the first segment by replacing values with background-tiles before segmentation is complete!!!
		//grid.prune();
		segment->tree().signedFloodFill();
		//segment->prune();

		tools::dilateVoxels(segment->tree(),1+round(handleThreshold/grid.voxelSize()[0])); //assumes that voxels are cubes
		for(FloatGrid::ValueOnIter it = segment->beginValueOn(); it.test(); ++it){
			if( it.isTileValue() ) it.setValueOff();
			if( it.isVoxelValue() && (it.getValue() == segment->background()) ){
				value = sourceAcc.getValue(it.getCoord());
				if( value > handleThreshold )
					it.setValue(value);
			}
			it.setValue(it.getValue()-handleThreshold);
		}


        segments.push_back(segment);
    }}catch(std::exception& e){
		coordList.clear(); segments.clear();
		throw std::logic_error(std::string("in mySegment: ").append(e.what()));
	}
#ifdef _DEBUG_SEGMENT_
	printf("set %lu voxels and %lu tiles\n", voxels_done, tiles_done);
	printf("segments has size %d, done %luK push ops, max queue size was %luK (%lu MB)\n", segments.size(),push_back_counter/1024,(max_queue_size)/1024,(max_queue_size*sizeof(Coord))/1024/1024);
#endif
	return segments;
}
#endif