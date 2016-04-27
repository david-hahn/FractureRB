/* 
 * File:   ElmerReader.h
 * Author: David
 *
 * Created on 19. November 2013, 14:16
 */

#ifndef READER_H
#define	READER_H

#include "types.h"

#include <string>

namespace FractureSim{

    enum ELEM_TYPE {
        NODE = 101,
        LINE = 202,
        TRI = 303,
        TET = 504
    }; // element type constants used in Elmer meshes

    /* Depending on which kind of file we are reading, the element type is
     * found in a different column.
     * The first column is always the ID of the element, and
     * the second column is the body ID to which the element belongs.
     * In .element files the element type is in the 3rd column, whereas
     * in .boundary files, we first have two columns for the parent element and
     * the element type in the 5th column.
     * All columns after the element type are the node IDs forming the element.
     */
    enum TYPE_COLUMN {
        ELEMENTS_FILE = 3,
        BOUNDARY_FILE = 5
    };
    
    /* Provides general reading methods for Elmer mesh files.
     * These files are text based and space-separated with each line
     * representing one node/element in the mesh.
     */
    class ElmerReader {
    public:   
        ElmerReader(std::string nodeFile_="", std::string elemFile_=""){
            nodeFile=nodeFile_;
            elemFile=elemFile_;
        }
        ElmerReader(const ElmerReader& ori){
            nodeFile=ori.nodeFile;
            elemFile=ori.elemFile;
        }
        virtual ~ElmerReader(){}
        
        void setNodeFile(std::string fname){
            nodeFile=fname;
        }
        void setElemFile(std::string fname){
            elemFile=fname;
        }
        
        /* Read elements from the input file
         * The output parameters are elems and bodyId (optional).
         * Previous contents will be deleted.
         * elems contains the node IDs of nodes forming each element
         * bodyIDs contains the body ID of each element
         * 
         * Only elements whose body ID is in the bodies set will be read.
         * If bodies is an empty set, all body IDs will be accepted.
         * Returns the number of elements read successfully or -1 on error
         * 
         * If strictlyListedBodies=false and the bodies set is empty, ALL bodies
         * will be read from the mesh.
         */
        int readElems(
            elem_map& elems, ELEM_TYPE elemType, std::set<unsigned int> bodies,
            int typeColumn=ELEMENTS_FILE, id_map* bodyIDs=NULL,
            elem_map* parentIDs=NULL, bool strictlyListedBodies=false
        );
        
        /* Read nodes from the input file
         * The output parameter is nodes. Previous contents will be deleted.
         * Only elements whose body ID is in the bodies set will be read.
         * Returns the number of nodes read successfully or -1 on error
         */
        int readNodes(node_map& nodes);

    protected:
        std::string nodeFile;
        std::string elemFile;
    };

    /* Provides methods to read a mesh that represents a boundary element model.
     * Hence the main element type (entries in .elements file) is a 
     * surface element. The .boundary file contains line elements
     * (e.g. specifying crack-tips)
     */
    class BEMReader : public ElmerReader{
    public:
        BEMReader() : ElmerReader(){}
        /* Constructs a BEMReader that loads surface elements from
         * meshFile.elements, nodes from meshFile.nodes and boundaries (lines)
         * from meshFile.boundary
         */
        BEMReader(std::string meshFile_) 
            : ElmerReader(meshFile_, meshFile_){
            meshFile=meshFile_;
            nodeFile.append(".nodes");
            elemFile.append(".elements");
        }
        BEMReader(const BEMReader& ori)
            : ElmerReader(ori){
            meshFile=ori.meshFile;
        }
        virtual ~BEMReader(){}
        
        /* Reads elems from .elements file matching elemType and elemBodies
         * (the body ID of each element is stored in bodyIDs),
         * bndrys from .boundary file matching bndryType and bndryBodies and
         * nodes from .nodes file
         */
        int readModel(
            elem_map& elems, id_map& bodyIDs, elem_map& bndrys,
            elem_map& bndryParents, node_map& nodes,
            const ELEM_TYPE elemType,  const std::set<unsigned int> elemBodies,
            const ELEM_TYPE bndryType, const std::set<unsigned int> bndryBodies);
        
    protected:
        std::string meshFile;
    };
    
	class RegionHandler{
	public:
		/* Read region definitions from the given file into regionDefs
		 * the file format is as follows:
		 * comment lines start with a '#' (as the very first character), empty lines are ignored
		 * region definitions start with a line "regions N", where N is the number of region definitions to follow
		 * a region definition has the format n (a,b,c,d) ... where n is the region-ID and
		 * the bracketed expression is a half-space condition for points (x,y,z) of the form ax + by + cz <= d
		 */
		int readRegionDefinitions(std::string filename); //node_map maps id to vector<double>

		/* Assign each element to the first region where at least one of its vertices matches all region-conditions
		 */
		int assignRegions(id_map& regions, node_map& nodes, elem_map& elems);
	protected:
		node_map regionDefs;
	};
 
}

#endif	/* READER_H */

