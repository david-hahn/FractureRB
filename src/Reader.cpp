/* 
 * File:   ElmerReader.cpp
 * Author: David
 * 
 * Created on 19. November 2013, 14:16
 */

#include "Reader.h"

#include <fstream>
#include <sstream>
#include <cstdio>

using namespace std;

namespace FractureSim{
	extern int readVCG(std::string filename, node_map& nodes, elem_map& elems); // prototype of helper function in Reader_VCG.cpp

    int ElmerReader::readElems(elem_map& elems,
            ELEM_TYPE elemType, std::set<unsigned int> bodies,
            int typeColumn, id_map* bodyIDs, elem_map* parentIDs, bool strictlyListedBodies
    ){
        unsigned int elemId,bodyId,buffer,column;
        bool flag;
        istringstream strstream;
        string line, token;
        std::vector<unsigned int> nodes; // store node IDs
        std::vector<unsigned int> parents; // store parents when reading boundary elements
        elems.clear();
        if(bodyIDs!=NULL) bodyIDs->clear();

        ifstream in(elemFile.c_str());
        if(!in.is_open()) return -1;
        while(in.good()){
            getline(in, line);
            strstream.clear();
            strstream.str(line);
            nodes.clear();
			parents.clear();
//            printf("parsing line \"%s\"\n", line.c_str());
//            printf("strstream holds \"%s\"\n", strstream.str().c_str());
            // first column is element ID
            getline(strstream, token, ' ');
            sscanf(token.c_str(), "%u", &elemId);
//            printf(" token 1 is \"%s\"\n",token.c_str());
            // second column is body id
            getline(strstream, token, ' ');
            sscanf(token.c_str(), "%u", &bodyId);
//            printf(" token 2 is \"%s\"\n",token.c_str());
//            printf("body id is %u", bodyId);
            flag= (!strictlyListedBodies && bodies.empty()) || (bodies.count(bodyId)>0); // only read further if body ID matches
//            printf(" -- %s found\n", flag?"":"not");
            column=3;
            while(getline(strstream, token, ' ') && flag){
                sscanf(token.c_str(), "%u", &buffer);
//                printf(" token %u is \"%s\"\n",column, token.c_str());
				if((parentIDs!=NULL) && (column < typeColumn)){ // columns 3 and 4 in a boundary file are the parent elements of the boundary element (face)
					parents.push_back(buffer);
				}else if(column == typeColumn){
                    // check the element type
                    flag = (buffer==elemType); // stop reading if the elment type is wrong
                }else if(column > typeColumn){
                    // store node ID
                    nodes.push_back(buffer);
                }
                ++column;
            }
            if(flag && !nodes.empty()){
                elems[elemId]=nodes;
                if(bodyIDs!=NULL) (*bodyIDs)[elemId]=bodyId;
				if(parentIDs!=NULL) (*parentIDs)[elemId]=parents;
//                printf("added element %u\n", elemId);
            }
        }
        in.close();
        return elems.size();
    }
    
    int ElmerReader::readNodes(node_map& nodes){
        string line;
        double coords[3]; // store node coords
        int id, buffer, tokens;
        char test;
        nodes.clear();
		
//		printf("reading nodes from \"%s\"\n",nodeFile.c_str());
        ifstream in(nodeFile.c_str());
        if(!in.is_open()) return -1;
//		printf("file opened successfully\n");
        while(in.good()){
            getline(in, line); test=0;
            // format for nodes is "node-id unused x y z"
            tokens=sscanf(line.c_str(),"%d %d %lf %lf %lf%c",
                    &id, &buffer, coords, coords+1, coords+2, &test);
            //printf("read %d tokens: %d %d %lf %lf %lf, test is %d\n",
            //        tokens, id, buffer, coords[0], coords[1], coords[2], test);
            if(tokens==5 && test==0){ //line format is correct
                nodes[id]=vector<double>(coords, coords+3);
            }
        }
        in.close();
//		printf("read %d nodes from %s\n", nodes.size(),nodeFile.c_str());
        return nodes.size();
    }
    
    
    int BEMReader::readModel(elem_map& elems, id_map& bodyIDs,
        elem_map& bndrys, elem_map& bndryParents, node_map& nodes,
        const ELEM_TYPE elemType,  const std::set<unsigned int> elemBodies,
        const ELEM_TYPE bndryType, const std::set<unsigned int> bndryBodies
    ){
		//if meshFile ends with ".stl", ".obj", ".ply" etc. use VCG to read raw-mesh (nodes & elems, nothing more)
		if( meshFile.find(".stl")!=meshFile.npos ||
			meshFile.find(".obj")!=meshFile.npos ||
			meshFile.find(".off")!=meshFile.npos ||
			meshFile.find(".ply")!=meshFile.npos
		){
			printf("\n%% reading raw triangle mesh from %s (requires --remesh option)\n%%\t",meshFile.c_str());
			return readVCG(meshFile,nodes,elems);
		}
		node_map nodes_in;
        //elem_map elems_in;
        int ret;
        // read elements and nodes
		ret=readElems(elems, elemType, elemBodies,ELEMENTS_FILE,&bodyIDs);
        if(ret<0) return ret;
        ret=readNodes(nodes_in);
        if(ret<0) return ret;
        // now switch to .boundary file
        string tmp=elemFile;
        elemFile=meshFile;
        elemFile.append(".boundary");
        ret=readElems(bndrys,bndryType,bndryBodies,BOUNDARY_FILE,NULL,&bndryParents,true);
        // restore state of the BEMReader object
        elemFile=tmp;
        //if(ret<0) return ret; // it's ok to not have a boundary file, if there's no pre-defined cracks

		// reading is complete, but for HyENA-BEM compatibility the nodes need to be numbered
		// in the same order as they appear in the element map.
        
		// run through the element map and decide new node numbers
		id_map fwd;// bkwd; // fwd[old_id]=new_id, bkwd[new_id]=old_id
		unsigned int new_id=NODE_BASE_INDEX;
		for(elem_map::iterator i = elems.begin(); i!=elems.end(); ++i){
			//run through all nodes of the element
			for(elem_map::mapped_type::iterator j = i->second.begin(); j!=i->second.end(); ++j){
				if(fwd.count(*j)==0){ // assing new number at first occurence of node
					fwd[*j]=new_id; //bkwd[new_id]=*j;
					new_id++;
				}
				(*j)=fwd[*j]; //update element
			}
		}
        nodes.clear();
		// copy from nodes_in to nodes while applying new numbering
		for(node_map::iterator i = nodes_in.begin(); i!= nodes_in.end(); ++i){
			nodes[fwd[i->first]] = i->second;
		}
		// apply new numbering to bndry elements
		for(elem_map::iterator i = bndrys.begin(); i!=bndrys.end(); ++i){
			for(elem_map::mapped_type::iterator j = i->second.begin(); j!=i->second.end(); ++j){
				(*j)=fwd[*j]; //update element
			}
		}

		return elems.size();
    }

	int RegionHandler::assignRegions(id_map& regions, node_map& nodes, elem_map& elems){
		regions.clear();
		Eigen::Vector3d n; // tmp storage for a region-constraint
		double d; // such that the constraint is evaluated on p as n.dot(p)<=d
		bool found, match;

		unsigned int i=0;
		for(elem_map::const_iterator it=elems.begin(); it!=elems.end(); ++it,++i){
			// check which region fits this tri
			Eigen::Vector3d // a,b,c are node coordinates
				a (nodes[it->second[0]][0], nodes[it->second[0]][1], nodes[it->second[0]][2]),
				b (nodes[it->second[1]][0], nodes[it->second[1]][1], nodes[it->second[1]][2]),
				c (nodes[it->second[2]][0], nodes[it->second[2]][1], nodes[it->second[2]][2]);
			found=false;
			for(node_map::iterator rd=regionDefs.begin(); rd!=regionDefs.end() && !found; ++rd){ //iterate region defs
				match=true;
				for(int j=0; j<rd->second.size() && match; j+=4){
					n[0]=rd->second[j  ]; n[1]=rd->second[j+1]; n[2]=rd->second[j+2];
					d   =rd->second[j+3];
					if( n.dot(a) > d || // using && means tri is added to region if at least 1 vertex matches
						n.dot(b) > d || // using || means tri is added to region if all of its vertices match
						n.dot(c) > d ) match=false;
				}
				if(match){
					found=true;
					regions[it->first]=rd->first;
				}
			}
		}
		return 0;
	}

	int RegionHandler::readRegionDefinitions(string filename){ //node_map maps id to vector<double>
		int nRegions=0; // number of regionDefs read
		bool regKwdFound=false; // look for the keyword "regions" followed by the number of regions to read
		string line, token;
		int ret, done=0;
		char check;
		double a,b,c,d;
		unsigned int nextId;

		regionDefs.clear();

		std::istringstream strstream;
		ifstream in(filename.c_str());
		if(!in.is_open()) return -1;

		while(in.good()){
			getline(in, line);
			if(line.empty() || line.at(0)=='#') continue; // comments have # in the first column
			strstream.clear();
			strstream.str(line);

			getline(strstream, token, ' ');
			if(!regKwdFound){
				if(token.compare("regions")==0){
					regKwdFound=true;
					getline(strstream, token, ' '); //next token is number of regions to read
					ret=sscanf(token.c_str(), "%u", &nRegions);
					if(ret!=1){
						in.close();
						//printf("invalid ret=%d should be 1 on token %s\n",ret, token.c_str());
						return -1;
					}
				}
			}else if(done<nRegions){ //reading next region until we have plenty
				//printf("processing '%s' done %d\n",line.c_str(),done);
				ret=sscanf(token.c_str(), "%u", &nextId);
				if(ret==1) while(getline(strstream, token, ' ')){ //token is now one condition of the region definition
					//printf("parsing token '%s'",token.c_str());
					ret=sscanf(token.c_str(), "(%lf,%lf,%lf,%lf%c",&a,&b,&c,&d,&check);
					if(ret==5 && check==')'){ // correct format
						regionDefs[nextId].push_back(a);
						regionDefs[nextId].push_back(b);
						regionDefs[nextId].push_back(c);
						regionDefs[nextId].push_back(d);
						//printf(" ... ok\n");
					}
					//else printf(" ... reject ret=%d, check='%c'\n",ret,check);
				}
				++done;
			}
		}
		//printf("\nregionDefs:\n");
		//printMap(regionDefs);
		// finally add an empty region def which will be the default
		if(regionDefs.empty()) nextId=ELEM_BASE_INDEX;
		else nextId = regionDefs.rbegin()->first +1;
		regionDefs[nextId].assign(0,0.0);

		return nRegions;
	}
}
