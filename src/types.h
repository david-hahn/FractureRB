#ifndef TYPES_H
#define	TYPES_H

#include <map>
#include <vector>
#include <set>
#include <Eigen/Dense>
#include <cstdio>
#include <cfloat>

namespace FractureSim{
	enum CONSTANTS{
        NODE_BASE_INDEX = 1, // index origin of nodes
        ELEM_BASE_INDEX = 1
    };
	
	enum FRACTURE_METHOD{
		/*0*/FULL_BEM, // update BEM solution every time step, compute SIFs via displacement correlation from opening displacements
		/*1*/INITIAL_BEM // only compute BEM solution once, ignore fractures, estimate SIFs from the regular stress field
	};
	
    enum CRACK_STATE{ // state of crack-tips
        /*0*/ACTIVE,   // can propagate, if fracture criterion matches
		/*1*/ACTIVE_A, // only node a is valid, b has intersected boundary
		/*2*/ACTIVE_B, // only node b is valid, a has intersected boundary
        /*3*/INACTIVE, // not allowed to propagate for any reason
        /*4*/UPDATE,   // modified ("dirty" state)
        /*5*/UPDATE_A, // modified, was ACTIVE_A before
        /*6*/UPDATE_B  // modified, was ACTIVE_B before
    };

	typedef std::map<unsigned int, double> value_map;
	typedef std::map<unsigned int, std::vector<double> > node_map;
	typedef std::map<unsigned int, std::vector<unsigned int> > elem_map;
	typedef std::map<unsigned int, unsigned int> id_map;
    typedef std::set<unsigned int> id_set;
    typedef std::map<unsigned int, CRACK_STATE> state_map;
	typedef std::map<unsigned int, Eigen::Vector3d> vect3d_map;
    
	typedef Eigen::MatrixXd matrix_type;
	typedef Eigen::VectorXd vector_type;
    
    inline void printMap(elem_map elems){
        for(FractureSim::elem_map::iterator i = elems.begin();
            i != elems.end(); ++i){
            for(int j=0; j<i->second.size(); ++j){
                printf("%u\t ", i->second[j]);
            }
            printf("%% %u\n", i->first);
        }
    }
    
    inline void printMap(node_map nodes){
        for(FractureSim::node_map::iterator i = nodes.begin();
            i != nodes.end(); ++i){
            for(int j=0; j<i->second.size(); ++j){
                printf("%lf\t ", i->second[j]);
            }
            printf("%% %u\n", i->first);
        }
    }

    inline void printMap(vect3d_map nodes){
        for(FractureSim::vect3d_map::iterator i = nodes.begin();
            i != nodes.end(); ++i){
            for(int j=0; j<3; ++j){
                printf("%12.4le\t ", i->second[j]);
            }
            printf("%% %u\n", i->first);
        }
    }

    inline id_set nodeSet(elem_map elems){
        id_set nodes;
        for(elem_map::iterator i = elems.begin(); i != elems.end(); ++i)
            for(elem_map::mapped_type::iterator j = i->second.begin();
                j != i->second.end(); ++j)
                nodes.insert(*j);
        return nodes;
    }

	// NEW FOR FractureRB:
	typedef std::pair<unsigned int, unsigned int> edge;
	typedef std::map<edge, unsigned int> edge_imap;

    inline id_set nodeSet(edge_imap edges){
        id_set nodes;
        for(edge_imap::iterator i = edges.begin(); i != edges.end(); ++i){
			nodes.insert(i->first.first);
			nodes.insert(i->first.second);
		}
        return nodes;
    }
}

#endif
