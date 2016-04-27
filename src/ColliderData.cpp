#include "ColliderData.h"

// make sure we have the proper destructor defs for the collision shapes and mesh interface
#include <btBulletDynamicsCommon.h>
#include <BulletCollision/Gimpact/btGImpactShape.h>

using namespace std;

namespace FractureSim{
	vector<ColliderData> ColliderData::storedData;

	void ColliderData::set(const ColliderData& other){
		shape  = other.shape;
		iface  = other.iface;
		tridata =other.tridata;
		vertdata=other.vertdata;
	}
	void ColliderData::setNull(){
		shape=NULL;
		iface=NULL;
		tridata=NULL;
		vertdata=NULL;
	}
	void ColliderData::deleteAll(bool delShape){
		if(shape && delShape) delete   shape;
		if(iface)    delete   iface;
		if(vertdata) delete[] vertdata;
		if(tridata)  delete[] tridata;
		setNull();
	}
	void ColliderData::storeData(const ColliderData& in){
		storedData.push_back(in);
		//printf("\n store!");
	}
	void ColliderData::clearStoredData(bool delShapes){
		for(int i=0; i<storedData.size(); ++i){
			storedData[i].deleteAll(delShapes);
			//printf("\n deleted stored data");
		}
		storedData.clear();
	}
}
