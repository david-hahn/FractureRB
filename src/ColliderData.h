#ifndef COLLIDERDATA_H
#define	COLLIDERDATA_H

#include <vector>

//fwd decl.
class btCollisionShape;
class btStridingMeshInterface;

namespace FractureSim{
	// a little helper class to keep track of collision meshes that we need to delete eventually
	class ColliderData{
	public:
		ColliderData(){ setNull(); }
		ColliderData(const ColliderData& ori){ set(ori); }
		~ColliderData(){}
		btCollisionShape* shape;
		btStridingMeshInterface* iface;
		char* tridata;
		char* vertdata;
		void set(const ColliderData& other);
		void setNull(); // sets stored pointers to NULL without deleting objects
		void deleteAll(bool delShape=true); // deletes all referenced objects, then calls setNull(),if delShape=false the shape object is assumed to have been deleted already
		static void clearStoredData(bool delShapes=true); // if delShapes=false the shape objects are assumed to have been deleted already
		static void storeData(const ColliderData& in);
	protected:
		static std::vector<ColliderData> storedData;
	};
}

#endif //COLLIDERDATA_H
