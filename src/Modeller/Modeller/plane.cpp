#include "lungModelling.h"

template <class T> bool sign(T value){
	if (value > 0){ return true; }
	else{ return false; }
}


bool plane::checkWhichSideOfPlaneIsPoint(Vector3f pointOffPlane){
	bool isCoNormal;
	Vector3f vec ;

	vec = pointOffPlane - this->pointOnPlane;
	isCoNormal = sign<double>(vec.dot(this->normal));
	return isCoNormal;
}

