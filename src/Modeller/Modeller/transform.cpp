#include "transform.h"


using namespace Eigen;
void dotObj::removeVertex(int i){

	this->vertices.erase(this->vertices.begin() + i);

}


void dotObj::removeNormal(int i){



	this->normals.erase(this->normals.begin() + i);


	return;
}

void dotObj::removeFace(int i){

	this->faces.erase(this->faces.begin() + i);

}


void dotObj::rotate(Vector3f & target){
	
	Vector3f point;

	//Normalize target up vector
	target.normalize();

	//find obj up vector...
	Vector3f up = Vector3f(0.0, 1.0, 0.0);
	
	//find vector vertical to both vectors
	Vector3f axis = up.cross(target);
	axis.normalize();
	//fin angle
	float angle = up.dot(target)/(up.norm()*target.norm());
	angle = acos(angle);
	//std::cout << angle << std::endl;
	//std::cout << axis << std::endl;
	//Declare rotation
	Quaternion<float> q;  
	q = AngleAxis<float>(angle, axis);
	
	for (int i = 0; i < this->vertices.size();i++){
		point = Vector3f(this->vertices.at(i).at(0), this->vertices.at(i).at(1), this->vertices.at(i).at(2));
		point = q*point;
		this->vertices.at(i).at(0) = point[0];
		this->vertices.at(i).at(1) = point[1];
		this->vertices.at(i).at(2) = point[2];
	
	}
	
	return;
}


void dotObj::scale(double scaleX, double scaleY, double scaleZ){
	for (int i = 0; i < this->vertices.size(); i++){
		this->vertices.at(i).at(0) = scaleX*this->vertices.at(i).at(0);
		this->vertices.at(i).at(1) = scaleY*this->vertices.at(i).at(1);
		this->vertices.at(i).at(2) = scaleZ*this->vertices.at(i).at(2);

	}


}

void dotObj::translate(Vector3f & target){


	//Find center of mass
	this->findCenterOfMass();
	

	//Diff center to target
	Vector3f diff = target - this->centerOfMass;

	

	

	//Add target to x,y,z
	for (int i = 0; i < this->vertices.size(); i++){
		this->vertices.at(i).at(0) = this->vertices.at(i).at(0) + diff[0];
		this->vertices.at(i).at(1) = this->vertices.at(i).at(1) + diff[1];
		this->vertices.at(i).at(2) = this->vertices.at(i).at(2) + diff[2];

	}

	return;
}

Vector3f dotObj::findCenterOfMass(void){
	
	

	float c0=0;
	float c1 = 0;
	float c2 = 0;

	for (int i = 0; i < this->vertices.size(); i++){
		c0+=this->vertices.at(i).at(0) ;
		c1+=this->vertices.at(i).at(1) ;
		c2+=this->vertices.at(i).at(2) ;

	}
	c0 = c1 / this->vertices.size();
	c1 = c1 / this->vertices.size();
	c2 = c1 / this->vertices.size();

	this->centerOfMass = Vector3f(c0, c1, c2);

	return Vector3f(c0, c1, c2);
}
