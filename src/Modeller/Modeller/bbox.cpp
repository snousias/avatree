#include "lungModelling.h"




bbox::bbox(void){

	maxY = -10000;
	maxX = -10000;
	maxZ = -10000;
	minY = 10000;
	minX = 10000;
	minZ = 10000;

}



void dotObj::getBoundingBox(void){

	for (int i = 0; i < this->vertices.size(); i++){
		if (this->boundingBox.maxX < this->vertices.at(i).at(0)){ this->boundingBox.maxX = this->vertices.at(i).at(0); }
		if (this->boundingBox.maxY < this->vertices.at(i).at(1)){ this->boundingBox.maxY = this->vertices.at(i).at(1); }
		if (this->boundingBox.maxZ < this->vertices.at(i).at(2)){ this->boundingBox.maxZ = this->vertices.at(i).at(2); }
	}

	for (int i = 0; i < this->vertices.size(); i++){
		if (this->boundingBox.minX > this->vertices.at(i).at(0)){ this->boundingBox.minX = this->vertices.at(i).at(0); }
		if (this->boundingBox.minY > this->vertices.at(i).at(1)){ this->boundingBox.minY = this->vertices.at(i).at(1); }
		if (this->boundingBox.minZ > this->vertices.at(i).at(2)){ this->boundingBox.minZ = this->vertices.at(i).at(2); }
	}

	return;
}