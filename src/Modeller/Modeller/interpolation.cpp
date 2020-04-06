#include "lungModelling.h"


void dotObj::interp(dotObj initial, double percentage){
	std::cout << "Interpolation" << std::endl;

	Vector3d diff, p, ip;
	for (int i = 0; i < this->vertices.size(); i++){
		diff[0] = this->vertices.at(i).at(0) - initial.vertices.at(i).at(0);
		diff[1] = this->vertices.at(i).at(1) - initial.vertices.at(i).at(1);
		diff[2] = this->vertices.at(i).at(2) - initial.vertices.at(i).at(2);
		diff = (1 - percentage)*diff;

		p[0] = initial.vertices.at(i).at(0);
		p[1] = initial.vertices.at(i).at(1);
		p[2] = initial.vertices.at(i).at(2);

		p = p + diff;

		this->vertices.at(i).at(0) = p[0];
		this->vertices.at(i).at(1) = p[1];
		this->vertices.at(i).at(2) = p[2];
	}

	return;
}

