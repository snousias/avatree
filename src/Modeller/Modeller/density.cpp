#include "lungModelling.h"


double dotObj::getDensity(int vertexIndex, float Radius){
	double density = 0.0;
	Vector3f examinedPoint, underInvestigation;

	examinedPoint[0] = this->vertices.at(vertexIndex).at(0);
	examinedPoint[1] = this->vertices.at(vertexIndex).at(1);
	examinedPoint[2] = this->vertices.at(vertexIndex).at(2);

	for (int i = 0; i < this->vertices.size(); i++){
		underInvestigation[0] = this->vertices.at(i).at(0);
		underInvestigation[1] = this->vertices.at(i).at(1);
		underInvestigation[2] = this->vertices.at(i).at(2);

		float r = (examinedPoint - underInvestigation).norm();

		if (r < Radius){
			density += (4.0 / (Pi*powf(Radius, 8)))*powf(((Radius*Radius) - (r*r)), 3);
		}
	}

	return density;
}



double dotObj::getDensity(int vertexIndex, float Radius, std::vector<int> Indices){
	double density = 0.0;
	Vector3f examinedPoint, underInvestigation;

	examinedPoint[0] = this->vertices.at(vertexIndex).at(0);
	examinedPoint[1] = this->vertices.at(vertexIndex).at(1);
	examinedPoint[2] = this->vertices.at(vertexIndex).at(2);

	for (int i = 0; i < Indices.size(); i++){
		underInvestigation[0] = this->vertices.at(Indices.at(i)).at(0);
		underInvestigation[1] = this->vertices.at(Indices.at(i)).at(1);
		underInvestigation[2] = this->vertices.at(Indices.at(i)).at(2);

		float r = (examinedPoint - underInvestigation).norm();
		if (r < Radius){
			density += (4.0 / (Pi*powf(Radius, 8)))*powf(((Radius*Radius) - (r*r)), 3);
		}
	}

	return density;
}

