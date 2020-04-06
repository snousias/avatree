#include "geodesicDistance.h"

double minimum(double &a, double &b){
	double result;
	if (a < b){ result = a; }
	else{ result = b; }
	return result;
}

double maximum(double &a, double &b){
	double result;
	if (a > b){ result = a; }
	else{ result = b; }
	return result;
}

void minMaxLoc(std::vector<double> &input, std::vector<int> &exclusionList, double &min, double &max, int &minLocation, int &maxLocation){
	min = 1000000;
	max = -1000000;
	std::vector<int>::iterator it;

	for (unsigned int i = 0; i < input.size(); i++){
		it = find(exclusionList.begin(), exclusionList.end(), i);
		if (it == exclusionList.end()){
			if (min > input.at(i)){ min = input.at(i); minLocation = i; }
		}
	}

	//Should be removed for speed purposes only
	
	for (unsigned int count = 0; count < input.size(); count++){
	it = find(exclusionList.begin(), exclusionList.end(), count);
	if (it == exclusionList.end()){
	if (max < input.at(count)){ max = input.at(count); maxLocation = count; }
	}
	}
	

	return;
}

void minMaxLoc(std::vector<double> &input, double &min, double &max, int &minLocation, int &maxLocation){
	min = 1000000;
	max = -1000000;

	for (unsigned int count = 0; count < input.size(); count++){
		if (min > input.at(count)){ min = input.at(count); minLocation = count; }
	}

	for (unsigned int count = 0; count < input.size(); count++){
		if (max < input.at(count)){ max = input.at(count); maxLocation = count; }
	}
	return;
}

double dotObj::getNorm(int p1, int p2){
	/*
	double result;
	Vector3d point1, point2;
	point1(0) = this->vertices.at(p1).at(0);
	point1(1) = this->vertices.at(p1).at(1);
	point1(2) = this->vertices.at(p1).at(2);
	point2(0) = this->vertices.at(p2).at(0);
	point2(1) = this->vertices.at(p2).at(1);
	point2(2) = this->vertices.at(p2).at(2);
	result = (point1 - point2).norm();
	*/

	return  sqrt(powf((this->vertices.at(p2).at(0) - this->vertices.at(p1).at(0)), 2) + powf((this->vertices.at(p2).at(1) - this->vertices.at(p1).at(1)), 2) + powf((this->vertices.at(p2).at(2) - this->vertices.at(p1).at(2)), 2));

	//return result;
}

//====================================================================================

std::vector<double> dotObj::getMinDists(int point1){
	double minVal = 0;
	double maxVal = 0;
	int minLocation = 0;
	int maxLocation = 0;
	unsigned int size;
	int i; int n;
	std::vector<double> weightsMatrix;
	std::vector<int> offMatrix;
	Vector3d p1, p2;
	double result;

	size = this->vertices.size();
	for (unsigned int k = 0; k < size; k++){
		weightsMatrix.push_back(1000.0);
	}
	if (this->adjacencyMatrixSimplified.size() == 0){
		this->adjacencyMatrixSimplified = this->getAdjacencySimplified();
	}
	//std::cout << this->adjacencyMatrixSimplified.size() << std::endl;
	i = point1;
	weightsMatrix.at(i) = 0;
	while (offMatrix.size() < size){
		for (unsigned int j = 1; j < this->adjacencyMatrixSimplified.at(i).size(); j++){
			n = this->adjacencyMatrixSimplified.at(i).at(j);
			result = getNorm(i, n) + weightsMatrix.at(i);
			//result = 1 + weightsMatrix.at(i);
			result = minimum(result, weightsMatrix.at(n));
			weightsMatrix.at(n) = result;
		}
		offMatrix.push_back(i);
		minMaxLoc(weightsMatrix, offMatrix, minVal, maxVal, i, maxLocation);
	}
	return weightsMatrix;
}

std::vector<double> dotObj::getGeodesicDistanceFromSeedPoint(int index){
	double minVal = 0;
	double maxVal = 0;
	int minLocation = 0;
	int maxLocation = 0;
	std::vector<double> W;
	W = this->getMinDists(index);
	minMaxLoc(W, minVal, maxVal, minLocation, maxLocation);
	W = this->getMinDists(maxLocation);
	return W;
}

std::vector<double>  dotObj::geodesicBasedWeighting(double theFreq){
	std::vector<double> W;
	W = this->getGeodesicDistanceFromSeedPoint(10);

	for (unsigned int j = 0; j < W.size(); j++){
		W.at(j) = sin(theFreq*W.at(j));
	}


	return W;
}
std::vector<double>  dotObj::geodesicBasedWeighting(void){
	std::vector<double> W;
	W = this->getGeodesicDistanceFromSeedPoint(10);
	return W;
}

std::vector<int>  dotObj::geodesicSelectorBetweenPoints(int a = -1, int b = -1){
	std::vector<int> list;
	std::vector<double> W1, W2;
	double dist;
	std::vector<int>::iterator it;


	if ((a >= 0) && (b >= 0)){
		W1 = this->getMinDists(a);
		W2 = this->getMinDists(b);
		dist = 1.2*W1.at(b);
		for (int i = 0; i < W1.size(); i++){
			if ((W1.at(i) <= dist) && (W2.at(i) <= dist)){
				it = find(list.begin(), list.end(), i);
				if (it == list.end()) { list.push_back(i); }
			
			}
		}
	}

	return list;
}

std::vector<int>  dotObj::geodesicSelectorBetweenPoints(std::vector<int> pointlist){
	std::vector<int> verticeslist;

	if (pointlist.size() == 1){ verticeslist = pointlist; }
	if (pointlist.size() == 2){ verticeslist = this->geodesicSelectorBetweenPoints(pointlist.at(0), pointlist.at(1)); }
	
	
	/*if (pointlist.size() >2){ 
		verticeslist = this->geodesicSelectorBetweenPoints(pointlist.at(0), pointlist.at(1));
		verticeslist = this->geodesicSelectorBetweenPoints(pointlist.at(1), pointlist.at(2)); 
		verticeslist = this->geodesicSelectorBetweenPoints(pointlist.at(0), pointlist.at(2));
	}*/


	return verticeslist;
}