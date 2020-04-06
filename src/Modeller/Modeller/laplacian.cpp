#include "lungModelling.h"



double dotObj::computeOmega(int firstIndex, int secondIndex){
	double omega = 0;
	std::vector<int> N;
	Vector3d firstPoint, secondPoint, NeighbourPoint, A, B;
	N = findCommonNeighbours(firstIndex, secondIndex);
	for (int i = 0; i < 3; i++){
		firstPoint[i] = this->vertices.at(firstIndex).at(i);
		secondPoint[i] = this->vertices.at(secondIndex).at(i);
	}
	if (N.size() == 2){
		for (unsigned int j = 0; j < N.size(); j++){
			for (int i = 0; i < 3; i++){
				NeighbourPoint[i] = this->vertices.at(N.at(j)).at(i);
			}

			A = firstPoint - NeighbourPoint;
			B = secondPoint - NeighbourPoint;
			//omega += (1 / tan(computeAngle(firstPoint - NeighbourPoint, secondPoint - NeighbourPoint)));
			omega += A.dot(B) / (A.cross(B)).norm();  //Division

			//if (omega > 10000)omega = 0;//Check


		}
	}
	else
	{
		omega = 1;
	}
	//omega = (1 / tan(computeAngle(firstPoint - NeighbourPoint, secondPoint - NeighbourPoint))) + (1 / tan(computeAngle(firstPoint - secondNeighbourPoint, secondPoint - secondNeighbourPoint)));
	return omega;
}



Matrix<double, Dynamic, Dynamic> dotObj::getLaplacian(std::string type){
	//Get laplacian matrix
	Matrix <double, Dynamic, Dynamic > Output;
	Matrix < double, Dynamic, Dynamic > L;
	Matrix < double, Dynamic, Dynamic > curvatureFlowL;
	int size = this->vertices.size();
	/*--------------------------L------------------------------*/
	Matrix<double, Dynamic, Dynamic>I, A, Coef, Omega;
	VectorXd D;
	I.setIdentity(size, size);

	if (this->adjacencyMatrix.size() == 0){
		std::cout << "Get Adjacency" << std::endl;
		this->adjacencyMatrix = this->getAdjacency();
	}

	if (type == "onoff"){
		A = this->adjacencyMatrix;
		D = A.rowwise().sum();
		//D = A*VectorXd::Ones(size);
		D = D.array().inverse().transpose();
		A = D.asDiagonal()*A;
		L = I - A;
	}
	/*--------------------Curvature-Flow-L------------------*/
	if (type == "curv"){
		A = this->adjacencyMatrix;
		Omega.setZero(size, size);

		for (int i = 0; i < size; i++){
			for (int j = 0; j < size; j++){
				if (this->adjacencyMatrix(i, j) == 1) {
					Omega(i, j) = computeOmega(i, j);
				}
			}
		}

		D = Omega.rowwise().sum();
		Coef = D.transpose().asDiagonal()*I;
		curvatureFlowL = -Coef + Omega;
	}

	if (type == "onoff"){
		Output = L;
	}

	if (type == "curv"){
		Output = curvatureFlowL;
	}

	return Output;
}
