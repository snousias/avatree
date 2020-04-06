#include "lungModelling.h"

//=================Skeletonization sparse=======================

void dotObj::skeletonizeSparse(std::string method, int iterations, double lamda){
	//This function performs sparse skeletonization for a number of iteration

	this->initializeSparse();
	for (int i = 0; i < iterations; i++){
		std::cout << "Mesh contraction iteration #" << i << std::endl;
		skelStepSparse(i, "laplacian.horizontal");
	}
	return;
}

void dotObj::skelStepSparse(int iteration, std::string method){
	SparseMatrix<double> LSparse;
	double currentRingArea;
	SparseMatrix<double> VectorAt;
	VectorAt = getRingAreasSparse();
	for (unsigned int i = 0; i < this->vertices.size(); i++)
	{
		currentRingArea = VectorAt.coeffRef(i, i);
		if (currentRingArea < 0.00001){ return; }
	}
	LSparse = this->getLaplacianSparse("curv");
	this->VSparse.setZero();
	convToMatSparse(this->VSparse, this->vertices, 0);
	this->VSparse = solveSparse(LSparse, WLSparse, WHSparse, VSparse);
	updateWLSparse();
	updateWHSparse();
	convToMatSparse(this->VSparse, this->vertices, 1);
	return;
}

SparseMatrix<double> dotObj::getLaplacianSparse(std::string type){
	double rowsum, currentOmega;
	SparseMatrix<double> Output;
	int size = this->vertices.size();
	SparseMatrix<double> L(size, size);
	SparseMatrix<double> Omega(size, size);
	L.setZero();

	if (this->adjacencyMatrixSparse.size() == 0){
		this->adjacencyMatrixSparse = this->getAdjacencySparse();
	}



	for (int i = 0; i < this->adjacencyMatrixSparse.rows(); i++){
		rowsum = 0;
		std::cout << i << std::endl;
		for (int j = 0; j < this->adjacencyMatrixSparse.cols(); j++){
			if (this->adjacencyMatrixSparse.coeffRef(i, j) == 1) {

				L.coeffRef(i, j) = computeOmegaSparse(i, j);
				rowsum += L.coeffRef(i, j);
			}
		}
		L.coeffRef(i, i) = -rowsum;
	}
	return L;
}

SparseMatrix<double> dotObj::getRingAreasSparse(void){
	SparseMatrix<double> getRingVec(this->vertices.size(), this->vertices.size());

	getRingVec.setZero();

	for (unsigned int i = 0; i < this->vertices.size(); i++){
		getRingVec.coeffRef(i, i) = getRingArea(i);
	}

	return getRingVec;
}

void dotObj::initializeSparse(void){
	this->VSparse = Eigen::SparseMatrix<double>(this->vertices.size(), 3);
	this->initialPerVertexRingAreasSparse = this->getRingAreasSparse();
	this->initialAverageFaceArea = this->getAverageFaceArea();
	this->setSL();
	this->setWHSparse();
	this->setWLSparse();
	this->adjacencyMatrixSparse = this->getAdjacencySparse();
	return;
}

void dotObj::setWHSparse(void){
	int size = this->vertices.size();
	this->WH0Sparse = SparseMatrix <double>(size, size);
	this->WH0Sparse.setIdentity();
	this->WHSparse = SparseMatrix <double>(size, size);
	this->WHSparse.setIdentity();
	return;
}

void dotObj::setWLSparse(void){
	int size = this->vertices.size();
	SparseMatrix<double> I(size, size);
	I.setIdentity();

	this->WLSparse = 0.001*sqrt(this->initialAverageFaceArea)*I;

	//===Special function for weighting
	/*
	#ifdef WEIGHTER
	for (unsigned int i = 0; i < this->geoDistPerVertexFromFurthestVertex.size(); i++){
	this->WL(i, i) = this->WL(i, i) - amplitude* WL(i, i)*this->geoDistPerVertexFromFurthestVertex.at(i);
	}
	#endif
	*/

	return;
}

void dotObj::updateWLSparse(void)
{
	this->WLSparse = this->SL*this->WLSparse;

	return;
}

void dotObj::updateWHSparse(void)
{
	double ratio;
	double ratioMax = 100000;
	double currentRingArea;
	SparseMatrix <double> At;
	At = getRingAreasSparse();
	for (unsigned int i = 0; i < this->vertices.size(); i++){
		currentRingArea = At.coeffRef(i, i);

		ratio = initialPerVertexRingAreasSparse.coeffRef(i, i) / currentRingArea;
		if (ratio < ratioMax){
			this->WHSparse.coeffRef(i, i) = this->WH0Sparse.coeffRef(i, i)*sqrt(ratio); //Division
		}
		else
		{
			this->WHSparse.coeffRef(i, i) = 1;
			this->WLSparse.coeffRef(i, i) = 1;
		}
	}
	return;
}

SparseMatrix<double> dotObj::solveSparse(SparseMatrix<double> &LSparse, SparseMatrix<double>  &WLSparse, SparseMatrix<double>  &WHSparse, SparseMatrix<double>  &VSparse)
{
	std::cout << "Solve linear system at time" << std::endl;


	//this->toc();

	SparseMatrix<double> A((LSparse.rows() + WHSparse.rows()), LSparse.cols());
	SparseMatrix<double> AUpper(LSparse.rows(), LSparse.cols());
	SparseMatrix<double> ALower(WHSparse.rows(), WHSparse.cols());
	SparseMatrix<double> B((LSparse.rows() + WHSparse.rows()), VSparse.cols());
	SparseMatrix<double> BUpper(LSparse.rows(), VSparse.cols());
	SparseMatrix<double> BLower(WHSparse.rows(), VSparse.cols());
	SparseMatrix<double> X(VSparse.rows(), VSparse.cols());
	X.setZero();
	//Z.setZero();

	AUpper = WLSparse*LSparse;
	ALower = WHSparse;
	BUpper.setZero();
	BLower = WHSparse*VSparse;



	for (int j = 0; j < A.cols(); j++){
		for (int i = 0; i < AUpper.rows(); i++){
			A.coeffRef(i, j) = AUpper.coeffRef(i, j);
		}

		for (int i = 0; i < ALower.rows(); i++){
			A.coeffRef(i + AUpper.rows() - 1, j) = ALower.coeffRef(i, j);
		}
	}

	for (int j = 0; j < B.cols(); j++){
		for (int i = 0; i < BUpper.rows(); i++){
			B.coeffRef(i, j) = BUpper.coeffRef(i, j);
		}

		for (int i = 0; i < BLower.rows(); i++){
			B.coeffRef(i + BUpper.rows() - 1, j) = BLower.coeffRef(i, j);
		}
	}

	std::cout << "Solve linear system at time" << std::endl;

	SparseMatrix<double> Asq = A.transpose()*A;
	SparseMatrix<double> Bsq = A.transpose()*B;
	SparseLU  <SparseMatrix<double>> solver;
	//SparseQR<SparseMatrix<double>, AMDOrdering<int>> solver;
	//Asq.makeCompressed();
	solver.compute(Asq);
	X = solver.solve(Bsq);


	return X;
}

double dotObj::computeOmegaSparse(int firstIndex, int secondIndex){
	double omega = 0;
	std::vector<int> N;
	Vector3d firstPoint, secondPoint, NeighbourPoint, A, B;

	N = findCommonNeighboursSparse(firstIndex, secondIndex);

	for (int i = 0; i < 3; i++){
		firstPoint[i] = this->vertices.at(firstIndex).at(i);
		secondPoint[i] = this->vertices.at(secondIndex).at(i);
	}

	//std::cout << N.size() << std::endl;

	if (N.size() == 2){
		for (unsigned int j = 0; j < N.size(); j++){
			for (int i = 0; i < 3; i++){
				NeighbourPoint[i] = this->vertices.at(N.at(j)).at(i);
			}

			A = firstPoint - NeighbourPoint;
			B = secondPoint - NeighbourPoint;
			//omega += (1 / tan(computeAngle(firstPoint - NeighbourPoint, secondPoint - NeighbourPoint)));
			omega += A.dot(B) / (A.cross(B)).norm();  //Division
		}
	}
	else
	{
		omega = 1;
	}

	//omega = (1 / tan(computeAngle(firstPoint - NeighbourPoint, secondPoint - NeighbourPoint))) + (1 / tan(computeAngle(firstPoint - secondNeighbourPoint, secondPoint - secondNeighbourPoint)));

	return omega;
}

std::vector<int> dotObj::findCommonNeighboursSparse(int p1, int p2){
	std::vector<int> N;
	int size = this->vertices.size();
	for (int i = 0; i < size; i++){
		if ((this->adjacencyMatrixSparse.coeffRef(p1, i) == 1) && (this->adjacencyMatrixSparse.coeffRef(p2, i) == 1) && (i != p1) && (i != p2)){
			N.push_back(i);
		}
	}
	return N;
}

void dotObj::convToMatSparse(SparseMatrix<double> &V, std::vector<std::vector<float>>&vert, int direction){
	if (direction == 0){
		for (unsigned int i = 0; i < this->vertices.size(); i++){
			for (int j = 0; j < 3; j++){							//this->vertices.at(0).size()
				V.coeffRef(i, j) = this->vertices.at(i).at(j);
			}
		}
	}
	else{
		for (unsigned int i = 0; i < this->vertices.size(); i++){
			for (int j = 0; j < 3; j++){
				this->vertices.at(i).at(j) = V.coeffRef(i, j);
			}
		}
	}
	return;
}