#include "lungModelling.h"

//----------------Initialize---------------------
void dotObj::initializeSkeletonizationProcess(void)
{
	this->V.setZero(this->vertices.size(), 3);
	this->initialPerVertexRingAreas = this->getRingAreas();
	this->initialAverageFaceArea = this->getAverageFaceArea();
	this->setSL(4);
	this->setWH();
	this->setWL();
	this->adjacencyMatrix = this->getAdjacency();

	return;
}



void dotObj::skeletonize(int iterations){
	this->initializeSkeletonizationProcess();
	for (int i = 0; i < iterations; i++){
		skelStep();
	}
	return;
}


void dotObj::skelStep(){
	Matrix<double, Dynamic, Dynamic> L;
	SparseMatrix<double> Lsp;
	double currentRingArea;
	MatrixXd At;
	At = getRingAreas();
	L = this->getLaplacian("curv");

	//Safety
	for (int i = 0; i < this->vertices.size(); i++){
		if (At(i, i) < 0.00000001)
		{
			return;
		}
	}


	//Lsp = this->getLaplacianSparse("curv");
	//L = MatrixXd(Lsp);
	this->V.setZero();
	convToMat(this->V, this->vertices, 0);
	V = solve(L, WL, WH, V);
	updateWL();
	updateWH();
	convToMat(this->V, this->vertices, 1);
	return;
}


void dotObj::setSL(void){
	this->SL = 3;
	return;
}
void dotObj::setSL(double value){
	this->SL = value;
	return;
}

void dotObj::setWL(void){
	MatrixXd I;
	int size = this->vertices.size();
	I.setIdentity(size, size);
	this->WL = 0.001*sqrt(this->initialAverageFaceArea)*I;
	double x;

	//===Special function for weighting
	if (customFunctionAmplitude >= 0){
		if (uniformContraction){
			std::cout << "Uniform contraction functionallity" << std::endl;
		}
		else{
			std::cout << "Non-Uniform contraction functionallity" << std::endl;
		}
		for (unsigned int i = 0; i < this->geoDistPerVertexFromFurthestVertex.size(); i++){
			if (this->uniformContraction){
				this->WL(i, i) = this->WL(i, i)*(1 - customFunctionAmplitude);
			}
			else{
				x = 2 * Pi*this->customFunctionFrequency*this->geoDistPerVertexFromFurthestVertex.at(i);
				this->WL(i, i) = this->WL(i, i)*(((1 - customFunctionAmplitude)*(0.6 + (0.4*sin(x)))));

				//this->WL(i, i) = (1 - customFunctionAmplitude)*(this->WL(i, i) - customFunctionAmplitude* WL(i, i)*(sin(2 * Pi*this->customFunctionFrequency*this->geoDistPerVertexFromFurthestVertex.at(i))));
				//	this->WL(i, i) = this->WL(i, i) - customFunctionAmplitude*WL(i, i)*this->geoDistPerVertexFromFurthestVertex.at(i); //WL(i, i)*this->geoDistPerVertexFromFurthestVertex.at(i);
			}
		}
	}

	return;
}

void dotObj::updateWL(void)
{
	this->WL = this->SL*this->WL;

	return;
}

void dotObj::setWH(void){
	//Set, initialize WH
	int size = this->vertices.size();
	this->WH0.setIdentity(size, size);
	this->WH.setIdentity(size, size);
	return;
}

void dotObj::updateWH(void)
{
	//Update WH as follows WH=WH0/sqrt(ratio) , ratio = initialPerVertexRingAreas(i, i) / currentRingArea

	double ratio;
	double ratioMax = 100000;
	double currentRingArea;
	MatrixXd At;
	At = getRingAreas();
	for (unsigned int i = 0; i < this->vertices.size(); i++){
		currentRingArea = At(i, i);

		ratio = initialPerVertexRingAreas(i, i) / currentRingArea;
		if (ratio < ratioMax){
			this->WH(i, i) = this->WH0(i, i)*sqrt(ratio); //Division
		}
		else
		{
			this->WH(i, i) = 1;
			this->WL(i, i) = 1;
		}
	}
	return;
}



Matrix<double, Dynamic, Dynamic> dotObj::solve(Matrix<double, Dynamic, Dynamic> & L, Matrix<double, Dynamic, Dynamic> & WL, Matrix<double, Dynamic, Dynamic> & WH, Matrix<double, Dynamic, Dynamic> &V)
{
	MatrixXd WLL = WL*L;
	MatrixXd A(WLL.rows() + WH.rows(), WLL.cols());
	A << WLL, WH;

	//MatrixXd Asquare = A.transpose()*A;

	MatrixXd WHV = WH*V;
	MatrixXd Z(WHV.rows(), WHV.cols());
	Z.setZero();
	//Z.setOnes();
	//Z = Z * 1000;
	MatrixXd B(WHV.rows() + WHV.rows(), WHV.cols());
	B << Z, WHV;



	//MatrixXd Bsquare = A.transpose()*B;

	//return (A.transpose()*A).lu().solve(A.transpose()*B);

	return (A.transpose()*A).ldlt().solve(A.transpose()*B);
	//return (Asquare.transpose()*Asquare).ldlt().solve(Asquare.transpose()*Bsquare);



}
