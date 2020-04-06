#include "lungModelling.h"

Vector3f PCA::getAverage(std::vector<std::vector<float>> &points){
	
	Vector3f average;
	average.setZero();
	
	for (unsigned int i = 0; i < points.size(); i++)
	{
		average[0] += points.at(i).at(0);
		average[1] += points.at(i).at(1);
		average[2] += points.at(i).at(2);
	}


	average[0] = average[0] / (points.size() - 1);
	average[1] = average[1] / (points.size() - 1);
	average[2] = average[2] / (points.size() - 1);
	
	return average;
}

Vector3f PCA::getEigenVector(int order){
	Vector3f theVector;

	return theVector;
}



void PCA::getPlane(int planeOne, int planeTwo){

	return;
}


void PCA::getPCA(void){

	Matrix3f C;
	Vector3f point;
	Vector3f average;
	double maximumValue = 0;
	double middleValue = 0;
	double minimumValue = 0;
	int maximumIndex = 0;
	int middleIndex = 0;
	int minimumIndex = 0;


	C << 0, 0, 0, 0, 0, 0, 0, 0, 0;

	this->centerOfMass = this->getAverage(this->setOfPoints);

	for (unsigned int i = 0; i < this->setOfPoints.size(); i++)
	{
		point[0] = this->setOfPoints.at(i).at(0);
		point[1] = this->setOfPoints.at(i).at(1);
		point[2] = this->setOfPoints.at(i).at(2);

		//Substruct the average each point
		point = point - this->centerOfMass;
		C = C + (point*point.transpose());
	}

	//Divide by number of points
	C = C / (this->setOfPoints.size() - 1);

	//Find Eigen vectors
	SelfAdjointEigenSolver<Matrix3f> es2(C);

	for (int i = 0; i < es2.eigenvalues().size(); i++)
	{
		this->eigenVectors.at(i)[0] = es2.eigenvectors().col(i)[0];
		this->eigenVectors.at(i)[1] = es2.eigenvectors().col(i)[1];
		this->eigenVectors.at(i)[2] = es2.eigenvectors().col(i)[2];
		this->eigenValues.at(i) = es2.eigenvalues()[i];
	}


	//Sort

	for (int i = 0; i < this->eigenValues.size(); i++)
	{
		if (this->eigenValues.at(i) > maximumValue){
			maximumValue = this->eigenValues.at(i);
			maximumIndex = i;
		}
	}

	for (int i = 0; i < this->eigenValues.size(); i++)
	{
		if (i != maximumIndex)
			if (this->eigenValues.at(i)> middleValue){
				middleValue = this->eigenValues.at(i);
				middleIndex = i;
			}
	}

	for (int i = 0; i < this->eigenValues.size(); i++)
	{
		if ((i != maximumIndex) && (i != middleIndex))
		{
			minimumIndex = i;
			minimumValue = this->eigenValues.at(i);
		}
	}

	this->eigenValuesSorted.at(2) = this->eigenValues.at(maximumIndex);
	this->eigenValuesSorted.at(1) = this->eigenValues.at(middleIndex);
	this->eigenValuesSorted.at(0) = this->eigenValues.at(minimumIndex);

	this->eigenVectorsSorted.at(2) = this->eigenVectors.at(maximumIndex);
	this->eigenVectorsSorted.at(1) = this->eigenVectors.at(middleIndex);
	this->eigenVectorsSorted.at(0) = this->eigenVectors.at(minimumIndex);


	return;
}


void PCA::getPCA(std::vector<std::vector<float>> &points){

	Matrix3f C;
	Vector3f point;
	Vector3f average;
	double maximumValue = 0;
	double middleValue = 0;
	double minimumValue = 0;
	int maximumIndex = 0;
	int middleIndex = 0;
	int minimumIndex = 0;


	C << 0, 0, 0, 0, 0, 0, 0, 0, 0;

	this->centerOfMass = this->getAverage(points);

	for (unsigned int i = 0; i < points.size(); i++)
	{
		point[0] = points.at(i).at(0);
		point[1] = points.at(i).at(1);
		point[2] = points.at(i).at(2);

		//Substruct the average each point
		point = point - this->centerOfMass;
		C = C + (point*point.transpose());
	}

	//Divide by number of points
	C = C /(points.size()-1);

	//Find Eigen vectors
	SelfAdjointEigenSolver<Matrix3f> es2(C);
	
	for (int i = 0; i < es2.eigenvalues().size(); i++)
	{
		this->eigenVectors.at(i)[0] = es2.eigenvectors().col(i)[0];
		this->eigenVectors.at(i)[1] = es2.eigenvectors().col(i)[1];
		this->eigenVectors.at(i)[2] = es2.eigenvectors().col(i)[2];
		this->eigenValues.at(i) = es2.eigenvalues()[i];
	}


	//Sort

	for (int i = 0; i < this->eigenValues.size(); i++)
	{
		if (this->eigenValues.at(i) > maximumValue){
			maximumValue = this->eigenValues.at(i);
			maximumIndex = i;
		}
	}

	for (int i = 0; i < this->eigenValues.size(); i++)
	{
		if (i != maximumIndex)
			if (this->eigenValues.at(i)> middleValue){
				middleValue = this->eigenValues.at(i);
				middleIndex = i;
		}
	}

	for (int i = 0; i < this->eigenValues.size(); i++)
	{
		if ((i != maximumIndex) && (i != middleIndex))
		{
			minimumIndex = i;
			minimumValue = this->eigenValues.at(i);
		}
	}
	
	this->eigenValuesSorted.at(2) = this->eigenValues.at(maximumIndex);
	this->eigenValuesSorted.at(1) = this->eigenValues.at(middleIndex);
	this->eigenValuesSorted.at(0) = this->eigenValues.at(minimumIndex);

	this->eigenVectorsSorted.at(2) = this->eigenVectors.at(maximumIndex);
	this->eigenVectorsSorted.at(1) = this->eigenVectors.at(middleIndex);
	this->eigenVectorsSorted.at(0) = this->eigenVectors.at(minimumIndex);


	return;
}

PCA::PCA(void)
{
	Vector3f initVec;
	initVec.setZero();

	//std::cout << "PCA Object is being created" << std::endl;
	this->eigenValues.push_back(0.0);
	this->eigenValues.push_back(0.0);
	this->eigenValues.push_back(0.0);
	this->eigenVectors.push_back(initVec);
	this->eigenVectors.push_back(initVec);
	this->eigenVectors.push_back(initVec);
	this->eigenValuesSorted.push_back(0.0);
	this->eigenValuesSorted.push_back(0.0);
	this->eigenValuesSorted.push_back(0.0);
	this->eigenVectorsSorted.push_back(initVec);
	this->eigenVectorsSorted.push_back(initVec);
	this->eigenVectorsSorted.push_back(initVec);
	centerOfMass = initVec;

}
