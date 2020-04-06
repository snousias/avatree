#include "obj.h"


void dotObj::mapNormals2Faces(void){
	std::cout << "Commence Mapping normals to corresponding vertices" << std::endl;
	//Find corresponding normal to each vertices using faces and store
	std::cout << "Initializing vertexToNormalMatch" << std::endl;

	for (unsigned int i = 0; i < this->faces.size(); i++){
		this->vertexToNormalMatch.push_back({ (float)i, 0 });
	}

	std::cout << "Commence Processing faces" << std::endl;
	for (unsigned int i = 0; i < this->faces.size(); i++){
		this->vertexToNormalMatch.at((int)this->faces.at(i).at(0) - 1).at(1) = this->faces.at(i).at(2) - 1;
		this->vertexToNormalMatch.at((int)this->faces.at(i).at(3) - 1).at(1) = this->faces.at(i).at(5) - 1;
		this->vertexToNormalMatch.at((int)this->faces.at(i).at(6) - 1).at(1) = this->faces.at(i).at(8) - 1;
	}
	std::cout << "Processing faces complete" << std::endl;
	std::cout << "Mapping Faces Complete" << std::endl;
}

void  dotObj::registerPart(dotObj &input,bool doNormals){ //child.registerPart(Parent)
	int k = 0;
	int buffer = 0;
	int toBeInsertedIntoFaces[9];
	bool tripleCheck[3];
	std::vector<int>::iterator itFinder;

	//GENERATE SELECTED VERTICES TABLE
	for (int i = 0; i < input.vertices.size(); i++){
		this->assemblyToChildPartVerticesMatch.push_back({ 0, 0 });
		itFinder = find(input.selectedVertices.begin(), input.selectedVertices.end(), i);
		if (itFinder != input.selectedVertices.end()){
			this->selectedVertices.push_back(i);
			this->vertices.push_back(input.vertices.at(i));
			this->assemblyToChildPartVerticesMatch.at(i).at(0) = i;
			this->assemblyToChildPartVerticesMatch.at(i).at(1) = k;
			k++;
		}
	}
	std::cout << "Register part:" << this->assemblyToChildPartVerticesMatch.size() << std::endl;

	//----------------------------------------------------------------------------------
	//GENERATE NORMALS
	if (doNormals){
		this->normals = input.normals;
	}
	//----------------------------------------------------------------------------------
	for (unsigned int i = 0; i < input.faces.size(); i++)
	{
		tripleCheck[0] = false; tripleCheck[1] = false; tripleCheck[2] = false;

		//FACE #0
		buffer = (int)input.faces.at(i).at(0) - 1;

		for (unsigned int iCheck = 0; iCheck < this->selectedVertices.size(); iCheck++)
		{
			if (this->selectedVertices.at(iCheck) == buffer)
			{
				tripleCheck[0] = true;
			}
		}

		if (tripleCheck[0] == true){
			//FACE #1
			buffer = (int)input.faces.at(i).at(3) - 1;

			for (unsigned int iCheck = 0; iCheck < this->selectedVertices.size(); iCheck++){ if (this->selectedVertices.at(iCheck) == buffer){ tripleCheck[1] = true; } }

			if (tripleCheck[1] == true){
				//FACE #2
				buffer = (int)input.faces.at(i).at(6) - 1;

				for (unsigned int iCheck = 0; iCheck < this->selectedVertices.size(); iCheck++){ if (this->selectedVertices.at(iCheck) == buffer){ tripleCheck[2] = true; } }

				if (tripleCheck[2] == true){
					toBeInsertedIntoFaces[0] = (int)this->assemblyToChildPartVerticesMatch.at((int)input.faces.at(i).at(0) - 1).at(1) + 1;
					toBeInsertedIntoFaces[1] = 0;
					toBeInsertedIntoFaces[2] = 0;// toBeInsertedIntoFaces[2] = (int)input.faces.at(i).at(2);
					toBeInsertedIntoFaces[3] = (int)this->assemblyToChildPartVerticesMatch.at((int)input.faces.at(i).at(3) - 1).at(1) + 1;
					toBeInsertedIntoFaces[4] = 0;
					toBeInsertedIntoFaces[5] = 0;// (int)input.faces.at(i).at(5);
					toBeInsertedIntoFaces[6] = (int)this->assemblyToChildPartVerticesMatch.at((int)input.faces.at(i).at(6) - 1).at(1) + 1;
					toBeInsertedIntoFaces[7] = 0;
					toBeInsertedIntoFaces[8] = 0;// (int)input.faces.at(i).at(8);

					this->faces.push_back({
						toBeInsertedIntoFaces[0],
						toBeInsertedIntoFaces[1],
						toBeInsertedIntoFaces[2],
						toBeInsertedIntoFaces[3],
						toBeInsertedIntoFaces[4],
						toBeInsertedIntoFaces[5],
						toBeInsertedIntoFaces[6],
						toBeInsertedIntoFaces[7],
						toBeInsertedIntoFaces[8]
					});
				}
			}
		}
	}
	return;
}

void  dotObj::getPositionsFromChild(dotObj &input){
	for (unsigned int i = 0; i < input.assemblyToChildPartVerticesMatch.size(); i++){
		if (((int)input.assemblyToChildPartVerticesMatch.at(i).at(0) > 0) && ((int)input.assemblyToChildPartVerticesMatch.at(i).at(1) > 0)){
			//std::cout << (int)input.assemblyToChildPartVerticesMatch.at(i).at(0) << "|"<< (int)input.assemblyToChildPartVerticesMatch.at(i).at(1) << std::endl; //Debug

			this->vertices.at((int)input.assemblyToChildPartVerticesMatch.at(i).at(0)).at(0) = input.vertices.at((int)input.assemblyToChildPartVerticesMatch.at(i).at(1)).at(0);
			this->vertices.at((int)input.assemblyToChildPartVerticesMatch.at(i).at(0)).at(1) = input.vertices.at((int)input.assemblyToChildPartVerticesMatch.at(i).at(1)).at(1);
			this->vertices.at((int)input.assemblyToChildPartVerticesMatch.at(i).at(0)).at(2) = input.vertices.at((int)input.assemblyToChildPartVerticesMatch.at(i).at(1)).at(2);
		}
	}

	return;
}


int dotObj::getPointIndex(float x, float y, float z){
	int theIndex = -2;
	for (unsigned int i = 0; i < this->vertices.size(); i++){
		if (this->functions.round2d(this->vertices.at(i).at(0)) == this->functions.round2d(x)){
			if (this->functions.round2d(this->vertices.at(i).at(1)) == this->functions.round2d(y)){
				if (this->functions.round2d(this->vertices.at(i).at(2)) == this->functions.round2d(z)){
					theIndex = i;
				}
			}
		}
	}

	if (theIndex == -2){
		std::cout << "No Point Found with index! Caution!" << std::endl;
	}

	return theIndex;
}


void dotObj::recalculateNormals(void){

	Vector3f v[3],e[2],n;

	this->normals.clear();
	this->normals.resize(0);
	for (int i = 0; i < this->faces.size(); i++){

		v[0] = Vector3f(this->vertices[this->faces[i][0]-1][0], this->vertices[this->faces[i][0]-1][1], this->vertices[this->faces[i][0]-1][2]);
		v[1] = Vector3f(this->vertices[this->faces[i][3]-1][0], this->vertices[this->faces[i][3]-1][1], this->vertices[this->faces[i][3]-1][2]);
		v[2] = Vector3f(this->vertices[this->faces[i][6]-1][0], this->vertices[this->faces[i][6]-1][1], this->vertices[this->faces[i][6]-1][2]);

		e[0] = v[2] - v[0];
		e[1] = v[1] - v[0];


		n = e[0].cross(e[1]);
		n.normalize();
		n = -n;
		if (isnan(n[0])){ n[0] = 0.0; }
		if (isnan(n[1])){ n[1] = 0.0; }
		if (isnan(n[2])){ n[2] = 0.0; }
		this->normals.push_back({ (float)n[0], (float)n[1], (float)n[2] });
		this->faces[i][2] = this->normals.size() ;
		this->faces[i][5] = this->normals.size() ;
		this->faces[i][8] = this->normals.size() ;

	}
	return;
}


double dotObj::getAverageFaceArea(void){
	double sum = 0;
	int a, b, c;
	for (unsigned int i = 0; i < this->faces.size(); i++)
	{
		a = this->faces.at(i).at(0) - 1;
		b = this->faces.at(i).at(3) - 1;
		c = this->faces.at(i).at(6) - 1;
		sum += computeTriangeArea(a, b, c);
	}

	return sum / this->faces.size(); //Division
}

double dotObj::computeTriangeArea(int firstIndex, int secondIndex, int thirdIndex){
	double area;
	double theta;
	Vector3d firstPoint, secondPoint, thirdPoint;
	for (int i = 0; i < 3; i++){
		firstPoint[i] = this->vertices.at(firstIndex).at(i);
		secondPoint[i] = this->vertices.at(secondIndex).at(i);
		thirdPoint[i] = this->vertices.at(thirdIndex).at(i);
	}
	theta = computeAngle((secondPoint - firstPoint), (thirdPoint - firstPoint));

	if (0.5*((double)((secondPoint - firstPoint).norm()))*((double)((thirdPoint - firstPoint).norm())) < 0.00001)
	{
		area = 0;
	}
	else
	{
		area = 0.5*((double)((secondPoint - firstPoint).norm()))*((double)((thirdPoint - firstPoint).norm()))*((double)sin(theta));
	}
	return area;
}

double dotObj::computeAngle(Vector3d a, Vector3d b){
	double angle;
	angle = acos((a.dot(b)) / (a.norm()*b.norm())); //Division

	return angle;
}

Matrix<double, Dynamic, Dynamic> dotObj::getRingAreas(void){
	//Get ring areas for each vertex of the obj
	MatrixXd I;
	I.setIdentity(this->vertices.size(), this->vertices.size());
	//#pragma omp parallel // shared(appendPoint,radiusMultiplier,v) private(p0,p1,proj)
	//#pragma omp for
	for (int i = 0; i < this->vertices.size(); i++){
		I(i, i) = getRingArea(i);
	}

	return I;
}

double dotObj::getRingArea(int vertex){
	//Get ring area of a certain vertex
	double area = 0;
	int a, b, c;
	for (int i = 0; i < this->faces.size(); i++){
		if ((this->faces.at(i).at(0) - 1 == vertex) || (this->faces.at(i).at(3) - 1 == vertex) || (this->faces.at(i).at(6) - 1 == vertex)){
			a = this->faces.at(i).at(0) - 1;
			b = this->faces.at(i).at(3) - 1;
			c = this->faces.at(i).at(6) - 1;
			area += computeTriangeArea(a, b, c);
		}
	}
	return area;
}

void dotObj::getDeltaCoordinates(Matrix<double, Dynamic, Dynamic> &M){
	Vector3f v;

	M.resize(this->vertices.size(), 3);
	M.setZero();

	for (int i = 0; i < this->vertices.size(); i++){
		if(i%10000==0)std::cout << i << std::endl;
		v= getDeltaCoordinate(i);
		M(i, 0) = v[0];
		M(i, 1) = v[1];
		M(i, 2) = v[2];

	}
	return;
}

void dotObj::getDeltaCoordinates(std::vector<Vector3f> &M){
	Vector3f v;
	M.resize(this->vertices.size());
	for (int i = 0; i < this->vertices.size(); i++){
		if (i % 10000 == 0)std::cout << i << std::endl;
		v = getDeltaCoordinate(i);
		M[i][0] = v[0];
		M[i][1] = v[1];
		M[i][2] = v[2];
	}
	return;
}





Vector3f dotObj::getDeltaCoordinate(int vertex){
		if (this->adjacencyMatrixSimplified.size() == 0){
		this->adjacencyMatrixSimplified = this->getAdjacencySimplified();
	}
	int v;
	std::vector<int> ns;
	Vector3f Ve,S;
	S.setZero();
	std::vector<float> Vo;
	for (int i = 1; i < this->adjacencyMatrixSimplified[vertex].size(); i++){
		ns=this->adjacencyMatrixSimplified[vertex];
		v = ns[i];
		Vo = this->vertices[v];
		Ve = functions.stdtoeigenvec(Vo);
		S = S + Ve;
	}
	
	S = S / (ns.size()-1);



	return S;
}

void dotObj::convToMat(Matrix<double, Dynamic, Dynamic> &V, std::vector<std::vector<float>>&vert, int direction){
	if (direction == 0){
		//#pragma omp parallel 
		//#pragma omp for
		for (int i = 0; i < this->vertices.size(); i++){
			for (int j = 0; j < 3; j++){
				V(i, j) = vert.at(i).at(j);
			}
		}
	}
	else{
		//#pragma omp parallel 
		//#pragma omp for
		for (int i = 0; i < this->vertices.size(); i++){
			for (int j = 0; j < 3; j++){
				vert.at(i).at(j) = V(i, j);
			}
		}
	}
	return;
}



std::vector<int> dotObj::findCommonNeighbours(int p1, int p2){
	//==== dotObj::findCommonNeighbours overview==============
	//This function returns a set of indices neighbouring the indices p1 & p2

	std::vector<int> N;
	int size = this->vertices.size();
	for (int i = 0; i < size; i++){
		if ((this->adjacencyMatrix(p1, i) == 1) && (this->adjacencyMatrix(p2, i) == 1) && (i != p1) && (i != p2)){
			N.push_back(i);
		}
	}
	return N;
}








Vector3f dotObj::calculateNormal(int faceIndexUnderInvestigation){
	Vector3f v[3], n;
	int i = faceIndexUnderInvestigation;
	v[0] = Vector3f(this->vertices[this->faces[i][0] - 1][0], this->vertices[this->faces[i][0] - 1][1], this->vertices[this->faces[i][0] - 1][2]);
	v[1] = Vector3f(this->vertices[this->faces[i][3] - 1][0], this->vertices[this->faces[i][3] - 1][1], this->vertices[this->faces[i][3] - 1][2]);
	v[2] = Vector3f(this->vertices[this->faces[i][6] - 1][0], this->vertices[this->faces[i][6] - 1][1], this->vertices[this->faces[i][6] - 1][2]);
	n = this->functions.getNormalFromTriangle(v[0], v[1], v[2]);
	return n;
}



//TODO
void dotObj::mergeCloseVertices(){
	//this->getBoundingBox();

	this->initCells(10);

	for (int i = 0; i < this->cells.cells.size(); i++){
	}

	return;
}
