#include "lungModelling.h"



//Create the adjacency matrix for the current dotObj
//Store the adjacency in a dense matrix
//Export matrix
//No data is stored in the dotObj
MatrixXd dotObj::getAdjacency(void){
	Matrix<double, Dynamic, Dynamic> Adj;
	int size = this->vertices.size();
	Adj.setZero(size, size);

	//For each face per 2 vetrices define a cell in the adjacency matrix
	for (unsigned int i = 0; i < this->faces.size(); i++){
		//------------------------------------------------------------------------
		Adj(this->faces.at(i).at(0) - 1, this->faces.at(i).at(3) - 1) = 1;
		Adj(this->faces.at(i).at(3) - 1, this->faces.at(i).at(6) - 1) = 1;
		Adj(this->faces.at(i).at(6) - 1, this->faces.at(i).at(0) - 1) = 1;
		//------------------------------------------------------------------------
		Adj(this->faces.at(i).at(0) - 1, this->faces.at(i).at(6) - 1) = 1;
		Adj(this->faces.at(i).at(6) - 1, this->faces.at(i).at(3) - 1) = 1;
		Adj(this->faces.at(i).at(3) - 1, this->faces.at(i).at(0) - 1) = 1;
	}

	return Adj;
}

//Function to pass all adjacency matrices to obj object
//Boolean switch controls which matrix to store
void dotObj::getAdjacencyMatrices(bool dense, bool vector, bool sparse){
	if (dense){ this->adjacencyMatrix = this->getAdjacency(); }

	if (vector){ this->adjacencyMatrixSimplified = this->getAdjacencySimplified(); }

	if (sparse){ this->adjacencyMatrixSparse = this->getAdjacencySparse(); }

	return;
}

//Adjacency Sparse matrix
//Create the adjacency matrix for the current dotObj
//Store the adjacency in a sparse matrix
//Export matrix
//No data is stored in the dotObj
SparseMatrix<double> dotObj::getAdjacencySparse(void){
	int size = this->vertices.size();
	SparseMatrix<double> Adj(size, size);
	Adj.setZero();
	std::vector<Triplet<double>> trip;
	for (int i = 0; i < faces.size(); i++){
		trip.push_back(Triplet<double>((double)this->faces.at(i).at(0) - 1, (double)this->faces.at(i).at(3) - 1, (double)1));
		trip.push_back(Triplet<double>((double)this->faces.at(i).at(3) - 1, (double)this->faces.at(i).at(6) - 1, (double)1));
		trip.push_back(Triplet<double>((double)this->faces.at(i).at(6) - 1, (double)this->faces.at(i).at(0) - 1, (double)1));
		//------------------------------------------------------------------------
		trip.push_back(Triplet<double>((double)this->faces.at(i).at(0) - 1, (double)this->faces.at(i).at(6) - 1, (double)1));
		trip.push_back(Triplet<double>((double)this->faces.at(i).at(6) - 1, (double)this->faces.at(i).at(3) - 1, (double)1));
		trip.push_back(Triplet<double>((double)this->faces.at(i).at(3) - 1, (double)this->faces.at(i).at(0) - 1, (double)1));
	}

	Adj.setFromTriplets(trip.begin(), trip.end());
	return Adj;
}


SparseMatrix<double> dotObj::getAdjacencySparse(std::vector<int> &selectionList){
	int size = this->vertices.size();
	SparseMatrix<double> Adj(size, size);
	Adj.setZero();
	std::vector<Triplet<double>> trip;
	for (int i = 0; i < faces.size(); i++){
		int a = this->faces.at(i).at(0) - 1;
		int b = this->faces.at(i).at(3) - 1;
		int c = this->faces.at(i).at(6) - 1;
		if (this->functions.valueExistsInVector(selectionList, a) &&
			this->functions.valueExistsInVector(selectionList, b) &&
			this->functions.valueExistsInVector(selectionList, c)){
			trip.push_back(Triplet<double>((double)this->faces.at(i).at(0) - 1, (double)this->faces.at(i).at(3) - 1, (double)1));
			trip.push_back(Triplet<double>((double)this->faces.at(i).at(3) - 1, (double)this->faces.at(i).at(6) - 1, (double)1));
			trip.push_back(Triplet<double>((double)this->faces.at(i).at(6) - 1, (double)this->faces.at(i).at(0) - 1, (double)1));
			//------------------------------------------------------------------------
			trip.push_back(Triplet<double>((double)this->faces.at(i).at(0) - 1, (double)this->faces.at(i).at(6) - 1, (double)1));
			trip.push_back(Triplet<double>((double)this->faces.at(i).at(6) - 1, (double)this->faces.at(i).at(3) - 1, (double)1));
			trip.push_back(Triplet<double>((double)this->faces.at(i).at(3) - 1, (double)this->faces.at(i).at(0) - 1, (double)1));
		}
	}

	Adj.setFromTriplets(trip.begin(), trip.end());
	return Adj;
}


//Weighted Adjacency Sparse matrix
//Create the weighted adjacency matrix for the current dotObj
//Store the adjacency in a sparse matrix
//Export matrix
//No data is stored in the dotObj
SparseMatrix<double> dotObj::getAdjacencySparseWeighted(void){
	int size = this->vertices.size();
	SparseMatrix<double> Adj(size, size);
	Adj.setZero();
	Vector3d p1, p2;
	double distancePoints;
	std::vector<Triplet<double>> trip;
	for (int i = 0; i < faces.size(); i++){
		for (int j = 0; j < 3; j++){
			p1[j] = this->vertices[this->faces.at(i).at(0) - 1][j];
			p2[j] = this->vertices[this->faces.at(i).at(3) - 1][j];
		}
		distancePoints = (p1 - p2).norm();
		trip.push_back(Triplet<double>((double)this->faces.at(i).at(0) - 1, (double)this->faces.at(i).at(3) - 1, distancePoints));
		trip.push_back(Triplet<double>((double)this->faces.at(i).at(3) - 1, (double)this->faces.at(i).at(0) - 1, distancePoints));
		//===================================
		for (int j = 0; j < 3; j++){
			p1[j] = this->vertices[this->faces.at(i).at(3) - 1][j];
			p2[j] = this->vertices[this->faces.at(i).at(6) - 1][j];
		}
		distancePoints = (p1 - p2).norm();
		trip.push_back(Triplet<double>((double)this->faces.at(i).at(3) - 1, (double)this->faces.at(i).at(6) - 1, distancePoints));
		trip.push_back(Triplet<double>((double)this->faces.at(i).at(6) - 1, (double)this->faces.at(i).at(3) - 1, distancePoints));
		//===================================
		for (int j = 0; j < 3; j++){
			p1[j] = this->vertices[this->faces.at(i).at(0) - 1][j];
			p2[j] = this->vertices[this->faces.at(i).at(6) - 1][j];
		}
		distancePoints = (p1 - p2).norm();
		trip.push_back(Triplet<double>((double)this->faces.at(i).at(0) - 1, (double)this->faces.at(i).at(6) - 1, distancePoints));
		trip.push_back(Triplet<double>((double)this->faces.at(i).at(6) - 1, (double)this->faces.at(i).at(0) - 1, distancePoints));
		//===================================
	}

	Adj.setFromTriplets(trip.begin(), trip.end());
	return Adj;
}

//Create the triplet list to be used as input for the adjacency matrix
std::vector<std::vector<int>> dotObj::getAdjacencyTripletList(void){
	std::vector<std::vector<int>> trip;
	for (unsigned int i = 0; i < this->faces.size(); i++){
		//------------------------------------------------------------------------
		trip.push_back({ this->faces.at(i).at(0) - 1, this->faces.at(i).at(3) - 1, 1 });
		trip.push_back({ this->faces.at(i).at(3) - 1, this->faces.at(i).at(6) - 1, 1 });
		trip.push_back({ this->faces.at(i).at(6) - 1, this->faces.at(i).at(0) - 1, 1 });
		//------------------------------------------------------------------------
		trip.push_back({ this->faces.at(i).at(0) - 1, this->faces.at(i).at(6) - 1, 1 });
		trip.push_back({ this->faces.at(i).at(6) - 1, this->faces.at(i).at(3) - 1, 1 });
		trip.push_back({ this->faces.at(i).at(3) - 1, this->faces.at(i).at(0) - 1, 1 });
	}

	//trip.erase(unique(trip.begin(), trip.end()), trip.end());

	return trip;
}

//Create simplified adjacency matrix
//Deploy a vector of vectors
// 0 4 3 2 1  Vector 0 is adjacent to 4,3,2,1
// 1 3 5 6 4
// 2 3 7 8 9
// 3 2 6 3 1
// 4 3 0 9 7
std::vector<std::vector<int>> dotObj::getAdjacencySimplified(void){
	int temp;
	int searchable;
	std::vector<std::vector<int>> adjacencyMatrixSimplified;
	std::vector<int>::iterator it;

	for (int i = 0; i < (int)vertices.size(); i++){
		adjacencyMatrixSimplified.push_back({ i + 1 });
	}

	for (int i = 0; i < (int)faces.size(); i++){
		//=====Select vertice in faces vector(first point)=====================================================
		temp = (int)faces.at(i).at(0) - 1;
		//Second point of triange
		searchable = (int)faces.at(i).at(3);
		it = find(adjacencyMatrixSimplified.at(temp).begin(), adjacencyMatrixSimplified.at(temp).end(), searchable); //Check if exists
		if (it == adjacencyMatrixSimplified.at(temp).end()){
			adjacencyMatrixSimplified.at(temp).push_back(searchable); //If not exist push
		}
		//Third point of triange
		searchable = (int)faces.at(i).at(6);
		it = find(adjacencyMatrixSimplified.at(temp).begin(), adjacencyMatrixSimplified.at(temp).end(), searchable); //Check if exists
		if (it == adjacencyMatrixSimplified.at(temp).end()){
			adjacencyMatrixSimplified.at(temp).push_back(searchable);//If not exist push
		}

		//=====Select vertice in faces vector(second point)=====================================================
		temp = (int)faces.at(i).at(3) - 1;

		//Second point of triange
		searchable = (int)faces.at(i).at(0);
		it = find(adjacencyMatrixSimplified.at(temp).begin(), adjacencyMatrixSimplified.at(temp).end(), searchable); //Check if exists
		if (it == adjacencyMatrixSimplified.at(temp).end()){
			adjacencyMatrixSimplified.at(temp).push_back(searchable);//If not exist push
		}

		//Third point of triange
		searchable = (int)faces.at(i).at(6);
		it = find(adjacencyMatrixSimplified.at(temp).begin(), adjacencyMatrixSimplified.at(temp).end(), searchable); //Check if exists
		if (it == adjacencyMatrixSimplified.at(temp).end()){
			adjacencyMatrixSimplified.at(temp).push_back(searchable);//If not exist push
		}

		//=====Select vertice in faces vector(third point)=====================================================
		temp = (int)faces.at(i).at(6) - 1;

		//Second point of triange
		searchable = (int)faces.at(i).at(3);
		it = find(adjacencyMatrixSimplified.at(temp).begin(), adjacencyMatrixSimplified.at(temp).end(), searchable);//Check if exists
		if (it == adjacencyMatrixSimplified.at(temp).end()){
			adjacencyMatrixSimplified.at(temp).push_back(searchable);//If not push
		}

		//Third point of triange
		searchable = (int)faces.at(i).at(0);
		it = find(adjacencyMatrixSimplified.at(temp).begin(), adjacencyMatrixSimplified.at(temp).end(), searchable);//Check if exists
		if (it == adjacencyMatrixSimplified.at(temp).end()){
			adjacencyMatrixSimplified.at(temp).push_back(searchable);//If not push
		}
	}

	for (unsigned int i = 0; i < adjacencyMatrixSimplified.size(); i++){
		for (unsigned int j = 0; j < adjacencyMatrixSimplified.at(i).size(); j++){
			adjacencyMatrixSimplified.at(i).at(j) = adjacencyMatrixSimplified.at(i).at(j) - 1;
		}
	}

	// 0 4 3 2 1  Vector 0 is adjacent to 4,3,2,1
	// 1 3 5 6 4
	// 2 3 7 8 9
	// 3 2 6 3 1
	// 4 3 0 9 7

	return adjacencyMatrixSimplified;
}




std::vector<int> dotObj::findAdjacentFaces(int faceIndexUnderInvestigation){
	std::vector<int> adjacentFaces;
	bool pass = false;
	int a = this->faces.at(faceIndexUnderInvestigation).at(0);
	int b = this->faces.at(faceIndexUnderInvestigation).at(3);
	int c = this->faces.at(faceIndexUnderInvestigation).at(6);

	for (int j = 0; j < this->faces.size(); j++){
		if (j != faceIndexUnderInvestigation){
			int d = this->faces.at(j).at(0);
			int e = this->faces.at(j).at(3);
			int f = this->faces.at(j).at(6);

			if ((a == d) || (a == e) || (a == f)){
				pass = true;
				for (int i = 0; i < adjacentFaces.size(); i++)
				{
					if ((pass) && (j == adjacentFaces.at(i)))
					{
						pass = false;
					}
				}
				if (pass)
				{
					adjacentFaces.push_back(j);
				}
				//adjacentFaces.push_back(j);
			}

			if ((b == d) || (b == e) || (b == f)){
				pass = true;
				for (int i = 0; i < adjacentFaces.size(); i++)
				{
					if ((pass) && (j == adjacentFaces.at(i)))
					{
						pass = false;
					}
				}
				if (pass)
				{
					adjacentFaces.push_back(j);
				}
				//adjacentFaces.push_back(j);
			}

			if ((c == d) || (c == e) || (c == f)){
				pass = true;
				for (int i = 0; i < adjacentFaces.size(); i++)
				{
					if ((pass) && (j == adjacentFaces.at(i)))
					{
						pass = false;
					}
				}
				if (pass)
				{
					adjacentFaces.push_back(j);
				}
				//adjacentFaces.push_back(j);
			}
		}
	}
	return adjacentFaces;
}

std::vector<int> dotObj::findAdjacentFacesPerVertex(int vertexIndexUnderInvestigation){
	std::vector<int> adjacentFaces;
	bool pass = false;
	int a = vertexIndexUnderInvestigation;

	for (int j = 0; j < this->faces.size(); j++){
		int d = this->faces.at(j).at(0);
		int e = this->faces.at(j).at(3);
		int f = this->faces.at(j).at(6);

		if ((a == d - 1) || (a == e - 1) || (a == f - 1)){
			pass = true;
			for (int i = 0; i < adjacentFaces.size(); i++)
			{
				if ((pass) && (j == adjacentFaces.at(i)))
				{
					pass = false;
				}
			}
			if (pass)
			{
				adjacentFaces.push_back(j);
			}
		}
	}
	return adjacentFaces;
}



std::vector<int> dotObj::searchImmediateNeighboors(int startingpoint){

	std::vector<int> neighboors;

	if (this->adjacencyMatrixSimplified.size() == 0){

		this->adjacencyMatrixSimplified = this->getAdjacencySimplified();
	}

	for (int i = 1; i < this->adjacencyMatrixSimplified[startingpoint].size(); i++){

		neighboors.push_back(this->adjacencyMatrixSimplified[startingpoint][i]);

	}


	return neighboors;
}

std::vector<int> dotObj::searchImmediateLineNeighboors(int startingpoint){

	std::vector<int> neighboors;

	if (this->lines.size() > 0){

		for (int i = 0; i < this->lines.size(); i++){

			if (this->lines.at(i).at(0) == startingpoint){ neighboors.push_back(this->lines.at(i).at(1)); }
			if (this->lines.at(i).at(1) == startingpoint){ neighboors.push_back(this->lines.at(i).at(0)); }
		}
	}
	return  neighboors;
}





