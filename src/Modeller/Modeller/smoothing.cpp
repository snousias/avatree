#include "lungModelling.h"


void  dotObj::smoothingSparse(std::string type, std::string laplacian, int iterations, double lamda, double mu, std::vector<int> & selection){
	std::cout << "Smoothing" << "," << type << "," << laplacian << std::endl;
	int size = this->vertices.size();
	SparseMatrix<double> L;
	//SparseVector<double> Points;
	std::cout << "Get Adjacency" << std::endl;
	if (this->adjacencyMatrixSparse.size() == 0){
		this->adjacencyMatrixSparse = this->getAdjacencySparse();
	}
	VectorXd Points;
	Points.resize(size);

	

	for (int k = 0; k < iterations; k++){
		float previousperc = 0;
		float perc = 100 * ((float)(k * 100) / (float)iterations) / 100;
		int currentperc = (int)perc;
		if (currentperc % 10 == 0){
			if (currentperc != previousperc){
				currentperc = previousperc;
				std::cout << "Smoothing iteration :" << k << std::endl;
			}
		}

		if (laplacian == "onoff")
		{
			VectorXd D = this->adjacencyMatrixSparse * VectorXd::Ones(size);
			D = D.array().inverse().transpose();
			SparseMatrix<double> Dmat(size, size);// D.asDiagonal();
			std::vector<Triplet<double>> trp;
			for (int i = 0; i < D.size(); i++){
				trp.push_back(Triplet<double>(i, i, D(i)));
			}
			Dmat.setFromTriplets(trp.begin(), trp.end());
			SparseMatrix<double> K = Dmat*this->adjacencyMatrixSparse;
			SparseMatrix<double> I(size, size);
			I.setIdentity();
			L = I - K;

		}
		else if (laplacian == "curv")
		{
			//L = getLaplacianSparse("curv");  //Not functional
			this->adjacencyMatrixSparse = this->getAdjacencySparseWeighted();
			VectorXd D = this->adjacencyMatrixSparse * VectorXd::Ones(size);
			D = D.array().inverse().transpose();
			D = D.array();
			SparseMatrix<double> Dmat(size, size);// D.asDiagonal();
			std::vector<Triplet<double>> trp;
			for (int i = 0; i < D.size(); i++){
				trp.push_back(Triplet<double>(i, i, D(i)));
			}
			Dmat.setFromTriplets(trp.begin(), trp.end());
			SparseMatrix<double> K = Dmat*this->adjacencyMatrixSparse;
			SparseMatrix<double> I(size, size);
			I.setIdentity();
			L = I - K;
		}
		else if (laplacian == "weighted")
		{
			this->adjacencyMatrixSparse = this->getAdjacencySparse();
			VectorXd D = this->adjacencyMatrixSparse * VectorXd::Ones(size);
			D = D.array().inverse().transpose();
			SparseMatrix<double> Dmat(size, size);// D.asDiagonal();
			std::vector<Triplet<double>> trp;
			for (int i = 0; i < D.size(); i++){
				trp.push_back(Triplet<double>(i, i, D(i)));
			}
			Dmat.setFromTriplets(trp.begin(), trp.end());

			SparseMatrix<double> K = Dmat*this->adjacencyMatrixSparse;
			SparseMatrix<double> I(size, size);
			I.setIdentity();
			L = I - K;
		}
		for (int j = 0; j < 3; j++){
			for (int i = 0; i < size; i++){
					Points(i) = this->vertices.at(i).at(j);
			}
			if (type == "laplacian"){
				Points = Points + lamda*L*Points;
			}
			if (type == "taubin"){
				if (k % 2 == 0){
					Points = Points + lamda*L*Points;
				}
				else{
					Points = Points + mu*L*Points;
				}
			}

			for (int i = 0; i < size; i++){
				if (this->functions.valueExistsInVector(selection, i)){
					this->vertices.at(i).at(j) = Points(i);
				}
			}
		}
	}
	return;
}


void  dotObj::smoothingSparse(std::string type, std::string laplacian, int iterations, double lamda, double mu){
	std::cout << "Smoothing" << "," << type << "," << laplacian << std::endl;

	int size = this->vertices.size();
	SparseMatrix<double> L;
	//SparseVector<double> Points;
	VectorXd Points;
	Points.resize(size);

	std::cout << "Get Adjacency" << std::endl;
	/*if (this->adjacencyMatrixSparse.size() == 0){
		this->adjacencyMatrixSparse = this->getAdjacencySparse();
	}*/
	if (laplacian == "onoff")
	{
		this->adjacencyMatrixSparse = this->getAdjacencySparse();
		VectorXd D = this->adjacencyMatrixSparse * VectorXd::Ones(size);
		D = D.array().inverse().transpose();
		SparseMatrix<double> Dmat(size, size);// D.asDiagonal();
		std::vector<Triplet<double>> trp;
		for (int i = 0; i < D.size(); i++){
			trp.push_back(Triplet<double>(i, i, D(i)));
		}
		Dmat.setFromTriplets(trp.begin(), trp.end());
		SparseMatrix<double> K = Dmat*this->adjacencyMatrixSparse;
		SparseMatrix<double> I(size, size);
		I.setIdentity();
		L = I - K;

	}
	if (laplacian == "curv")
	{
		//L = getLaplacianSparse("curv");  //Not functional
		this->adjacencyMatrixSparse = this->getAdjacencySparseWeighted();
		VectorXd D = this->adjacencyMatrixSparse * VectorXd::Ones(size);
		D = D.array().inverse().transpose();
		D = D.array();
		SparseMatrix<double> Dmat(size, size);// D.asDiagonal();
		std::vector<Triplet<double>> trp;
		for (int i = 0; i < D.size(); i++){
			trp.push_back(Triplet<double>(i, i, D(i)));
		}
		Dmat.setFromTriplets(trp.begin(), trp.end());
		SparseMatrix<double> K = Dmat*this->adjacencyMatrixSparse;
		SparseMatrix<double> I(size, size);
		I.setIdentity();
		L = I - K;



	}
	if (laplacian == "weighted")
	{
		this->adjacencyMatrixSparse = this->getAdjacencySparse();
		VectorXd D = this->adjacencyMatrixSparse * VectorXd::Ones(size);
		D = D.array().inverse().transpose();
		SparseMatrix<double> Dmat(size, size);// D.asDiagonal();
		std::vector<Triplet<double>> trp;
		for (int i = 0; i < D.size(); i++){
			trp.push_back(Triplet<double>(i, i, D(i)));
		}
		Dmat.setFromTriplets(trp.begin(), trp.end());

		SparseMatrix<double> K = Dmat*this->adjacencyMatrixSparse;
		SparseMatrix<double> I(size, size);
		I.setIdentity();
		L = I - K;
	}




	for (int k = 0; k < iterations; k++){
		float previousperc = 0;
		float perc = 100 * ((float)(k * 100) / (float)iterations) / 100;
		int currentperc = (int)perc;
		if (currentperc % 10 == 0){
			if (currentperc != previousperc){
				currentperc = previousperc;
				std::cout << "Smoothing iteration :" << k << std::endl;
			}
		}

		

		for (int j = 0; j < 3; j++){
			for (int i = 0; i < size; i++){
				//Points.coeffRef(i)= this->vertices.at(i).at(j);
				Points(i) = this->vertices.at(i).at(j);
			}
			if (type == "laplacian"){
				Points = Points + lamda*L*Points;
			}
			if (type == "taubin"){
				if (k % 2 == 0){
					Points = Points + lamda*L*Points;
				}
				else{
					Points = Points + mu*L*Points;
				}
			}

			for (int i = 0; i < size; i++){
				this->vertices.at(i).at(j) = Points(i);
			}
		}
	}

	return;
}

void  dotObj::smoothing(std::string type, std::string laplacian, int iterations, double lamda, double mu){
	std::cout << "Smoothing" << "," << type << "," << laplacian << std::endl;
	int size = this->vertices.size();
	MatrixXd L;
	VectorXd Points;
	Points.resize(size);


	if (this->adjacencyMatrix.size() == 0){
		std::cout << "Get Adjacency" << std::endl;
		this->adjacencyMatrix = this->getAdjacency();
	}
	std::cout << "Get Laplacian" << std::endl;

	if (laplacian == "onoff")
	{
		L = getLaplacian("onoff");
	}
	else if (laplacian == "curv")
	{
		L = getLaplacian("curv");
	}

	for (int k = 0; k < iterations; k++){
		std::cout << 100 * ((float)(k * 100) / (float)iterations) / 100 << "%" << std::endl;

		for (int j = 0; j < 3; j++){
			for (int i = 0; i < size; i++){
				Points(i) = this->vertices.at(i).at(j);
			}
			if (type == "laplacian"){
				Points = Points + lamda*L*Points;
			}
			if (type == "taubin"){
				Points = Points + lamda*L*Points;
				Points = Points + mu*L*Points;
			}

			for (int i = 0; i < size; i++){
				this->vertices.at(i).at(j) = Points(i);
			}
		}
	}

	return;
}

void dotObj::taubin(int iterations, double lamda, double mu){

	int size = this->vertices.size();
	Matrix<double, Dynamic, Dynamic> L;
	VectorXd Points;
	Points.setZero(size);

	if (this->adjacencyMatrix.size() == 0){
		this->adjacencyMatrix = this->getAdjacency();
	}

	L = getLaplacian("onoff");
	//L = getLaplacian("curv");

	for (int k = 0; k < iterations; k++){

		for (int j = 0; j < 3; j++){

			for (int i = 0; i < size; i++){
				Points[i] = this->vertices.at(i).at(j);
			}


			Points = Points - lamda*L*Points;
			Points = Points + mu*L*Points;


			for (int i = 0; i < size; i++){
				this->vertices.at(i).at(j) = Points[i];
			}
		}
	}
	return;
}



