#include "lungModelling.h"

//======Connectivity surgery=======================================

MatrixXf dotObj::getK(int i, int j){
	MatrixXf K(3, 4);
	Vector3f vi(this->vertices.at(i).at(0), this->vertices.at(i).at(1), this->vertices.at(i).at(2));
	Vector3f vj(this->vertices.at(j).at(0), this->vertices.at(j).at(1), this->vertices.at(j).at(2));
	Vector3f a;
	Vector3f b;
	a = vj - vi;
	//vi.normalize();
	a.normalize();
	b = vi.cross(a);
	//b.normalize();
	K << 0, -a(2), a(1), -b(0), a(2), 0, -a(0), -b(1), -a(1), a(0), 0, -b(2);
	return K;
}

MatrixXf dotObj::getQ(int vertexIndex){
	MatrixXf Q(4, 4), K;
	Q.setZero();
	if (this->adjacencyMatrixSimplified.size() == 0){ this->adjacencyMatrixSimplified = this->getAdjacencySimplified(); }

	//for (int i = 0; i < this->adjacencyMatrixSimplified.size(); i++){
	for (int j = 1; j < this->adjacencyMatrixSimplified.at(vertexIndex).size(); j++){ //Neighboor indices start from the second element of the line
		int v2 = this->adjacencyMatrixSimplified.at(vertexIndex).at(j);
		K = getK(vertexIndex, v2);
		Q = Q + K.transpose()*K;
	}
	//}
	return Q;
}


/*
void dotObj::getQs(void){
MatrixXf Q(4, 4);
for (int j = 0; j < this->vertices.size(); j++){
Q = this->getQ(j);
Qs.push_back(Q);
//std::cout << Q << std::endl;
}

return;
}*/

//======Shape cost for edge i,j=============
float dotObj::getErrorMetricOfVertexToPoint(int i, Vector3f p){
	Vector4f point;
	float errorMetrix;
	point[0] = p[0];
	point[1] = p[1];
	point[2] = p[2];
	point[3] = 1.0;
	MatrixXf Q(4, 4);
	Q = getQ(i);
	errorMetrix = point.transpose()*Q*point;

	return errorMetrix;
	//return p.transpose()*this->Qs.at(i)p;
}

float dotObj::getShapeCost(int i, int j){
	float shapeCost;

	Vector3f vj(this->vertices.at(j).at(0), this->vertices.at(j).at(1), this->vertices.at(j).at(2));
	//vj.normalize();
	shapeCost = getErrorMetricOfVertexToPoint(i, vj) + getErrorMetricOfVertexToPoint(j, vj);

	return shapeCost;
}




//======Sampling cost for edge i,j=============

float dotObj::getSamplingCost(int i, int j){
	float samplingCostSUM,samplingCost = 0.0;

	Vector3f vi, vj, vk;
	vi[0] = this->vertices.at(i).at(0);
	vi[1] = this->vertices.at(i).at(1);
	vi[0] = this->vertices.at(i).at(2);
	//vi.normalize();
	vj[0] = this->vertices.at(j).at(0);
	vj[1] = this->vertices.at(j).at(1);
	vj[2] = this->vertices.at(j).at(2);
	//vj.normalize();
	float multip=(vi - vj).norm();
	int v2 = 0;
	int vertexIndex = i;

	for (int j = 1; j < this->adjacencyMatrixSimplified.at(vertexIndex).size(); j++){ //Neighboor indices start from the second element of the line
		v2 = this->adjacencyMatrixSimplified.at(vertexIndex).at(j);
		vk[0] = this->vertices.at(v2).at(0);
		vk[1] = this->vertices.at(v2).at(1);
		vk[2] = this->vertices.at(v2).at(2);
		//vk.normalize();
		samplingCost += (vi - vk).norm();
	}
	samplingCostSUM = multip *samplingCost;
	return samplingCostSUM;
}









float dotObj::computeTotalCostForEdge(int i, int j){
	float totalCost = 0;

	float shapeCost = this->getShapeCost(i, j);
	float samplingCost = this->getSamplingCost(i, j);
	totalCost = shapeCost + 0.1*samplingCost;

	return totalCost;
}







void dotObj::halfEdgeCollapse(int from, int to){
	std::vector<int> toBeDeleted;

	

	for (int i = 0; i < this->faces.size(); i++){
		if ((this->faces.at(i).at(0) == from+1) || (this->faces.at(i).at(3) == from+1) || (this->faces.at(i).at(6) == from+1)){
			if ((this->faces.at(i).at(0) == to+1) || (this->faces.at(i).at(3) == to+1) || (this->faces.at(i).at(6) == to+1)){
					    toBeDeleted.push_back(i);
						//facesToBeDeleted.push_back(i);
					this->vertices.at(from) = this->vertices.at(to);
			}
			else
			{
				//if (this->faces.at(i).at(0) == from){ this->faces.at(i).at(0) = to; }
				//if (this->faces.at(i).at(3) == from){ this->faces.at(i).at(3) = to; }
				//if (this->faces.at(i).at(6) == from){ this->faces.at(i).at(6) = to; }
			}
		}
	}


	
	for (int j = 0; j < toBeDeleted.size(); j++){
		if (this->faces.size()>0){
			//this->faces.erase(this->faces.begin() + toBeDeleted.at(j));
		}
	}

	return;
}

void dotObj::searchAndRemoveMinimum(void){
	if (this->edges.size() == 0){ this->getEdges(); }
	float minCost = 0.0;
	int minCostVertexIndex = 0;
	float temp;
	for (int i = 0; i < this->edges.size(); i++){
		bool start = true;
		for (int j = 0; j < this->edgesToBeDeleted.size(); j++){ if (i == edgesToBeDeleted.at(j)){ start = false; } }
		if (start){
			if (i == 0){
 				minCost = computeTotalCostForEdge(this->edges.at(i).at(0), this->edges.at(i).at(1));
				minCostVertexIndex = i;
			}
			else{
				temp = computeTotalCostForEdge(this->edges.at(i).at(0), this->edges.at(i).at(1));
				if (temp < minCost){ minCost = temp; minCostVertexIndex = i; }
			}
		}
	}

	edgesToBeDeleted.push_back(minCostVertexIndex);
	halfEdgeCollapse(this->edges.at(minCostVertexIndex).at(0), this->edges.at(minCostVertexIndex).at(1));
	
	return;
}

//----------------------------------------------------------------

/*
void dotObj::cleanUp(void){
std::vector<int> importantVectors;

for (int i = 0; i<this->faces.size(); i++){
importantVectors.push_back(this->faces.at(i).at(0));
importantVectors.push_back(this->faces.at(i).at(3));
importantVectors.push_back(this->faces.at(i).at(6));
}
sort(importantVectors.begin(), importantVectors.end());
importantVectors.erase(unique(importantVectors.begin(), importantVectors.end()), importantVectors.end());

std::vector<std::vector<float>> newVertices;

for (int j = 0; j<importantVectors.at(j); j++){
newVertices.push_back(this->vertices.at(importantVectors.at(j) - 1));
}

for (int i = 0; i<this->faces.size(); i++){
for (int j = 0; j<importantVectors.at(j); j++){
if (this->faces.at(i).at(0) == importantVectors.at(j)){ this->faces.at(i).at(0) = j + 1; }
if (this->faces.at(i).at(3) == importantVectors.at(j)){ this->faces.at(i).at(0) = j + 1; }
if (this->faces.at(i).at(6) == importantVectors.at(j)){ this->faces.at(i).at(0) = j + 1; }
}
}
}

bool dotObj::faceHasVertex(int faceIndex,int vertexIndex){
bool thereIs=false;

if ((this->faces.at(faceIndex).at(0) == vertexIndex) || (this->faces.at(faceIndex).at(3) == vertexIndex) || (this->faces.at(faceIndex).at(6) == vertexIndex)){
thereIs = true;
}

return thereIs;
}

*/