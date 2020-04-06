#include "lungModelling.h"


std::vector<int> dotObj::getFaceIndexContainingVertexIndex(int vertexIndex){

	std::vector<int> result;
	
	for (int k = 0; k < this->lines.size(); k++){



	}
		
	return result;
}




std::vector<int> dotObj::getLineIndexContainingVertexIndex(int vertexIndex){
	
	std::vector<int> result;

	for (int k = 0; k < this->lines.size(); k++){
		if ((this->lines[k][0] == vertexIndex + 1) || (this->lines[k][1] == vertexIndex + 1)){
			result.push_back(k);
		}
		
	}
	
	return result;
}


