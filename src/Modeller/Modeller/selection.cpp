#include "lungModelling.h"



void dotObj::selector(int startingPointIndex, int numberOfPoints){
	int startingpoint = startingPointIndex;
	std::vector<int>::iterator it;

	it = find(selectedVertices.begin(), selectedVertices.end(), startingpoint);
	if (it == selectedVertices.end()){
		selectedVertices.push_back(startingpoint);
		this->searchpoints(selectedVertices.size() - 1, numberOfPoints);
	}
	else{
		this->searchpoints(it - selectedVertices.begin(), numberOfPoints);
	}

	for (unsigned int i = 0; i < selectedVertices.size(); i++){
		if ((int)selectedVertices.size() < numberOfPoints){
			this->searchpoints(i, numberOfPoints);
		}
	}
	std::cout << selectedVertices.size() << " points Selected" << std::endl;


	/*
	std::vector<int> kkkk;
	for (int m = 0; m < this->faces.size(); m++){
		for (int k = 0; k < selectedVertices.size(); k++){
			if (this->faces[m][0]-1 == selectedVertices[k]){
				kkkk.push_back(m);
			}
			else if (this->faces[m][3]-1 == selectedVertices[k]){
				kkkk.push_back(m);
			}
			else if (this->faces[m][6]-1 == selectedVertices[k]){
				kkkk.push_back(m);
			}
		}
	}
	selectedVertices = kkkk;
	*/




	return;
}



void dotObj::searchpoints(int j, int numberOfPoints){
	//Search for adjacent points
	//Select point to be investigated by the index of the point on the faces array
	std::vector<int>::iterator it;
	for (unsigned int i = 0; i < this->faces.size(); i++){
		if ((this->faces.at(i).at(0) - 1 == selectedVertices.at(j)) || (this->faces.at(i).at(3) - 1 == selectedVertices.at(j)) || (this->faces.at(i).at(6) - 1 == selectedVertices.at(j))){


			//selectedVertices.push_back(this->faces.at(i).at(0));
			//selectedVertices.push_back(this->faces.at(i).at(3));
			//selectedVertices.push_back(this->faces.at(i).at(6));

			it = find(selectedVertices.begin(), selectedVertices.end(), this->faces.at(i).at(0) - 1);
			if (it == selectedVertices.end()){
				selectedVertices.push_back(this->faces.at(i).at(0) - 1);
				if ((int)selectedVertices.size() == numberOfPoints){
					return;
				}
			}

			it = find(selectedVertices.begin(), selectedVertices.end(), this->faces.at(i).at(3) - 1);
			if (it == selectedVertices.end()){
				selectedVertices.push_back(this->faces.at(i).at(3) - 1);
				if ((int)selectedVertices.size() == numberOfPoints){
					return;
				}
			}

			it = find(selectedVertices.begin(), selectedVertices.end(), this->faces.at(i).at(6) - 1);
			if (it == selectedVertices.end()){
				selectedVertices.push_back(this->faces.at(i).at(6) - 1);
				if ((int)selectedVertices.size() == numberOfPoints){
					return;
				}
			}
		}
		else{
			//Do Nothing ... for now!
		}
	}
	return;
}




void dotObj::searchpoints(int j){
	//Search for adjacent points
	//Select point to be investigated by the index of the point on the faces array
	std::vector<int>::iterator it;
	for (unsigned int i = 0; i < this->faces.size(); i++){
		if ((this->faces.at(i).at(0) - 1 == selectedVertices.at(j)) || (this->faces.at(i).at(3) - 1 == selectedVertices.at(j)) || (this->faces.at(i).at(6) - 1 == selectedVertices.at(j))){


			//selectedVertices.push_back(this->faces.at(i).at(0));
			//selectedVertices.push_back(this->faces.at(i).at(3));
			//selectedVertices.push_back(this->faces.at(i).at(6));

			it = find(selectedVertices.begin(), selectedVertices.end(), this->faces.at(i).at(0) - 1);
			if (it == selectedVertices.end()){
				selectedVertices.push_back(this->faces.at(i).at(0) - 1);
			}

			it = find(selectedVertices.begin(), selectedVertices.end(), this->faces.at(i).at(3) - 1);
			if (it == selectedVertices.end()){
				selectedVertices.push_back(this->faces.at(i).at(3) - 1);
			}

			it = find(selectedVertices.begin(), selectedVertices.end(), this->faces.at(i).at(6) - 1);
			if (it == selectedVertices.end()){
				selectedVertices.push_back(this->faces.at(i).at(6) - 1);
			}
		}
		else{
			//Do Nothing ... for now!
		}
	}
	return;
}



void dotObj::extendSelection(){
	std::vector<int> randomSelection;
	std::vector<int> newSelection;

	for (int i = 0; i < this->selectedVertices.size(); i++)
	{
		randomSelection = searchImmediateNeighboors(this->selectedVertices.at(i));
		//std::cout << this->selectedVertices.at(i) << std::endl;
		newSelection.insert(newSelection.end(), randomSelection.begin(), randomSelection.end());
	}

	this->selectedVertices.insert(this->selectedVertices.end(), newSelection.begin(), newSelection.end());
	sort(this->selectedVertices.begin(), this->selectedVertices.end());
	this->selectedVertices.erase(unique(this->selectedVertices.begin(), this->selectedVertices.end()), this->selectedVertices.end());


	std::unique(this->selectedVertices.begin(), this->selectedVertices.end());

	return;
}




//TODO
void dotObj::selectorBasedOnCenterline(int startingPointIndex, double distance){
	std::vector<float> raw = this->getSDFPropertyMap(0.000002, 10, false);  //(4, 25, true, 4, 0.0001);
	int v = 0;
	bool found = false;
	for (int i = 0; i < this->faces.size(); i++){
		if (!found){
			if ((this->faces[i][0] == startingPointIndex + 1) || (this->faces[i][3] == startingPointIndex + 1) || (this->faces[i][6] == startingPointIndex + 1))
			{
				v = i;

				found = true;
			}
		}
	}
	std::cout << raw[v] << std::endl;
	raw[v];
	Vector3f n = this->calculateNormal(v);
	n.normalize();
	int a = this->faces[v][0] - 1;
	int b = this->faces[v][3] - 1;
	int c = this->faces[v][6] - 1;
	Vector3f v1, v2, v3, vf;
	v1[0] = this->vertices[a][0];
	v1[1] = this->vertices[a][1];
	v1[2] = this->vertices[a][2];
	v2[0] = this->vertices[b][0];
	v2[1] = this->vertices[b][1];
	v2[2] = this->vertices[b][2];
	v3[0] = this->vertices[c][0];
	v3[1] = this->vertices[c][1];
	v3[2] = this->vertices[c][2];
	Vector3f centroid = (v1 + v2 + v3) / 3;

	float radius = (raw[v] / 2.0);

	Vector3f rot = -radius* n;

	Vector3f center = centroid - rot;

	std::cout << v1[0] << "-" << v1[1] << "-" << v1[2] << "-" << std::endl;
	std::cout << centroid[0] << "-" << centroid[1] << "-" << centroid[2] << "-" << std::endl;
	std::cout << center[0] << "-" << center[1] << "-" << center[2] << "-" << std::endl;

	//Vector3f rot2 = -2 * rot;
	for (int i = 0; i < this->vertices.size(); i++){
		vf[0] = this->vertices[i][0];
		vf[1] = this->vertices[i][1];
		vf[2] = this->vertices[i][2];

		//if ((vf - center).norm() < 2.5 * radius){
		if ((vf - center).norm() < distance* radius){
			selectedVertices.push_back(i);
		}
	}

	return;
}

void dotObj::selectorBasedOnGeodesic(int startingPointIndex, double distance){
	std::vector<double> W;
	W = this->getMinDists(startingPointIndex);

	for (int i = 0; i < W.size(); i++){
		if (W[i] < distance){
			selectedVertices.push_back(i);
		}
	}

	return;
}



void dotObj::getNeighboursByFaceTri(int faceIndex, int numOfNeighbours) {

	this->findAdjacentFaces(faceIndex);




	return;
}