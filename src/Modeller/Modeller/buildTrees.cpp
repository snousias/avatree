#include "lungModelling.h"



void dotObj::buildInitialTrees(dotObj &mod1, dotObj &mod2, int gens){
	int init = 1;
	float maxz = 0;
	float maxzIndex = 0;
	std::vector<int> neighBours;
	std::vector<int> temp;
	std::vector<int> examined;
	int index;
	bool exists;
	bool exists2;
	int index2;
	int pointer = 0;
	dotObj tempModel, null;
	int l2, l1;
	int selection;
	int toErase;
	int toEraseIndex;



	std::cout << "Build Initial Trees" << std::endl;

	////======FIRST BRANCH====================
	selection = 0;
	tempModel = null;
	neighBours.clear();
	temp.clear();
	examined.clear();
	for (int i = 0; i < this->vertices.size(); i++){
		if (this->vertices[i][2]>maxz){
			maxz = this->vertices[i][2];
			maxzIndex = i;
		}
	}
	index = maxzIndex + 1;
	index2 = index;
	examined.push_back(index);
	tempModel.vertices.push_back(this->vertices[index2 - 1]);
	std::cout << index << std::endl;

	for (int i = 0; i < this->lines.size(); i++){
		if (this->lines[i][1] == index){
			neighBours.push_back(this->lines[i][0]);
			index2 = this->lines[i][0];
			tempModel.vertices.push_back(this->vertices[index2 - 1]);
			for (int g = 0; g < tempModel.vertices.size(); g++){
				if (tempModel.vertices[g] == this->vertices[index - 1]){ l1 = g + 1; }
				if (tempModel.vertices[g] == this->vertices[index2 - 1]){ l2 = g + 1; }
			}
			tempModel.lines.push_back({ l1, l2 });
			std::cout << this->lines[i][0] << std::endl;
		}
		else if (this->lines[i][0] == index){
			neighBours.push_back(this->lines[i][1]);
			index2 = this->lines[i][1];
			tempModel.vertices.push_back(this->vertices[index2 - 1]);
			for (int g = 0; g < tempModel.vertices.size(); g++){
				if (tempModel.vertices[g] == this->vertices[index - 1]){ l1 = g + 1; }
				if (tempModel.vertices[g] == this->vertices[index2 - 1]){ l2 = g + 1; }
			}
			tempModel.lines.push_back({ l1, l2 });
			std::cout << this->lines[i][1] << std::endl;
		}
	}

	for (int k = 0; k < 1; k++){
		temp = neighBours;
		neighBours.clear();
		for (int i = 0; i < temp.size(); i++){
			index = temp[i];
			examined.push_back(index);
			for (int i = 0; i < this->lines.size(); i++){
				if (this->lines[i][1] == index){
					exists = false;
					exists = this->functions.valueExistsInVector(examined, this->lines[i][0]);
					if (!exists){
						exists2 = false;
						exists2 = this->functions.valueExistsInVector(neighBours, this->lines[i][0]);
						if (!exists2){
							neighBours.push_back(this->lines[i][0]);
							index2 = this->lines[i][0];
							tempModel.vertices.push_back(this->vertices[index2 - 1]);
							for (int g = 0; g < tempModel.vertices.size(); g++){
								if (tempModel.vertices[g] == this->vertices[index - 1]){ l1 = g + 1; }
								if (tempModel.vertices[g] == this->vertices[index2 - 1]){ l2 = g + 1; }
							}
							tempModel.lines.push_back({ l1, l2 });
							std::cout << this->lines[i][0] << std::endl;
						}
					}
				}
				else if (this->lines[i][0] == index){
					exists = false;
					exists = this->functions.valueExistsInVector(examined, this->lines[i][1]);
					if (!exists){
						exists2 = false;
						exists2 = this->functions.valueExistsInVector(neighBours, this->lines[i][1]);
						if (!exists2){
							neighBours.push_back(this->lines[i][1]);
							index2 = this->lines[i][1];
							tempModel.vertices.push_back(this->vertices[index2 - 1]);
							for (int g = 0; g < tempModel.vertices.size(); g++){
								if (tempModel.vertices[g] == this->vertices[index - 1]){ l1 = g + 1; }
								if (tempModel.vertices[g] == this->vertices[index2 - 1]){ l2 = g + 1; }
							}
							tempModel.lines.push_back({ l1, l2 });
							std::cout << this->lines[i][1] << std::endl;
						}
					}
				}
			}
		}
	}

	toEraseIndex = selection;
	for (int g = 0; g < tempModel.vertices.size(); g++){
		if (tempModel.vertices[g] == this->vertices[neighBours[toEraseIndex] - 1]){ toErase = g; }
	}
	neighBours.erase(neighBours.begin() + toEraseIndex);

	tempModel.vertices.erase(tempModel.vertices.begin() + toErase);
	tempModel.vertices.erase(tempModel.vertices.begin());
	tempModel.lines.clear();
	tempModel.lines.push_back({ 1, 2 });

	for (int k = 0; k < gens; k++){
		temp = neighBours;
		neighBours.clear();
		for (int i = 0; i < temp.size(); i++){
			index = temp[i];
			examined.push_back(index);
			for (int i = 0; i < this->lines.size(); i++){
				if (this->lines[i][1] == index){
					exists = false;
					exists = this->functions.valueExistsInVector(examined, this->lines[i][0]);
					if (!exists){
						exists2 = false;
						exists2 = this->functions.valueExistsInVector(neighBours, this->lines[i][0]);
						if (!exists2){
							neighBours.push_back(this->lines[i][0]);
							index2 = this->lines[i][0];
							tempModel.vertices.push_back(this->vertices[index2 - 1]);
							for (int g = 0; g < tempModel.vertices.size(); g++){
								if (tempModel.vertices[g] == this->vertices[index - 1]){ l1 = g + 1; }
								if (tempModel.vertices[g] == this->vertices[index2 - 1]){ l2 = g + 1; }
							}
							tempModel.lines.push_back({ l1, l2 });
							std::cout << this->lines[i][0] << std::endl;
						}
					}
				}
				else if (this->lines[i][0] == index){
					exists = false;
					exists = this->functions.valueExistsInVector(examined, this->lines[i][1]);
					if (!exists){
						exists2 = false;
						exists2 = this->functions.valueExistsInVector(neighBours, this->lines[i][1]);
						if (!exists2){
							neighBours.push_back(this->lines[i][1]);
							index2 = this->lines[i][1];
							tempModel.vertices.push_back(this->vertices[index2 - 1]);
							for (int g = 0; g < tempModel.vertices.size(); g++){
								if (tempModel.vertices[g] == this->vertices[index - 1]){ l1 = g + 1; }
								if (tempModel.vertices[g] == this->vertices[index2 - 1]){ l2 = g + 1; }
							}
							tempModel.lines.push_back({ l1, l2 });
							std::cout << this->lines[i][1] << std::endl;
						}
					}
				}
			}
		}
	}
	mod1 = tempModel;

	//======SECOND BRANCH====================
	selection = 1;
	tempModel = null;
	neighBours.clear();
	temp.clear();
	examined.clear();
	for (int i = 0; i < this->vertices.size(); i++){
		if (this->vertices[i][2]>maxz){
			maxz = this->vertices[i][2];
			maxzIndex = i;
		}
	}
	index = maxzIndex + 1;
	index2 = index;
	examined.push_back(index);
	tempModel.vertices.push_back(this->vertices[index2 - 1]);
	std::cout << index << std::endl;
	for (int i = 0; i < this->lines.size(); i++){
		if (this->lines[i][1] == index){
			neighBours.push_back(this->lines[i][0]);
			index2 = this->lines[i][0];
			tempModel.vertices.push_back(this->vertices[index2 - 1]);
			for (int g = 0; g < tempModel.vertices.size(); g++){
				if (tempModel.vertices[g] == this->vertices[index - 1]){ l1 = g + 1; }
				if (tempModel.vertices[g] == this->vertices[index2 - 1]){ l2 = g + 1; }
			}
			tempModel.lines.push_back({ l1, l2 });
			std::cout << this->lines[i][0] << std::endl;
		}
		else if (this->lines[i][0] == index){
			neighBours.push_back(this->lines[i][1]);
			index2 = this->lines[i][1];
			tempModel.vertices.push_back(this->vertices[index2 - 1]);
			for (int g = 0; g < tempModel.vertices.size(); g++){
				if (tempModel.vertices[g] == this->vertices[index - 1]){ l1 = g + 1; }
				if (tempModel.vertices[g] == this->vertices[index2 - 1]){ l2 = g + 1; }
			}
			tempModel.lines.push_back({ l1, l2 });
			std::cout << this->lines[i][1] << std::endl;
		}
	}

	for (int k = 0; k < 1; k++){
		temp = neighBours;
		neighBours.clear();
		for (int i = 0; i < temp.size(); i++){
			index = temp[i];
			examined.push_back(index);
			for (int i = 0; i < this->lines.size(); i++){
				if (this->lines[i][1] == index){
					exists = false;
					exists = this->functions.valueExistsInVector(examined, this->lines[i][0]);
					if (!exists){
						exists2 = false;
						exists2 = this->functions.valueExistsInVector(neighBours, this->lines[i][0]);
						if (!exists2){
							neighBours.push_back(this->lines[i][0]);
							index2 = this->lines[i][0];
							tempModel.vertices.push_back(this->vertices[index2 - 1]);
							for (int g = 0; g < tempModel.vertices.size(); g++){
								if (tempModel.vertices[g] == this->vertices[index - 1]){ l1 = g + 1; }
								if (tempModel.vertices[g] == this->vertices[index2 - 1]){ l2 = g + 1; }
							}
							tempModel.lines.push_back({ l1, l2 });
							std::cout << this->lines[i][0] << std::endl;
						}
					}
				}
				else if (this->lines[i][0] == index){
					exists = false;
					exists = this->functions.valueExistsInVector(examined, this->lines[i][1]);
					if (!exists){
						exists2 = false;
						exists2 = this->functions.valueExistsInVector(neighBours, this->lines[i][1]);
						if (!exists2){
							neighBours.push_back(this->lines[i][1]);
							index2 = this->lines[i][1];
							tempModel.vertices.push_back(this->vertices[index2 - 1]);
							for (int g = 0; g < tempModel.vertices.size(); g++){
								if (tempModel.vertices[g] == this->vertices[index - 1]){ l1 = g + 1; }
								if (tempModel.vertices[g] == this->vertices[index2 - 1]){ l2 = g + 1; }
							}
							tempModel.lines.push_back({ l1, l2 });
							std::cout << this->lines[i][1] << std::endl;
						}
					}
				}
			}
		}
	}

	toEraseIndex = selection;
	for (int g = 0; g < tempModel.vertices.size(); g++){
		if (tempModel.vertices[g] == this->vertices[neighBours[toEraseIndex] - 1]){ toErase = g; }
	}
	neighBours.erase(neighBours.begin() + toEraseIndex);
	tempModel.vertices.erase(tempModel.vertices.begin() + toErase);
	tempModel.vertices.erase(tempModel.vertices.begin());
	tempModel.lines.clear();
	tempModel.lines.push_back({ 1, 2 });

	for (int k = 0; k < gens; k++){
		temp = neighBours;
		neighBours.clear();
		for (int i = 0; i < temp.size(); i++){
			index = temp[i];
			examined.push_back(index);
			for (int i = 0; i < this->lines.size(); i++){
				if (this->lines[i][1] == index){
					exists = false;
					exists = this->functions.valueExistsInVector(examined, this->lines[i][0]);
					if (!exists){
						exists2 = false;
						exists2 = this->functions.valueExistsInVector(neighBours, this->lines[i][0]);
						if (!exists2){
							neighBours.push_back(this->lines[i][0]);
							index2 = this->lines[i][0];
							tempModel.vertices.push_back(this->vertices[index2 - 1]);
							for (int g = 0; g < tempModel.vertices.size(); g++){
								if (tempModel.vertices[g] == this->vertices[index - 1]){ l1 = g + 1; }
								if (tempModel.vertices[g] == this->vertices[index2 - 1]){ l2 = g + 1; }
							}
							tempModel.lines.push_back({ l1, l2 });
							std::cout << this->lines[i][0] << std::endl;
						}
					}
				}
				else if (this->lines[i][0] == index){
					exists = false;
					exists = this->functions.valueExistsInVector(examined, this->lines[i][1]);
					if (!exists){
						exists2 = false;
						exists2 = this->functions.valueExistsInVector(neighBours, this->lines[i][1]);

						if (!exists2){
							neighBours.push_back(this->lines[i][1]);
							index2 = this->lines[i][1];
							tempModel.vertices.push_back(this->vertices[index2 - 1]);
							for (int g = 0; g < tempModel.vertices.size(); g++){
								if (tempModel.vertices[g] == this->vertices[index - 1]){ l1 = g + 1; }
								if (tempModel.vertices[g] == this->vertices[index2 - 1]){ l2 = g + 1; }
							}
							tempModel.lines.push_back({ l1, l2 });
							std::cout << this->lines[i][1] << std::endl;
						}
					}
				}
			}
		}
	}
	mod2 = tempModel;

	return;
}

