#include "lungModelling.h"


std::string multiBodyObj::exportToFile(std::string name){

	std::ofstream outfile;
	std::string output = name+ ".obj";
	
	outfile.open(output);

	

	for (unsigned int j = 0; j < this->scene.size(); j++)
	{
		int verticesCounter = 0;
		int normalsCounter = 0;
		outfile << "o Object" + std::to_string(j) << std::endl;
		
		
		for (unsigned int i = 0; i < this->scene.at(j).vertices.size(); i++)
		{
			outfile << "v" << " " << this->scene.at(j).vertices.at(i).at(0) << " " << this->scene.at(j).vertices.at(i).at(1) << " " << this->scene.at(j).vertices.at(i).at(2) << std::endl;
			
		}

		for (unsigned int i = 0; i << this->scene.at(j).normals.size(); i++)
		{
			outfile << "vn" << " " << this->scene.at(j).normals.at(i).at(0) << " " << this->scene.at(j).normals.at(i).at(1) << " " << this->scene.at(j).normals.at(i).at(2) << std::endl;
			
		}
		if (j > 0){ 
			for (int k = j - 1; k>=0; k--){
				verticesCounter += this->scene.at(k).vertices.size();
				normalsCounter += this->scene.at(k).normals.size();
			}
		}


		
		
		
		for (unsigned int i = 0; i < this->scene.at(j).faces.size(); i++)
		{
			outfile << "f" << " " << this->scene.at(j).faces.at(i).at(0) + verticesCounter << "//" << this->scene.at(j).faces.at(i).at(2) + normalsCounter << " " << this->scene.at(j).faces.at(i).at(3) + verticesCounter << "//" << this->scene.at(j).faces.at(i).at(5) + normalsCounter << " " << this->scene.at(j).faces.at(i).at(6) + verticesCounter << "//" << this->scene.at(j).faces.at(i).at(8) + normalsCounter << std::endl;
		}
	

		




	}

	outfile.close();

	std::cout << "Multi OBJ Export Complete" << std::endl;

	return output;

}


void multiBodyObj::addObject(std::string partName, dotObj &input){
	
	this->scene.push_back(input);

	return;
}