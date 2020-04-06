#include "lungModelling.h"



void dotObj::getNeighbouringVerticesOnSkeleton(void){

	this->neighbouringVerticesPerVertexOnSkeleton.clear();
	this->neighbouringLinesPerVertexOnSkeleton.clear();

	this->neighbouringVerticesPerVertexOnSkeleton.resize(this->vertices.size());
	this->neighbouringLinesPerVertexOnSkeleton.resize(this->vertices.size());

	for (int i = 0; i < this->lines.size(); i++){
		neighbouringVerticesPerVertexOnSkeleton[lines[i][0] - 1].push_back(lines[i][1] - 1);
		neighbouringLinesPerVertexOnSkeleton[lines[i][0] - 1].push_back(i);



		neighbouringVerticesPerVertexOnSkeleton[lines[i][1] - 1].push_back(lines[i][0] - 1);
		neighbouringLinesPerVertexOnSkeleton[lines[i][1] - 1].push_back(i);

		

	}

	return;
}


std::vector<std::vector<int> > dotObj::getNeighboursPerVertex(void){
	std::vector<std::vector<int> > neighboursPerVertex;
	neighboursPerVertex.resize(this->vertices.size());
	for (int i = 0; i < this->vertices.size(); i++){
		for (int j = 0; j < this->lines.size(); j++){
			if (this->lines[j][0] - 1 == i)
			{
				neighboursPerVertex[i].push_back(this->lines[j][1] - 1);
			}
			if (this->lines[j][1] - 1 == i)
			{
				neighboursPerVertex[i].push_back(this->lines[j][0] - 1);
			}
		}
	}

	return neighboursPerVertex;
}


void dotObj::removeVertexFromSkeleton(int vertexIndex){

	
	if (this->neighbouringVerticesPerVertexOnSkeleton.size() == 0){
		this->getNeighbouringVerticesOnSkeleton();
	}
	for (int j = 0; j < this->lines.size(); j++){
		if (this->lines[j][0] > vertexIndex + 1){ this->lines[j][0] --; }
		if (this->lines[j][1] > vertexIndex + 1){ this->lines[j][1] --; }
	}

	if (neighbouringVerticesPerVertexOnSkeleton[vertexIndex][0] > vertexIndex){ neighbouringVerticesPerVertexOnSkeleton[vertexIndex][0] --; }
	if (neighbouringVerticesPerVertexOnSkeleton[vertexIndex][1] > vertexIndex){ neighbouringVerticesPerVertexOnSkeleton[vertexIndex][1] --; }

	std::sort(neighbouringLinesPerVertexOnSkeleton[vertexIndex].begin(), neighbouringLinesPerVertexOnSkeleton[vertexIndex].end());
	this->lines.erase(this->lines.begin() + neighbouringLinesPerVertexOnSkeleton[vertexIndex][1]);
	this->lines.erase(this->lines.begin() + neighbouringLinesPerVertexOnSkeleton[vertexIndex][0]);

	this->lines.push_back({ neighbouringVerticesPerVertexOnSkeleton[vertexIndex][0] + 1, neighbouringVerticesPerVertexOnSkeleton[vertexIndex][1] + 1 });
	this->vertices.erase(this->vertices.begin() + vertexIndex);


	return;
}

void dotObj::removeVertexFromSkeleton(std::vector<int> toDelete){
	std::sort(toDelete.begin(), toDelete.end());
	toDelete.erase(std::unique(toDelete.begin(), toDelete.end()), toDelete.end());

	std::vector<int> linesToDelete;
	for (int j = toDelete.size() - 1; j >= 0; j--){
		for (int i = 0; i < this->lines.size(); i++){
			if ((lines[i][0] == toDelete[j] + 1) || (lines[i][1] == toDelete[j] + 1)){
				linesToDelete.push_back(i);
			}
		}
	}
	std::sort(linesToDelete.begin(), linesToDelete.end());
	linesToDelete.erase(std::unique(linesToDelete.begin(), linesToDelete.end()), linesToDelete.end());

	std::sort(linesToDelete.begin(), linesToDelete.end());
	for (int j = linesToDelete.size() - 1; j >= 0; j--){
		std::cout << "Remove Line : " << linesToDelete[j] << ":" << lines[linesToDelete[j]][0] << "-" << lines[linesToDelete[j]][1] << std::endl;
		this->lines.erase(this->lines.begin() + linesToDelete[j]);
	}



	for (int j = toDelete.size() - 1; j >= 0; j--){
		this->vertices.erase(this->vertices.begin() + toDelete[j]);

		for (int i = 0; i < this->lines.size(); i++){
			if (lines[i][0] > toDelete[j] + 1) {
				lines[i][0]--;
			}
			if (lines[i][1] > toDelete[j] + 1){
				lines[i][1]--;
			}
		}
	}

	


	return;
}



void dotObj::removeLineFromSkeleton(int vertexIndex){






}




std::vector<int> dotObj::getOneRingNeighboursOnSkeleton(int skeletonVertexIndex){
	if (this->neighbouringVerticesPerVertexOnSkeleton.size() == 0){
		this->getNeighbouringVerticesOnSkeleton();
	}

	std::vector<int> list;
	list = this->neighbouringVerticesPerVertexOnSkeleton[skeletonVertexIndex];

	return list;
}

std::vector<int> dotObj::getSecondRingNeighboursOnSkeleton(int skeletonVertexIndex){
	if (this->neighbouringVerticesPerVertexOnSkeleton.size() == 0){
		this->getNeighbouringVerticesOnSkeleton();
	}
	std::vector<int> listFirstOrder;
	std::vector<int> temp;
	std::vector<int> listSecondOrder;
	listFirstOrder = this->neighbouringVerticesPerVertexOnSkeleton[skeletonVertexIndex];
	int v;
	for (int i = 0; i < listFirstOrder.size(); i++){
		v = listFirstOrder[i];
		temp = this->neighbouringVerticesPerVertexOnSkeleton[v];
		for (int j = 0; j < temp.size(); j++){
			if (temp[j] != skeletonVertexIndex){
				listSecondOrder.push_back(temp[j]);
			}
		}
	}

	return listSecondOrder;
}






//Extract MFC skeleton from MFC and store to a dotObj object
dotObj dotObj::mcfskel(void)
{
	dotObj skel;
	Triangle_mesh tmesh;

	std::stringstream input = toOFF();
	//toOFF("input.off");
	//std::ifstream input("input.off");


	input >> tmesh;


	Skeleton skeleton;
	Skeletonization mcs(tmesh);
	//mcs.set_is_medially_centered(true);

	std::cout << "quality_speed_tradeoff : "<<mcs.quality_speed_tradeoff() << std::endl;
	std::cout << "is_medially_centered : " << mcs.is_medially_centered() << std::endl;
	std::cout << "medially_centered_speed_tradeoff : " << mcs.medially_centered_speed_tradeoff() << std::endl;

	mcs.set_quality_speed_tradeoff(0.4);
	mcs.set_min_edge_length(0.4);
	//mcs.set_is_medially_centered(false);


	// 1. Contract the mesh by mean curvature flow.
	//mcs.contract_geometry();
	// 2. Collapse short edges and split bad triangles.
	//mcs.collapse_edges();
	//mcs.split_faces();
	// 3. Fix degenerate vertices.
	//mcs.detect_degeneracies();
	// Perform the above three steps in one iteration.
	//mcs.contract();
	// Iteratively apply step 1 to 3 until convergence.
	mcs.contract_until_convergence();
	// Convert the contracted mesh into a curve skeleton and
	// get the correspondent surface points
	mcs.convert_to_skeleton(skeleton);
	std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
	std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";
	// Output all the edges of the skeleton.
	//std::ofstream output("skel.obj");

	int v0index, v1index;
	std::vector<float> v0, v1;

	BOOST_FOREACH(Skeleton_edge e, boost::edges(skeleton))
	{
		const Point& s = skeleton[source(e, skeleton)].point;
		const Point& t = skeleton[target(e, skeleton)].point;
		//output << "2 " << s << " " << t << "\n";

		v0 = { (float)s.x(), (float)s.y(), (float)s.z() };
		v1 = { (float)t.x(), (float)t.y(), (float)t.z() };

		bool v0exists = false;
		bool v1exists = false;

		for (int i = 0; i < skel.vertices.size(); i++){
			if (!v0exists){
				if (skel.vertices.at(i) == v0){
					v0index = i;
					v0exists = true;
				}
			}
			if (!v1exists){
				if (skel.vertices.at(i) == v1){
					v1index = i;
					v1exists = true;
				}
			}
		}

		if (!v0exists){
			skel.vertices.push_back(v0);
			v0index = skel.vertices.size() - 1;
		}

		if (!v1exists){
			skel.vertices.push_back(v1);
			v1index = skel.vertices.size() - 1;
		}
		skel.lines.push_back({ v0index + 1, v1index + 1 });
	}

	return skel;
}



std::vector<std::vector<int>> dotObj::pathFinder(std::vector<std::vector<int> > &neighboursPerVertex){
	std::vector<std::vector<int> > paths;
	std::vector<int>path, neighbourhood;
	int lastElement;
	bool found, doContinue;

	for (int i = 0; i < neighboursPerVertex.size(); i++){
		if (neighboursPerVertex[i].size()>2){
			for (int j = 0; j < neighboursPerVertex[i].size(); j++){
				paths.resize(paths.size() + 1);
				paths[paths.size() - 1].push_back(i);
				paths[paths.size() - 1].push_back(neighboursPerVertex[i][j]);
			}
		}
	}

	for (int i = 0; i < paths.size(); i++){
		path = paths[i]; //Get a path
		lastElement = path[path.size() - 1]; //Get the last element of a path
		neighbourhood = neighboursPerVertex[lastElement];
		if (neighbourhood.size() == 1){
		}
		if (neighbourhood.size() == 2){
			for (int temp = 0; temp < neighbourhood.size(); temp++){
				int  currentNeighbourhoodVertex = neighbourhood[temp];
				//-------------------------------------------------------------
				found = false;
				for (int temp2 = 0; temp2 < path.size(); temp2++){
					if (currentNeighbourhoodVertex == path[temp2]){
						found = true;
					}
				}
				//-------------------------------------------------------------
				if (!found){
					paths[i].push_back(currentNeighbourhoodVertex);
				}
			}
		}
		if (neighbourhood.size() == 3){
		}
	}

	doContinue = true;
	//for (int iterations = 0; iterations < 10; iterations++){
	while (doContinue){
		doContinue = false;
		for (int i = 0; i < paths.size(); i++){
			path = paths[i]; //Get a path
			lastElement = path[path.size() - 1]; //Get the last element of a path
			neighbourhood = neighboursPerVertex[lastElement];
			if (neighbourhood.size() == 1){
			}
			if (neighbourhood.size() == 2){
				for (int temp = 0; temp < neighbourhood.size(); temp++){
					int  currentNeighbourhoodVertex = neighbourhood[temp];
					//-------------------------------------------------------------
					found = false;
					for (int temp2 = 0; temp2 < path.size(); temp2++){
						if (currentNeighbourhoodVertex == path[temp2]){
							found = true;
						}
					}
					//-------------------------------------------------------------
					if (!found){
						paths[i].push_back(currentNeighbourhoodVertex);
						doContinue = true;
					}
				}
			}
			if (neighbourhood.size() == 3){
			}
		}
	}
	return paths;
}

std::vector<int> dotObj::getBrachesToDeleteBasedOnPathSize(std::vector<std::vector<int>> &paths, std::vector<std::vector<int>> &neighboursPerVertex, int lim){
	std::vector<int> toDelete;
	std::vector<int> path, neighbourhood;
	int vert;
	bool doProcessPath;
	//Find which to delete
	for (int i = 0; i < paths.size(); i++){
		doProcessPath = false;
		path = paths[i]; //Get a path
		if ((path.size() < lim)){
			for (int pathMemberIndex = 0; pathMemberIndex < path.size(); pathMemberIndex++){
				if (neighboursPerVertex[path[pathMemberIndex]].size() == 1)
				{
					doProcessPath = true;
				}
			}

			for (int pathMemberIndex = 1; pathMemberIndex < path.size(); pathMemberIndex++){
				if (neighboursPerVertex[path[pathMemberIndex]].size() >2)
				{
					//doProcessPath = false;
				}
			}

			if (doProcessPath){
				for (int pathMemberIndex = 1; pathMemberIndex < path.size(); pathMemberIndex++){
					toDelete.push_back(path[pathMemberIndex]);
				}
			}
		}
	}

	return toDelete;
}

std::vector<std::vector<int>> dotObj::analyzePathsPerVertex(void){
	//This function aims to build a bifurcating graph
	//This function receives as exports the graph of the 1-D graph and
	//gets the last element of each path. Assuming that the beginning of each path is
	//a 3-n vertex. If the last element is a 1-n vertex it is taken as a starting point
	//For each 3-n point get the connecting paths.

	std::vector<std::vector<int>> neighboursPerVertex;
	std::vector<std::vector<int>> paths;
	std::vector<std::vector<int>> result;
	std::vector<std::vector<int>> mesograph;
	std::vector<int> path, tempVar, tempVar2, r;

	int lastindex;
	int bifurcationPoint;
	int numOfNeighbours;
	neighboursPerVertex = getNeighboursPerVertex();
	paths = pathFinder(neighboursPerVertex);

	for (int i = 0; i < paths.size(); i++){
		path = paths[i];
		lastindex = path[path.size() - 1];
		numOfNeighbours = neighboursPerVertex[lastindex].size();
		if (numOfNeighbours == 1){
			tempVar = path;
			reverse(tempVar.begin(), tempVar.end());
			bifurcationPoint = tempVar[tempVar.size() - 1];
			for (int j = 0; j < paths.size(); j++){
				tempVar2 = paths[j];
				if ((tempVar2[0] == bifurcationPoint) && (tempVar2[tempVar2.size() - 1] != tempVar[0])){
					tempVar2.erase(tempVar2.begin());
					r = tempVar;
					r.insert(r.end(), tempVar2.begin(), tempVar2.end());
					mesograph.push_back(r);
				}
			}
		}
	}

	bool doContinue = true;

	while (doContinue){
		doContinue = false;

		result = mesograph;
		mesograph.clear();
		for (int i = 0; i < result.size(); i++){
			path = result[i];
			lastindex = path[path.size() - 1];
			numOfNeighbours = neighboursPerVertex[lastindex].size();

			if (numOfNeighbours == 3){
				doContinue = true;
				for (int j = 0; j < paths.size(); j++){
					tempVar2 = paths[j];
					if ((tempVar2[0] == lastindex) && (tempVar2[1] != path[path.size() - 2])) {
						tempVar2.erase(tempVar2.begin());
						r = path;
						r.insert(r.end(), tempVar2.begin(), tempVar2.end());
						mesograph.push_back(r);
					}
				}
			}
			if (numOfNeighbours == 1){
				mesograph.push_back(path);
			}
		}
	}

	result = mesograph;

	return result;
}

void dotObj::pathFinderTwoStep(std::vector<std::vector<int>> &neighboursPerVertex, std::vector<std::vector<int>> &paths){
	neighboursPerVertex = getNeighboursPerVertex();
	paths = pathFinder(neighboursPerVertex);
	return;
}

std::vector<std::vector<int>> dotObj::pathFinderTwoStep(void){
	std::vector<std::vector<int>> neighboursPerVertex;
	std::vector<std::vector<int>> paths;
	neighboursPerVertex = getNeighboursPerVertex();
	paths = pathFinder(neighboursPerVertex);
	return paths;
}
void dotObj::mcfskelRefineStepOne(void){
	std::vector<std::vector<int>> neighboursPerVertex;
	std::vector<std::vector<int>> paths;
	std::vector<int> toDelete;

	neighboursPerVertex = getNeighboursPerVertex();
	paths = pathFinder(neighboursPerVertex);
	toDelete = getBrachesToDeleteBasedOnPathSize(paths, neighboursPerVertex, 6);

	this->removeVertexFromSkeleton(toDelete);

	return;
}

void dotObj::mcfskelRefineStepOne(int lim){
	std::vector<std::vector<int>> neighboursPerVertex;
	std::vector<std::vector<int>> paths;
	std::vector<int> toDelete;

	neighboursPerVertex = getNeighboursPerVertex();
	paths = pathFinder(neighboursPerVertex);
	toDelete = getBrachesToDeleteBasedOnPathSize(paths, neighboursPerVertex, lim);
	this->removeVertexFromSkeleton(toDelete);
	
	return;
}




//Refine MFC Skeleton so that the skeleton consists only of nodes with 3 or 1 neighbours.
void dotObj::mcfskelRefine(void){
	std::cout << "Refine skeleton" << std::endl;
	for (int i = 0; i < this->vertices.size(); i++){
		this->getNeighbouringVerticesOnSkeleton();
		if ((neighbouringVerticesPerVertexOnSkeleton[i].size() == 2) && (neighbouringLinesPerVertexOnSkeleton[i].size() == 2))
		{
			
			this->removeVertexFromSkeleton(i);
			//std::vector<int> toDelete;
			//toDelete.resize(1);
			//toDelete[0] = i;
			//this->removeVertexFromSkeleton(toDelete);
			i = 0;
		}
	}

	for (int i = 0; i < this->lines.size(); i++){
		std::sort(this->lines[i].begin(), this->lines[i].end());
	
	}
	std::sort(this->lines.begin(), this->lines.end(),
		[](const std::vector<int>& a, const std::vector<int>& b) {
		return a[0] < b[0];
	});



	return;
}







/*
//Refine MFC Skeleton so that the skeleton consists only of nodes with 3 or 1 neighbours.
void dotObj::mcfskelRefine(void){

	for (int i = 0; i < this->vertices.size(); i++){

		std::vector<int> neighbouringLines;
		std::vector<int> neighbours;
		for (int j = 0; j <this->lines.size(); j++){
			if (this->lines[j][0] == i + 1)
			{

				neighbouringLines.push_back(j);
				neighbours.push_back(this->lines[j][1]);
			}
			if (this->lines[j][1] == i + 1)
			{

				neighbouringLines.push_back(j);
				neighbours.push_back(this->lines[j][0]);
			}
		}

		if ((neighbouringLines.size() == 2) && (neighbours.size() == 2))
		{
			for (int j = 0; j < this->lines.size(); j++){
				if (this->lines[j][0] > i + 1){ this->lines[j][0] --; }
				if (this->lines[j][1] > i + 1){ this->lines[j][1] --; }
			}

			if (neighbours[0] > i + 1){ neighbours[0] --; }
			if (neighbours[1] > i + 1){ neighbours[1] --; }

			std::sort(neighbouringLines.begin(), neighbouringLines.end());
			this->lines.erase(this->lines.begin() + neighbouringLines[1]);
			this->lines.erase(this->lines.begin() + neighbouringLines[0]);

			this->lines.push_back({ neighbours[0], neighbours[1] });
			this->vertices.erase(this->vertices.begin() + i);
			i = 0;
		}
		neighbouringLines.clear();
		neighbours.clear();
	}

	for (int i = 0; i < this->lines.size(); i++){
		int numOfNeighbours = 0;
		for (int j = 0; j < this->lines.size(); j++){
			if (this->lines[i][1] == this->lines[j][0]){ numOfNeighbours++; }
		}
		if (numOfNeighbours > 2){
			std::cout << "caution 3 neighbours" << std::endl;
		}
	}

	return;
}


*/



dotObj dotObj::mcfskel(bool contract, bool convertToSkeleton)  //<<==TODO
{
	dotObj skel;
	std::stringstream input = toOFF();
	Triangle_mesh tmesh;
	input >> tmesh;
	Skeleton skeleton;
	Skeletonization mcs(tmesh);


	// 1. Contract the mesh by mean curvature flow.
	//mcs.contract_geometry();
	// 2. Collapse short edges and split bad triangles.
	//mcs.collapse_edges();
	//mcs.split_faces();
	// 3. Fix degenerate vertices.
	//mcs.detect_degeneracies();
	// Perform the above three steps in one iteration.
	//mcs.contract();
	// Iteratively apply step 1 to 3 until convergence.

	if (contract){
		//mcsf.contract_until_convergence();
		for (int i = 0; i < 2; i++){
			std::cout << "Contract & Detect degeneracies : Iteration " << i << std::endl;
			mcs.contract_geometry();
			mcs.detect_degeneracies();
		}
	}

	if (!convertToSkeleton){
		std::vector<float> ContractedPoint;
		ContractedPoint.resize(3);
		CGAL::Mean_curvature_flow_skeletonization<Triangle_mesh>::Meso_skeleton j = mcs.meso_skeleton();
		std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(j) << "\n";
		for (CGAL::Mean_curvature_flow_skeletonization<Triangle_mesh>::Meso_skeleton::Vertex_iterator v2 = j.vertices_begin(); v2 != j.vertices_end(); ++v2){
			//std::cout << v2->point() << std::endl;
			ContractedPoint[0] = v2->point().x();
			ContractedPoint[1] = v2->point().y();
			ContractedPoint[2] = v2->point().z();
			skel.vertices.push_back(ContractedPoint);
		}
	}



	// Convert the contracted mesh into a curve skeleton and
	// get the correspondent surface points
	if (convertToSkeleton){
		mcs.contract_until_convergence();
		mcs.convert_to_skeleton(skeleton);
		std::cout << "Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
		std::cout << "Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";
		// Output all the edges of the skeleton.
		//std::ofstream output("skel.obj");

		int v0index, v1index;
		std::vector<float> v0, v1;

		BOOST_FOREACH(Skeleton_edge e, boost::edges(skeleton))
		{
			const Point& s = skeleton[source(e, skeleton)].point;
			const Point& t = skeleton[target(e, skeleton)].point;
			//output << "2 " << s << " " << t << "\n";

			v0 = { (float)s.x(), (float)s.y(), (float)s.z() };
			v1 = { (float)t.x(), (float)t.y(), (float)t.z() };

			bool v0exists = false;
			bool v1exists = false;

			for (int i = 0; i < skel.vertices.size(); i++){
				if (!v0exists){
					if (skel.vertices.at(i) == v0){
						v0index = i;
						v0exists = true;
					}
				}
				if (!v1exists){
					if (skel.vertices.at(i) == v1){
						v1index = i;
						v1exists = true;
					}
				}
			}

			if (!v0exists){
				skel.vertices.push_back(v0);
				v0index = skel.vertices.size() - 1;
			}

			if (!v1exists){
				skel.vertices.push_back(v1);
				v1index = skel.vertices.size() - 1;
			}
			skel.lines.push_back({ v0index + 1, v1index + 1 });
		}
	}
	return skel;
}









void dotObj::mcfskelRefineStepTwo(void){


	std::vector<int> indicesToExamine;
	int debug = 0;
	Vector3f p;
	int i = 0;
	std::vector<std::vector<int>> n;


	n = this->getNeighboursPerVertex();


	for (i = 0; i < n.size(); i++){
		if (n[i].size() >3){
			indicesToExamine.push_back(i);
		}
	}





	if (indicesToExamine.size() > 0){
		i = indicesToExamine[0];
		functions.printVectorToFile(indicesToExamine, "../../../data/indicesToExamine.txt");
		functions.printVectorToFile(n[i], "n_i.txt");
		std::cout << "Debug point:" << debug++ << std::endl;

		this->vertices.push_back(this->vertices[i]);
		int lastVertexIndex = this->vertices.size() - 1;
		std::cout << "Debug point:" << debug++ << std::endl;
		this->lines.push_back({ lastVertexIndex + 1, i + 1 });


		for (int k = 0; k < lines.size(); k++){
			if ((lines[k][0] == n[i][2] + 1) && (lines[k][1] == i + 1)){
				lines[k][1] = lastVertexIndex + 1;

			}
			if ((lines[k][0] == n[i][3] + 1) && (lines[k][1] == i + 1)){
				lines[k][1] = lastVertexIndex + 1;

			}
			if ((lines[k][1] == n[i][2] + 1) && (lines[k][0] == i + 1)){
				lines[k][0] = lastVertexIndex + 1;

			}
			if ((lines[k][1] == n[i][3] + 1) && (lines[k][0] == i + 1)){
				lines[k][0] = lastVertexIndex + 1;

			}
		}

	}
	return;
}



