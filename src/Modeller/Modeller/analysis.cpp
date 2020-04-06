#include "lungModelling.h"
#include <CGAL/IO/write_ply_points.h>
// Property map associating a facet with an integer as id to an
// element in a vector stored internally

void dotObj::segmentBySkeleton(void) {
	while (this->sdf_property_map.size() > 0) { this->sdf_property_map.pop_back(); }
	while (this->segment_property_map_per_vertex.size() > 0) { this->segment_property_map_per_vertex.pop_back(); }
	while (this->segment_property_map.size() > 0) { this->segment_property_map.pop_back(); }

	//toOFF("input.off");
	//std::ifstream input1("input.off");
	//std::ifstream input2("input.off");

	std::stringstream input1 = toOFF();
	std::stringstream input2 = toOFF();

	std::vector<float> temp;
	std::vector<float> temp2;

	std::cout << "Extract sdf values" << std::endl;
	// create and read Polyhedron----------------------------------------------------------------------------------
	Polyhedron mesh;
	input1 >> mesh;
	// assign id field for each facet--------------------------------------------------------------------------------
	std::size_t facet_id = 0;
	for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin();
		facet_it != mesh.facets_end(); ++facet_it, ++facet_id) {
		facet_it->id() = facet_id;
	}
	// create a property-map for SDF values-----------------------------------------------------------------------
	std::vector<double> sdf_values_2(mesh.size_of_facets());
	Facet_with_id_pmap<double> sdf_property_map_2(sdf_values_2);
	//Caution
	CGAL::sdf_values(mesh, sdf_property_map_2, 0.005, 25, true);//<============
	// access SDF values (with constant-complexity)
	for (Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it) {
		temp.push_back(sdf_property_map_2[facet_it]);
	}

	std::cout << "Extract sdf based on distances to skeleton" << std::endl;
	PolyhedronSimpleCartesianItemsWithID3 tmesh;
	input2 >> tmesh;
	// extract the skeleton
	Skeleton_PS skeleton;
	CGAL::extract_mean_curvature_flow_skeleton(tmesh, skeleton);
	// init the polyhedron simplex indices
	CGAL::set_halfedgeds_items_id(tmesh);
	//for each input vertex compute its distance to the skeleton
	std::vector<double> distances(num_vertices(tmesh));
	BOOST_FOREACH(Skeleton_vertex_PS v, boost::vertices(skeleton))
	{
		const Point& skel_pt = skeleton[v].point;
		BOOST_FOREACH(vertex_descriptor_PS mesh_v, skeleton[v].vertices)
		{
			const Point& mesh_pt = mesh_v->point();
			distances[mesh_v->id()] = std::sqrt(CGAL::squared_distance(skel_pt, mesh_pt));
		}
	}
	// create a property-map for sdf values
	std::vector<double> sdf_values(num_faces(tmesh));
	Facet_with_id_pmap_s<double> sdf_property_map(sdf_values);
	// compute sdf values with skeleton
	int counter = 0;
	BOOST_FOREACH(face_descriptor_PS f, CGAL::faces(tmesh))
	{
		double dist = 0;
		BOOST_FOREACH(halfedge_descriptor_PS hd, halfedges_around_face(halfedge(f, tmesh), tmesh))
			dist += distances[target(hd, tmesh)->id()];
		sdf_property_map[f] = dist / 3.;
	}

	std::cout << "Segmentation from sdf values" << std::endl;
	// create a property-map for segment-ids (it is an adaptor for this case)
	std::vector<std::size_t> segment_ids(num_faces(tmesh));
	Facet_with_id_pmap_s<std::size_t> segment_property_map(segment_ids);
	const std::size_t number_of_clusters = 4;       // use 4 clusters in soft clustering
	const double smoothing_lambda = 0.26;  // importance of surface features, suggested to be in-between [0,1]
	// segment the mesh using default parameters
	CGAL::segmentation_from_sdf_values(tmesh, sdf_property_map, segment_property_map, number_of_clusters, smoothing_lambda, false);
	std::cout << "Number of segments: " << CGAL::segmentation_from_sdf_values(tmesh, sdf_property_map, segment_property_map) << "\n";
	BOOST_FOREACH(face_descriptor_PS f, CGAL::faces(tmesh))
	{
		this->segment_property_map.push_back(segment_property_map[f]);
		this->sdf_property_map.push_back(sdf_property_map[f]);
	}
	//WritePerVertex
	for (int k = 0; k < this->vertices.size(); k++) {
		segment_property_map_per_vertex.push_back(0);
	}
	for (int k = 0; k < this->faces.size(); k++) {
		int v0, v1, v2, group;
		group = this->segment_property_map.at(k);
		v0 = this->faces.at(k).at(0) - 1;
		v1 = this->faces.at(k).at(3) - 1;
		v2 = this->faces.at(k).at(6) - 1;
		segment_property_map_per_vertex.at(v0) = group;
		segment_property_map_per_vertex.at(v1) = group;
		segment_property_map_per_vertex.at(v2) = group;
	}

	return;
}

void dotObj::segmentByGeneration(void) {
	//dotObj tt = *(this);
	//tt.smoothingSparse("taubin", "weighted", 100, -0.52, 0.5);
	//tt.exportToFile("../../../data/smoothed");
	//dotObj contractedGeometry = tt.mcfskel(true, false);
	//dotObj contractedGeometry = this->mcfskel(true, false);
	//contractedGeometry.exportToFile("../../../data/SkeletonizationTest");
	//dotObj skel = tt.mcfskel();

	dotObj skel = this->mcfskel();
	skel.exportToFile("skeletonInitial");
	skel.mcfskelRefineStepOne(1);
	skel.mcfskelRefineStepOne(2);
	skel.mcfskelRefineStepOne(3);
	skel.mcfskelRefineStepOne(8);
	skel.mcfskelRefineStepTwo();
	skel.exportToFile("skeletonInitialRef");
	//skel.mcfskelRefine();
	this->segmentByGeneration(skel);

	return;
}

void dotObj::segmentByGeneration(dotObj & skel) {
	this->graphBasedAnalysis(skel);

	/*
	dotObj * treePCExp;
	simulation * sim;
	treePCExp = new dotObj();
	sim = new simulation();
	sim->generate3DPointCloudFromTree(120, 0.7, skel, *treePCExp, this->mAnalysis.LocalDiameterPerSkeletonEdge, "ok");
	treePCExp->exportToFile("../../../data/tree-0.5");
	delete sim;*/

	//=============================================================
	//Graph 2 OBJ test
	ggraph * rebuilt;
	dotObj * test3;
	gnode * g;
	rebuilt = new ggraph();
	*rebuilt = this->mAnalysis.graph;
	test3 = new dotObj;
	rebuilt->init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0]->doStop = true;
	rebuilt->init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0]->doStop = true;
	rebuilt->ggraph2ModelV(test3);
	test3->exportToFile("trachea");
	//=============================================================
	rebuilt = new ggraph();
	*rebuilt = this->mAnalysis.graph;
	test3 = new dotObj;
	g = new gnode();
	g = rebuilt->init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr;
	if (g->isLeft) {
		g = rebuilt->init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr;
	}
	else {
		g = rebuilt->init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr;
	}
	rebuilt->init = g;
	rebuilt->ggraph2Model(test3, true);
	test3->exportToFile("left");
	//=============================================================

	rebuilt = new ggraph();
	*rebuilt = this->mAnalysis.graph;
	test3 = new dotObj;
	g = new gnode();
	g = rebuilt->init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr;
	if (g->isRight) {
		g = rebuilt->init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr;
	}
	else {
		g = rebuilt->init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr;
	}
	rebuilt->init = g;
	rebuilt->ggraph2Model(test3, true);
	test3->exportToFile("right");

	this->segmentationApply();
	this->updateSegmentProperties();

	return;
}

void dotObj::segmentationApply(void) {
	ggraph g = this->mAnalysis.graph;
	std::vector<int> input = this->mAnalysis.MeshVertex2SkeletonVertex;	//Mesh2SkeletonVertices

	this->vMap.map.clear();
	this->vMap.map.resize(this->vertices.size());

	this->fMap.map.clear();
	this->fMap.map.resize(this->faces.size());

	std::cout << "Segmentation mode is " << this->segmentationMode << std::endl;

	//v.resize(input.size());
	for (int i = 0; i < input.size(); i++) {
		for (int j = 0; j < this->mAnalysis.SkeletonVertex2Generation.size(); j++) {
			if (input[i] == j) {
				//this->vMap.map[i].generation = this->mAnalysis.SkeletonVertex2Generation[j]; // + 1 ;
				this->vMap.map[i].generation = g.nodes[input[i]].generation;
				this->vMap.map[i].segmentIndex = g.nodes[input[i]].segmentIndex;
				this->vMap.map[i].isLeft = g.nodes[input[i]].isLeft;
				this->vMap.map[i].isRight = g.nodes[input[i]].isRight;
			}
		}
	}

	for (int i = 0; i < this->faces.size(); i++) {
		this->fMap.map[i].generation = (this->vMap.map[this->faces[i][0] - 1].generation + this->vMap.map[this->faces[i][3] - 1].generation + this->vMap.map[this->faces[i][6] - 1].generation) / 3;
		this->fMap.map[i].segmentIndex = this->vMap.map[this->faces[i][0] - 1].segmentIndex;
		this->fMap.map[i].isLeft = this->vMap.map[this->faces[i][0] - 1].isLeft;
		this->fMap.map[i].isRight = this->vMap.map[this->faces[i][0] - 1].isRight;
	}

#if TRUE
	if (this->segmentationMode == "NORMALS") {
		std::cout << "Detecting outlets based on normals" << std::endl;

		for (int x = 0; x < g.nodes.size(); x++) {
			if (g.nodes[x].isTerminal) {
				Vector3f direction = g.nodes[x].position - g.nodes[x].previousNodePtr->position;
				direction.normalize();
				for (int i = 0; i < faces.size(); i++) {
					Vector3f a = Vector3f(vertices[faces[i][0] - 1][0], vertices[faces[i][0] - 1][1], vertices[faces[i][0] - 1][2]);
					Vector3f b = Vector3f(vertices[faces[i][3] - 1][0], vertices[faces[i][3] - 1][1], vertices[faces[i][3] - 1][2]);
					Vector3f c = Vector3f(vertices[faces[i][6] - 1][0], vertices[faces[i][6] - 1][1], vertices[faces[i][6] - 1][2]);
					Vector3f n = (a - b).cross(b - c);
					n.normalize();

					if ((g.nodes[x].position - a).norm() < 10.0 && (n.dot(direction) > 0.99990)) {
						//this->fMap.map[i].isOutlet = true;
						//this->fMap.map[i].generation = 10;

						for (int x = 0; x < g.nodes.size(); x++) {
							if (g.nodes[x].isTerminal) {
								if ((g.nodes[x].position - g.nodes[input[faces[i][0] - 1]].position).norm() < 3.0) {
									this->fMap.map[i].isOutlet = true;
									this->fMap.map[i].generation = 10;
									this->fMap.map[i].segmentIndex = this->vMap.map[faces[i][0] - 1].segmentIndex;
									this->fMap.map[i].isLeft = this->vMap.map[faces[i][0] - 1].isLeft;
									this->fMap.map[i].isRight = this->vMap.map[faces[i][0] - 1].isRight;
								}
							}
						}
					}
				}
			}

			if (g.nodes[x].isInlet) {
				Vector3f direction = g.nodes[x].position - g.nodes[x].nextNodesPtr[0]->position;
				direction.normalize();
				for (int i = 0; i < faces.size(); i++) {
					Vector3f a = Vector3f(vertices[faces[i][0] - 1][0], vertices[faces[i][0] - 1][1], vertices[faces[i][0] - 1][2]);
					Vector3f b = Vector3f(vertices[faces[i][3] - 1][0], vertices[faces[i][3] - 1][1], vertices[faces[i][3] - 1][2]);
					Vector3f c = Vector3f(vertices[faces[i][6] - 1][0], vertices[faces[i][6] - 1][1], vertices[faces[i][6] - 1][2]);
					Vector3f n = (a - b).cross(b - c);
					n.normalize();

					if ((g.nodes[x].position - a).norm() < 10.0 && (n.dot(direction) > 0.99990)) {
						for (int x = 0; x < g.nodes.size(); x++) {
							if (g.nodes[x].isInlet) {
								if ((g.nodes[x].position - g.nodes[input[faces[i][0] - 1]].position).norm() < 8.0) {
									this->fMap.map[i].isInlet = true;
									this->fMap.map[i].generation = 20;
									this->fMap.map[i].segmentIndex = this->vMap.map[faces[i][0] - 1].segmentIndex;
								}
							}
						}
					}
				}
			}
		}
	}
#endif

#if TRUE

	if (this->segmentationMode != "NORMALS") {
		std::cout << "Detecting outlets based on SDF" << std::endl;
		std::cout << "Shape Diameter Function 0.01 rad" << std::endl;
		this->getSDF(0.01, 20, false, false, 4, 0.2, false);
		std::vector<float> sdf1 = this->sdf_property_map;
		std::cout << "Shape Diameter Function 1.50 rad" << std::endl;
		this->getSDF(2.5, 20, false, false, 4, 0.2, false);
		std::vector<float> sdf2 = this->sdf_property_map;

		std::cout << "SDF Complete" << std::endl;

		std::vector<bool> v4;
		v4.resize(faces.size());
		for (int i = 0; i < faces.size(); i++) {
			v4[i] = false;
		}
		for (int i = 0; i < faces.size(); i++) {
			if (sdf1[i] / sdf2[i] > 3.5) {
				//if (abs(sdf1[i] - sdf2[i])>0.5){
				v4[i] = true;
				/*int v0 = this->faces.at(i).at(0) - 1;
				int v1 = this->faces.at(i).at(3) - 1;
				int v2 = this->faces.at(i).at(6) - 1;
				v4[v0] = true;
				v4[v1] = true;
				v4[v2] = true;*/
			}
		}

		std::cout << "SDF Complete" << std::endl;

		for (int x = 0; x < g.nodes.size(); x++) {
			if (g.nodes[x].isTerminal) {
				for (int i = 0; i < faces.size(); i++) {
					if (v4[i]) {
						for (int x = 0; x < g.nodes.size(); x++) {
							if (g.nodes[x].isTerminal) {
								if ((g.nodes[x].position - g.nodes[input[faces[i][0] - 1]].position).norm() < 3.0) {
									this->fMap.map[i].isOutlet = true;
									this->fMap.map[i].generation = 10;
									this->fMap.map[i].segmentIndex = this->vMap.map[faces[i][0] - 1].segmentIndex;
									this->fMap.map[i].isLeft = this->vMap.map[faces[i][0] - 1].isLeft;
									this->fMap.map[i].isRight = this->vMap.map[faces[i][0] - 1].isRight;
								}
							}
						}
					}
				}
			}

			if (g.nodes[x].isInlet) {
				for (int i = 0; i < faces.size(); i++) {
					if (v4[i]) {
						for (int x = 0; x < g.nodes.size(); x++) {
							if (g.nodes[x].isInlet) {
								if ((g.nodes[x].position - g.nodes[input[faces[i][0] - 1]].position).norm() < 8.0) {
									this->fMap.map[i].isInlet = true;
									this->fMap.map[i].generation = 20;
									this->fMap.map[i].segmentIndex = this->vMap.map[faces[i][0] - 1].segmentIndex;
								}
							}
						}
					}
				}
			}
		}
	}

	/*for (int i = 0; i < input.size(); i++){
		if (v4[i]){
		v[i] = 10;
		}
		for (int x = 0; x < g.nodes.size(); x++){
		if (g.nodes[x].isTerminal){
		if ((g.nodes[x].position - g.nodes[input[i]].position).norm() < 3.0){
		if (v4[i]){
		v[i] = 10;
		}
		}
		}
		if (g.nodes[x].isInlet){
		if ((g.nodes[x].position - g.nodes[input[i]].position).norm() < 10.0){
		if (v4[i]){
		v[i] = 20;
		}
		}
		}
		}
		if (mAnalysis.neighboursPerVertex[input[i]].size() == 1){
		if ((v4[i]) && input[i] != this->mAnalysis.inlet){ v[i] = 10; }
		if ((v4[i]) && input[i] == this->mAnalysis.inlet){ v[i] = 20; }
		}
		}*/

#endif

	return;
}

void dotObj::updateSegmentProperties(void) {
	this->segment_property_map_per_vertex.clear();
	this->segment_property_map_per_vertex.resize(this->vertices.size());

	this->segment_property_map.clear();
	this->segment_property_map.resize(this->faces.size());

	this->segment_property_name.clear();
	this->segment_property_name.resize(this->faces.size());

	for (int i = 0; i < this->faces.size(); i++) {
		this->segment_property_map[i] = this->fMap.map[i].generation * 1000 + this->fMap.map[i].segmentIndex;

		std::string LorR = "";
		if (this->fMap.map[i].isLeft)LorR = "L";
		if (this->fMap.map[i].isRight)LorR = "R";

		if ((this->fMap.map[i].generation != 10) && (this->fMap.map[i].generation != 20))
			this->segment_property_name[i] = "gen_" + std::to_string(this->fMap.map[i].generation - 1) + "_" + std::to_string(this->fMap.map[i].segmentIndex) + "_" + LorR;
		if ((this->fMap.map[i].generation == 10))
			this->segment_property_name[i] = "outlet_" + std::to_string(this->fMap.map[i].segmentIndex) + "_" + LorR;
		if ((this->fMap.map[i].generation == 20))
			this->segment_property_name[i] = "inlet";
	}

	if (this->segment_property_map.size() > 0) {
		for (int k = 0; k < this->segment_property_map_per_vertex.size(); k++) {
			this->segment_property_map_per_vertex[k] = 0;
		}

		for (int k = 0; k < this->faces.size(); k++) {
			int v0, v1, v2, group;
			group = this->segment_property_map.at(k);
			v0 = this->faces.at(k).at(0) - 1;
			v1 = this->faces.at(k).at(3) - 1;
			v2 = this->faces.at(k).at(6) - 1;
			this->segment_property_map_per_vertex.at(v0) = group;
			this->segment_property_map_per_vertex.at(v1) = group;
			this->segment_property_map_per_vertex.at(v2) = group;
		}
	}

	return;
}

void dotObj::graphBased1DModelAnalysis(bool LR, bool generations,int inlet) {
	//this->splitSkeletonEdges(0.4);
	this->mAnalysis.parseTerminals = true;
	std::cout << "Locate neighbours" << std::endl;
	this->mAnalysis.neighboursPerVertex = this->getNeighboursPerVertex();
	std::cout << "Find paths" << std::endl;
	this->mAnalysis.pathFinder(this->mAnalysis.neighboursPerVertex);
	distance d = this->mAnalysis.getMaxPath();
	if (inlet < 0) {
		this->mAnalysis.inlet = this->mAnalysis.paths[d.index].nodes[mAnalysis.paths[d.index].nodes.size() - 1];
	}else {
		this->mAnalysis.inlet = inlet;
	}


	//std::cout << "Analyze paths" << std::endl;
	//this->mAnalysis.pathAnalyzer(*this);
	//this->mAnalysis.inlet = this->mAnalysis.fullpaths[0][0];
	//this->mAnalysis.pathAnalyzer(*this);

	std::cout << "Generate graph" << std::endl;
	this->mAnalysis.generateGraph(*this, this->mAnalysis.inlet);

	this->mAnalysis.parseTerminals = true;
	std::cout << "Analyze graph" << std::endl;
	this->mAnalysis.analyzeGraph();
	std::cout << "Generate graph features" << std::endl;
	if (LR) {
		this->mAnalysis.generateGraphFeaturesLRDiscrimination();
	}
	if (generations) {
		this->mAnalysis.generateGraphFeaturesGeneration();
	}
	return;
}

void dotObj::graphBasedPruning(void) {
	this->mAnalysis.parseTerminals = true;
	this->mAnalysis.neighboursPerVertex = this->getNeighboursPerVertex();

	this->mAnalysis.pathFinder(this->mAnalysis.neighboursPerVertex);
	distance d = this->mAnalysis.getMaxPath();
	this->mAnalysis.inlet = this->mAnalysis.paths[d.index].nodes[mAnalysis.paths[d.index].nodes.size() - 1];
	this->mAnalysis.pathAnalyzer(*this);
	this->mAnalysis.inlet = this->mAnalysis.fullpaths[0][0];
	this->mAnalysis.pathAnalyzer(*this);
	this->mAnalysis.generateGraph(*this, this->mAnalysis.inlet);
	this->mAnalysis.parseTerminals = true;
	this->mAnalysis.analyzeGraph();

	bool remove = false;
	for (int i = 0; i < this->mAnalysis.graph.nodes.size(); i++) {
		if (this->mAnalysis.graph.nodes[i].isTerminal) {
			gnode * g;
			remove = false;
			g = this->mAnalysis.graph.nodes[i].previousBifurcationPtr;
			for (int j = 0; j < g->nextBifurcationPtr.size(); j++) {
				if (!g->nextBifurcationPtr[j]->isTerminal) {
					remove = true;
				}
			}
			if (remove) {
				for (int j = g->nextNodesPtr.size() - 1; j >= 0; j--) {
					if (g->nextNodesPtr[j]->nextBifurcationPtr.size() > 0) {
						if (g->nextNodesPtr[j]->nextBifurcationPtr[0]->isTerminal) {
							g->nextNodesPtr.erase(g->nextNodesPtr.begin() + j);
						}
					}
				}
			}
		}
	}
	return;
}

std::vector<float> dotObj::locate_flat_tips_based_on_calculated_normal(bool doFaces = true) {
	std::vector<int> adjacentFaces;
	std::vector<float> sdf_result;
	std::vector<Vector3f> adjacentFacesNormals;
	float sum;
	Vector3f n, f;

	sdf_result.resize(this->faces.size());

#pragma omp parallel  shared(sdf_result)
#pragma omp for
	for (int i = 0; i < this->segment_property_map.size(); i++) {
		sdf_result[i] = 0.5;
	}

	if (doFaces) {
#pragma omp parallel  shared(sdf_result) private(f,n,sum,adjacentFaces)
#pragma omp for
		for (int i = 0; i < this->faces.size(); i++) {
			adjacentFaces = findAdjacentFaces(i);
			f = calculateNormal(i);
			sum = 0;
			for (int j = 0; j < adjacentFaces.size(); j++) {
				n = calculateNormal(adjacentFaces[j]);
				sum += functions.cosineAngleBetweenVectors(f, n);
			}
			sum = sum / adjacentFaces.size();
			if (sum == 1.0) {
				for (int k = 0; k < adjacentFaces.size(); k++) {
					sdf_result[adjacentFaces[k]] = 1.0;
				}
			}
		}
	}
	if (!doFaces) {
		for (int i = 0; i < this->vertices.size(); i++) {
			adjacentFaces = findAdjacentFacesPerVertex(i);
			sum = 0;
			f = calculateNormal(adjacentFaces[0]);
			for (int j = 0; j < adjacentFaces.size(); j++) {
				n = calculateNormal(adjacentFaces[j]);
				sum += functions.cosineAngleBetweenVectors(f, n);
			}
			sum = sum / adjacentFaces.size();
			if (sum == 1.0) {
				//std::cout << "Vector: " << this->vertices[i][0] << "-" << this->vertices[i][1] << "-" << this->vertices[i][2] << std::endl;
				for (int k = 0; k < adjacentFaces.size(); k++) {
					sdf_result[adjacentFaces[k]] = 101.0;
				}
			}
		}
	}

	return sdf_result;
}

//Find edges
//Read faces in the dotObj and return edges
//This function stores data within the object
void dotObj::getEdges(void) {
	this->edges.clear();
	for (unsigned int i = 0; i < this->faces.size(); i++) {
		this->edges.push_back({ (int)this->faces.at(i).at(0) - 1, (int)this->faces.at(i).at(3) - 1 });
		this->edges.push_back({ (int)this->faces.at(i).at(3) - 1, (int)this->faces.at(i).at(6) - 1 });
		this->edges.push_back({ (int)this->faces.at(i).at(6) - 1, (int)this->faces.at(i).at(0) - 1 });
	}
	return;
}

analysis::analysis(void) {
	this->parseTerminals = true;
}

//Initializations
//
void analysis::initialize(dotObj * geom, dotObj * skel) {
	this->MeshVertex2SkeletonVertex.resize(geom->vertices.size());
	this->MeshVertex2SkeletonEdge.resize(geom->vertices.size());
	this->SkeletonEdge2MeshVertices.resize(skel->lines.size());
	this->LocalDiameterPerMeshVertex.resize(geom->vertices.size());
	this->LocalDiameterPerSkeletonEdge.resize(skel->lines.size());

	return;
}

void analysis::pathFinder(std::vector<std::vector<int>> &neighboursPerVertex) {
	std::vector<int>neighbourhood;
	gpath path;

	int lastElement;
	bool found, doContinue;

	for (int i = 0; i < neighboursPerVertex.size(); i++) {
		if (neighboursPerVertex[i].size() > 2) {
			for (int j = 0; j < neighboursPerVertex[i].size(); j++) {
				this->paths.resize(paths.size() + 1);
				this->paths[this->paths.size() - 1].nodes.push_back(i);
				this->paths[this->paths.size() - 1].nodes.push_back(neighboursPerVertex[i][j]);
			}
		}
	}

	for (int i = 0; i < this->paths.size(); i++) {
		path = this->paths[i]; //Get a path
		lastElement = path.nodes[path.nodes.size() - 1]; //Get the last element of a path
		neighbourhood = neighboursPerVertex[lastElement];
		if (neighbourhood.size() == 1) {
		}
		if (neighbourhood.size() == 2) {
			for (int j = 0; j < neighbourhood.size(); j++) {
				found = false;
				for (int k = 0; k < path.nodes.size(); k++) {
					if (neighbourhood[j] == path.nodes[k]) {
						found = true;
					}
				}
				if (!found) {
					paths[i].nodes.push_back(neighbourhood[j]);
				}
			}
		}
		if (neighbourhood.size() == 3) {
		}
	}

	doContinue = true;

	while (doContinue) {
		doContinue = false;
		for (int i = 0; i < paths.size(); i++) {
			path = paths[i]; //Get a path
			lastElement = path.nodes[path.nodes.size() - 1]; //Get the last element of a path
			neighbourhood = neighboursPerVertex[lastElement];
			if (neighbourhood.size() == 1) {
			}
			if (neighbourhood.size() == 2) {
				for (int j = 0; j < neighbourhood.size(); j++) {
					found = false;
					for (int k = 0; k < path.nodes.size(); k++) {
						if (neighbourhood[j] == path.nodes[k]) {
							found = true;
						}
					}
					//-------------------------------------------------------------
					if (!found) {
						paths[i].nodes.push_back(neighbourhood[j]);
						doContinue = true;
					}
				}
			}
			if (neighbourhood.size() == 3) {
			}
		}
	}

	this->branches.clear();
	for (int i = 0; i < this->paths.size(); i++) {
		this->branches.push_back(this->paths[i].nodes);
	}
	return;
}

void analysis::pathAnalyzer(dotObj & skel) {
	//This function aims to build a bifurcating graph
	//This function receives as exports the graph of the 1-D graph and
	//gets the last element of each path. Assuming that the beginning of each path is
	//a 3-n vertex. If the last element is a 1-n vertex it is taken as a starting point
	//For each 3-n point get the connecting paths.

	std::vector<std::vector<int>> m, q;
	std::vector<int> p, t1, t2, r;

	int lastindex;
	int bifurcationPoint;
	int numOfNeighbours;

	this->neighboursPerVertex.clear();
	this->paths.clear();
	this->branches.clear();
	this->fullpaths.clear();
	this->fullpathsSFInlet.clear();

	this->neighboursPerVertex = skel.getNeighboursPerVertex();
	this->pathFinder(neighboursPerVertex);

	for (int i = 0; i < paths.size(); i++) {
		p = paths[i].nodes;
		lastindex = p[p.size() - 1];
		numOfNeighbours = neighboursPerVertex[lastindex].size();
		if (numOfNeighbours == 1) {
			t1 = p;
			reverse(t1.begin(), t1.end());
			bifurcationPoint = t1[t1.size() - 1];
			for (int j = 0; j < paths.size(); j++) {
				t2 = paths[j].nodes;
				if ((t2[0] == bifurcationPoint) && (t2[t2.size() - 1] != t1[0])) {
					t2.erase(t2.begin());
					r = t1;
					r.insert(r.end(), t2.begin(), t2.end());
					m.push_back(r);
				}
			}
		}
	}

	bool doContinue = true;
	while (doContinue) {
		doContinue = false;

		q = m;
		m.clear();
		for (int i = 0; i < q.size(); i++) {
			p = q[i];
			lastindex = p[p.size() - 1];
			numOfNeighbours = neighboursPerVertex[lastindex].size();

			if (numOfNeighbours == 3) {
				doContinue = true;
				for (int j = 0; j < paths.size(); j++) {
					t2 = paths[j].nodes;
					if ((t2[0] == lastindex) && (t2[1] != p[p.size() - 2])) {
						t2.erase(t2.begin());
						r = p;
						r.insert(r.end(), t2.begin(), t2.end());
						m.push_back(r);
					}
				}
			}
			if (numOfNeighbours == 1) {
				m.push_back(p);
			}
		}
	}

	this->fullpaths = m;
	///Possible Errors
	this->fullpathsSFInlet.clear();
	for (int i = 0; i < this->fullpaths.size(); i++) {
		if (this->fullpaths[i][0] == this->inlet) {
			this->fullpathsSFInlet.push_back(this->fullpaths[i]);
		}
	}
	return;
}

void dotObj::graphBasedAnalysis(dotObj & skel) {
	dotObj * geom = (this);
	skel.getNeighbouringVerticesOnSkeleton();
	analysis a;
	this->mAnalysis = a;
	this->mAnalysis.initialize(geom, &skel);
	this->mAnalysis.model = (this);
	this->mAnalysis.skel = &skel;

	std::cout << "Mesh 2 Skeleton Vertices" << std::endl;
	for (int i = 0; i < geom->vertices.size(); i++) {
		Vector3f x = Vector3f(geom->vertices[i][0], geom->vertices[i][1], geom->vertices[i][2]);
		distance va = findMinimumDistance(geom->vertices[i], skel.vertices);
		distance vb = findMinimumDistance(geom->vertices[i], skel.vertices, skel.neighbouringVerticesPerVertexOnSkeleton[va.index]);
		int selectedLine = skel.neighbouringLinesPerVertexOnSkeleton[va.index][vb.nestedIndex];
		Vector3f p_a = Vector3f(skel.vertices[va.index][0], skel.vertices[va.index][1], skel.vertices[va.index][2]);
		Vector3f p_b = Vector3f(skel.vertices[vb.index][0], skel.vertices[vb.index][1], skel.vertices[vb.index][2]);
		Vector3f proj = functions.projPointOnVector(x, p_a, p_b);

		float projectionDistance;
		projectionDistance = (proj - x).norm();
		projectionDistance = (p_a - x).norm();
		this->mAnalysis.MeshVertex2SkeletonVertex[i] = va.index;
		this->mAnalysis.MeshVertex2SkeletonEdge[i] = selectedLine;
		this->mAnalysis.LocalDiameterPerMeshVertex[i] = projectionDistance;
		this->mAnalysis.SkeletonEdge2MeshVertices[selectedLine].push_back(i);
	}

	//SDF Based
	std::cout << "Before SDF" << std::endl;
	geom->sdf_property_map = geom->getSDFPropertyMap(1, 25, false);  //(4, 25, true, 4, 0.0001);
	std::cout << "After SDF" << std::endl;
	geom->sdf_property_map_per_vertex.resize(geom->vertices.size());
	for (int jk = 0; jk < geom->sdf_property_map_per_vertex.size(); jk++) {
		geom->sdf_property_map_per_vertex[jk] = 0.0;
	}
	std::cout << "After SDF" << std::endl;
	for (int jk = 0; jk < geom->sdf_property_map.size(); jk++) {
		geom->sdf_property_map_per_vertex[(geom->faces[jk][0] - 1)] = geom->sdf_property_map[jk];
		geom->sdf_property_map_per_vertex[(geom->faces[jk][3] - 1)] = geom->sdf_property_map[jk];
		geom->sdf_property_map_per_vertex[(geom->faces[jk][6] - 1)] = geom->sdf_property_map[jk];
		//geom->faces[jk][0];
		//geom->faces[jk][3];
		//geom->faces[jk][6];
	}
	this->mAnalysis.LocalDiameterPerMeshVertex = geom->sdf_property_map_per_vertex;

	for (int i = 0; i < this->mAnalysis.LocalDiameterPerSkeletonEdge.size(); i++) {
		this->mAnalysis.LocalDiameterPerSkeletonEdge[i] = 0.0;
	}
	for (int i = 0; i < this->mAnalysis.SkeletonEdge2MeshVertices.size(); i++) {
		for (int j = 0; j < this->mAnalysis.SkeletonEdge2MeshVertices[i].size(); j++) {
			this->mAnalysis.LocalDiameterPerSkeletonEdge[i] += this->mAnalysis.LocalDiameterPerMeshVertex[this->mAnalysis.SkeletonEdge2MeshVertices[i][j]];
		}
		this->mAnalysis.LocalDiameterPerSkeletonEdge[i] = this->mAnalysis.LocalDiameterPerSkeletonEdge[i] / this->mAnalysis.SkeletonEdge2MeshVertices[i].size();
	}

	this->mAnalysis.neighboursPerVertex = skel.getNeighboursPerVertex();
	this->mAnalysis.pathFinder(this->mAnalysis.neighboursPerVertex);
	distance d = this->mAnalysis.getMaxPath();
	this->mAnalysis.inlet = this->mAnalysis.paths[d.index].nodes[this->mAnalysis.paths[d.index].nodes.size() - 1];
	int inlet = this->mAnalysis.inlet;
	std::cout << "Inlet is:" << inlet << std::endl;
	this->mAnalysis.pathAnalyzer(skel);
	this->mAnalysis.getSkeletonVertex2Generation(skel);
	std::cout << "Export skeleton generations" << std::endl;

	this->mAnalysis.generateGraph(skel, inlet);
	this->mAnalysis.skel = &skel;
	this->mAnalysis.model = (this);

	this->mAnalysis.analyzeGraph();
	this->mAnalysis.generateGraphFeatures();
	std::cout << "Graph based analysis complete" << std::endl;
	return;
}

void analysis::generateGraph(dotObj &skel, int inlet) {
	std::cout << "Generate graph" << std::endl;
	if (this->neighboursPerVertex.empty()) {
		this->neighboursPerVertex = skel.getNeighboursPerVertex();
	}
	std::vector<int> neighbours, buffer;;
	std::vector<int> checked;

	buffer.push_back(inlet);
	this->graph.nodes.resize(skel.vertices.size());

	while (!buffer.empty()) {
		int v = buffer.back();
		buffer.pop_back();
		gnode * node = new gnode();
		node->index = v;
		node->position = Vector3f(skel.vertices[v][0], skel.vertices[v][1], skel.vertices[v][2]);
		checked.push_back(v);
		neighbours = this->neighboursPerVertex[v];

		for (int i = 0; i < neighbours.size(); i++) {
			if (skel.functions.valueExistsInVector(checked, neighbours[i])) {
				node->previous = neighbours[i];
				node->previousNodePtr = &this->graph.nodes[node->previous];
			}
			else {
				node->next.push_back(neighbours[i]);
				node->nextNodesPtr.push_back(&this->graph.nodes[node->next.back()]);
				buffer.push_back(neighbours[i]);
			}
		}
		this->graph.nodes[v] = (*node);
		delete node;
	}

	this->graph.init = &this->graph.nodes[inlet];

	return;
}

void analysis::analyzeGraph(void) {
	//This function aims to build a bifurcating graph
	//This function receives as exports the graph of the 1-D graph and
	//gets the last element of each path. Assuming that the beginning of each path is
	//a 3-n vertex. If the last element is a 1-n vertex it is taken as a starting point
	//For each 3-n point get the connecting paths.

	std::vector<std::vector<int>> m, q;
	std::vector<int> p, t1, t2, r;

	int numOfNeighbours;

	for (int i = 0; i < paths.size(); i++) {
		p = paths[i].nodes;
		int lastindex = p[p.size() - 1];
		numOfNeighbours = neighboursPerVertex[lastindex].size();
		if ((numOfNeighbours == 1) && (lastindex == inlet)) {
			t1 = p;
			reverse(t1.begin(), t1.end());

			this->graph.nodes[inlet].isInlet = true;
			this->inletPtr = &this->graph.nodes[inlet];
			this->graph.init = &this->graph.nodes[inlet];

			this->graph.nodes[t1.back()].isBifurcation = true;
			for (int k = 0; k < t1.size() - 1; k++) {
				this->graph.nodes[t1[k]].nextBifurcation.push_back(t1.back());
				this->graph.nodes[t1[k]].nextBifurcationPtr.push_back(&this->graph.nodes[t1.back()]);
				this->graph.nodes[t1[k]].nextMedianPtr = &this->graph.nodes[t1[floor(t1.size() / 2.0)]];
				this->graph.nodes[t1[k + 1]].previousBifurcation = (t1.front());
				this->graph.nodes[t1[k + 1]].previousBifurcationPtr = &this->graph.nodes[t1.front()];
			}
			this->graph.nodes[t1[floor(t1.size() / 2.0)]].isMedian = true;

			for (int j = 0; j < paths.size(); j++) {
				t2 = paths[j].nodes;

				if ((t2.front() == t1.back()) && (t2.back() != t1.front())) {
					this->graph.nodes[t2.front()].isBifurcation = true;
					if (neighboursPerVertex[t2.back()].size() > 2) {
						for (int k = 0; k < t2.size() - 1; k++) {
							this->graph.nodes[t2[k]].nextBifurcation.push_back(t2.back());
							this->graph.nodes[t2[k]].nextBifurcationPtr.push_back(&this->graph.nodes[t2.back()]);
							this->graph.nodes[t2[k]].nextMedianPtr = &this->graph.nodes[t2[floor(t2.size() / 2.0)]];
							this->graph.nodes[t2[k + 1]].previousBifurcation = (t2.front());
							this->graph.nodes[t2[k + 1]].previousBifurcationPtr = &this->graph.nodes[t2.front()];
						}
						this->graph.nodes[t2[floor(t2.size() / 2) - 1]].isMedian = true;
						this->graph.nodes[t2.back()].isBifurcation = true;
					}
					if (neighboursPerVertex[t2.back()].size() == 1) {
						for (int k = 0; k < t2.size() - 1; k++) {
							if (this->parseTerminals) {
								this->graph.nodes[t2[k]].nextBifurcation.push_back(t2.back());
								this->graph.nodes[t2[k]].nextBifurcationPtr.push_back(&this->graph.nodes[t2.back()]);
							}

							this->graph.nodes[t2[k]].nextMedianPtr = &this->graph.nodes[t2[floor(t2.size() / 2.0)]];
							this->graph.nodes[t2[k + 1]].previousBifurcation = (t2.front());
							this->graph.nodes[t2[k + 1]].previousBifurcationPtr = &this->graph.nodes[t2.front()];
						}
						this->graph.nodes[t2[floor(t2.size() / 2) - 1]].isMedian = true;
						this->graph.nodes[t2.back()].isTerminal = true;
						this->graph.nodes[t2.back()].nextBifurcation.clear();
						this->graph.nodes[t2.back()].nextBifurcationPtr.clear();
					}

					t2.erase(t2.begin());
					r = t1;
					r.insert(r.end(), t2.begin(), t2.end());
					m.push_back(r);
				}
			}
		}
	}

	bool doContinue = true;
	while (doContinue) {
		doContinue = false;

		q = m;
		m.clear();
		for (int i = 0; i < q.size(); i++) {
			p = q[i];
			int lastindex = p[p.size() - 1];
			numOfNeighbours = neighboursPerVertex[lastindex].size();

			if (numOfNeighbours == 3) {
				doContinue = true;
				for (int j = 0; j < paths.size(); j++) {
					t2 = paths[j].nodes;
					//if ((t2.front() == lastindex) && (t2.back() != p.front())){
					if ((t2[0] == lastindex) && (t2[1] != p[p.size() - 2])) {
						this->graph.nodes[t2.front()].isBifurcation = true;
						if (neighboursPerVertex[t2.back()].size() == 3) {
							for (int k = 0; k < t2.size() - 1; k++) {
								this->graph.nodes[t2[k]].nextBifurcation.push_back(t2.back());
								this->graph.nodes[t2[k]].nextBifurcationPtr.push_back(&this->graph.nodes[t2.back()]);
								this->graph.nodes[t2[k]].nextMedianPtr = &this->graph.nodes[t2[floor(t2.size() / 2.0)]];
								this->graph.nodes[t2[k + 1]].previousBifurcation = (t2.front());
								this->graph.nodes[t2[k + 1]].previousBifurcationPtr = &this->graph.nodes[t2.front()];
							}
							this->graph.nodes[t2[floor(t2.size() / 2) - 1]].isMedian = true;
							this->graph.nodes[t2.back()].isBifurcation = true;
						}

						if (neighboursPerVertex[t2.back()].size() == 1) {
							for (int k = 0; k < t2.size() - 1; k++) {
								if (this->parseTerminals) {
									this->graph.nodes[t2[k]].nextBifurcation.push_back(t2.back());
									this->graph.nodes[t2[k]].nextBifurcationPtr.push_back(&this->graph.nodes[t2.back()]);
								}
								this->graph.nodes[t2[k]].nextMedianPtr = &this->graph.nodes[t2[floor(t2.size() / 2.0)]];
								this->graph.nodes[t2[k + 1]].previousBifurcation = (t2.front());
								this->graph.nodes[t2[k + 1]].previousBifurcationPtr = &this->graph.nodes[t2.front()];
							}
							this->graph.nodes[t2[floor(t2.size() / 2) - 1]].isMedian = true;
							this->graph.nodes[t2.back()].isTerminal = true;
							this->graph.nodes[t2.back()].nextBifurcation.clear();
							this->graph.nodes[t2.back()].nextBifurcationPtr.clear();
						}

						t2.erase(t2.begin());
						r = p;
						r.insert(r.end(), t2.begin(), t2.end());
						m.push_back(r);
					}
				}
			}
			if (numOfNeighbours == 1) {
				m.push_back(p);
			}
		}
	}

	return;
}

/*
void analysis::analyzeGraph2(void) {
	//This function aims to build a bifurcating graph
	//This function receives as exports the graph of the 1-D graph and
	//gets the last element of each path. Assuming that the beginning of each path is
	//a 3-n vertex. If the last element is a 1-n vertex it is taken as a starting point
	//For each 3-n point get the connecting paths.

	std::vector<std::vector<int>> m, q;
	std::vector<int> p, t1, t2, r;

	int numOfNeighbours;

	for (int i = 0; i < paths.size(); i++) {
		p = paths[i].nodes;
		int lastindex = p[p.size() - 1];
		numOfNeighbours = neighboursPerVertex[lastindex].size();
		if ((numOfNeighbours == 1) && (lastindex == inlet)) {
			t1 = p;
			reverse(t1.begin(), t1.end());

			this->graph.nodes[inlet].isInlet = true;
			this->inletPtr = &this->graph.nodes[inlet];
			this->graph.init = &this->graph.nodes[inlet];

			this->graph.nodes[t1.back()].isBifurcation = true;
			for (int k = 0; k < t1.size() - 1; k++) {
				this->graph.nodes[t1[k]].nextBifurcation.push_back(t1.back());
				this->graph.nodes[t1[k]].nextBifurcationPtr.push_back(&this->graph.nodes[t1.back()]);
				this->graph.nodes[t1[k]].nextMedianPtr = &this->graph.nodes[t1[floor(t1.size() / 2.0)]];
				this->graph.nodes[t1[k + 1]].previousBifurcation = (t1.front());
				this->graph.nodes[t1[k + 1]].previousBifurcationPtr = &this->graph.nodes[t1.front()];
			}
			this->graph.nodes[t1[floor(t1.size() / 2.0)]].isMedian = true;

			for (int j = 0; j < paths.size(); j++) {
				t2 = paths[j].nodes;

				if ((t2.front() == t1.back()) && (t2.back() != t1.front())) {
					this->graph.nodes[t2.front()].isBifurcation = true;
					if (neighboursPerVertex[t2.back()].size() > 2) {
						for (int k = 0; k < t2.size() - 1; k++) {
							this->graph.nodes[t2[k]].nextBifurcation.push_back(t2.back());
							this->graph.nodes[t2[k]].nextBifurcationPtr.push_back(&this->graph.nodes[t2.back()]);
							this->graph.nodes[t2[k]].nextMedianPtr = &this->graph.nodes[t2[floor(t2.size() / 2.0)]];
							this->graph.nodes[t2[k + 1]].previousBifurcation = (t2.front());
							this->graph.nodes[t2[k + 1]].previousBifurcationPtr = &this->graph.nodes[t2.front()];
						}
						this->graph.nodes[t2[floor(t2.size() / 2) - 1]].isMedian = true;
						this->graph.nodes[t2.back()].isBifurcation = true;
					}
					if (neighboursPerVertex[t2.back()].size() == 1) {
						for (int k = 0; k < t2.size() - 1; k++) {
							if (this->parseTerminals) {
								this->graph.nodes[t2[k]].nextBifurcation.push_back(t2.back());
								this->graph.nodes[t2[k]].nextBifurcationPtr.push_back(&this->graph.nodes[t2.back()]);
							}

							this->graph.nodes[t2[k]].nextMedianPtr = &this->graph.nodes[t2[floor(t2.size() / 2.0)]];
							this->graph.nodes[t2[k + 1]].previousBifurcation = (t2.front());
							this->graph.nodes[t2[k + 1]].previousBifurcationPtr = &this->graph.nodes[t2.front()];
						}
						this->graph.nodes[t2[floor(t2.size() / 2) - 1]].isMedian = true;
						this->graph.nodes[t2.back()].isTerminal = true;
						this->graph.nodes[t2.back()].nextBifurcation.clear();
						this->graph.nodes[t2.back()].nextBifurcationPtr.clear();
					}

					t2.erase(t2.begin());
					r = t1;
					r.insert(r.end(), t2.begin(), t2.end());
					m.push_back(r);
				}
			}
		}
	}

	bool doContinue = true;
	while (doContinue) {
		doContinue = false;

		q = m;
		m.clear();
		for (int i = 0; i < q.size(); i++) {
			p = q[i];
			int lastindex = p[p.size() - 1];
			numOfNeighbours = neighboursPerVertex[lastindex].size();

			if (numOfNeighbours == 3) {
				doContinue = true;
				for (int j = 0; j < paths.size(); j++) {
					t2 = paths[j].nodes;
					//if ((t2.front() == lastindex) && (t2.back() != p.front())){
					if ((t2[0] == lastindex) && (t2[1] != p[p.size() - 2])) {
						this->graph.nodes[t2.front()].isBifurcation = true;
						if (neighboursPerVertex[t2.back()].size() == 3) {
							for (int k = 0; k < t2.size() - 1; k++) {
								this->graph.nodes[t2[k]].nextBifurcation.push_back(t2.back());
								this->graph.nodes[t2[k]].nextBifurcationPtr.push_back(&this->graph.nodes[t2.back()]);
								this->graph.nodes[t2[k]].nextMedianPtr = &this->graph.nodes[t2[floor(t2.size() / 2.0)]];
								this->graph.nodes[t2[k + 1]].previousBifurcation = (t2.front());
								this->graph.nodes[t2[k + 1]].previousBifurcationPtr = &this->graph.nodes[t2.front()];
							}
							this->graph.nodes[t2[floor(t2.size() / 2) - 1]].isMedian = true;
							this->graph.nodes[t2.back()].isBifurcation = true;
						}

						if (neighboursPerVertex[t2.back()].size() == 1) {
							for (int k = 0; k < t2.size() - 1; k++) {
								if (this->parseTerminals) {
									this->graph.nodes[t2[k]].nextBifurcation.push_back(t2.back());
									this->graph.nodes[t2[k]].nextBifurcationPtr.push_back(&this->graph.nodes[t2.back()]);
								}
								this->graph.nodes[t2[k]].nextMedianPtr = &this->graph.nodes[t2[floor(t2.size() / 2.0)]];
								this->graph.nodes[t2[k + 1]].previousBifurcation = (t2.front());
								this->graph.nodes[t2[k + 1]].previousBifurcationPtr = &this->graph.nodes[t2.front()];
							}
							this->graph.nodes[t2[floor(t2.size() / 2) - 1]].isMedian = true;
							this->graph.nodes[t2.back()].isTerminal = true;
							this->graph.nodes[t2.back()].nextBifurcation.clear();
							this->graph.nodes[t2.back()].nextBifurcationPtr.clear();
						}

						t2.erase(t2.begin());
						r = p;
						r.insert(r.end(), t2.begin(), t2.end());
						m.push_back(r);
					}
				}
			}
			if (numOfNeighbours == 1) {
				m.push_back(p);
			}
		}
	}

	return;
}

*/

void analysis::exportGraphFeatures(std::string outfile) {
	gnode * g1 = new gnode();
	g1 = this->graph.init->nextBifurcationPtr[0];
	std::vector<gnode*> A;
	A.push_back(&this->graph.nodes[this->inlet]);
	std::ofstream mFile;
	mFile.open(outfile);
	while (!A.empty()) {
		gnode* q = A.back();
		std::vector<gnode*> B = q->nextBifurcationPtr;
		std::vector<gnode*> M = q->nextNodesPtr;
		gnode* P = q->previousBifurcationPtr;

		std::vector<gnode*> C;
		if (q->previousBifurcationPtr != nullptr) {
			C = P->nextBifurcationPtr;
		}
		//
		// && !C.empty()
		A.pop_back();
		Vector3f b1, b2, vpos = q->position;
		float theta = 0.0;
		if (B.size() == 2) {
			b1 = B[0]->position - vpos;
			b2 = B[1]->position - vpos;
			theta = acos(b1.dot(b2) / (b1.norm()*b2.norm())) * 360 / (2 * Pi);
		}
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty() && !M.empty()) {
				//Process
				float len = (B[i]->position - vpos).norm();

				std::string LR = M[i]->isLeft ? "R" : "L";
				mFile << "gen" << "," << M[i]->HorsfieldOrder << "," << M[i]->StrahlerOrderStage_1 << "," << M[i]->generation - 1 << "," << M[i]->diameter << "," << len << "," << theta << "," << LR << std::endl;
				//A.back()->position;
				/*if (B[i]->previousNodePtr->isLeft != B[i]->previousNodePtr->isRight) {
					B[i]->isLeft = B[i]->previousNodePtr->isLeft;
					B[i]->isRight = B[i]->previousNodePtr->isRight;
				}*/
			}
		}
		A.insert(A.end(), B.begin(), B.end());
	}
	mFile.close();
	return;
}

void analysis::generateGraphFeatures(void) {
	gnode * g1 = new gnode();

	g1 = this->graph.init->nextBifurcationPtr[0];

	std::vector < gnode *> g2 = g1->nextNodesPtr;

	if ((g2[0]->nextBifurcationPtr[0]->position - g1->position).norm() > (g2[1]->nextBifurcationPtr[0]->position - g1->position).norm()) {
		g2[0]->isRight = true;
		g2[1]->isLeft = true;
	}
	else {
		g2[0]->isLeft = true;
		g2[1]->isRight = true;
	}

	std::vector<gnode*> A;
	A.push_back(&this->graph.nodes[this->inlet]);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextNodesPtr;
		A.pop_back();
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				//Process
				if (B[i]->previousNodePtr->isLeft != B[i]->previousNodePtr->isRight) {
					B[i]->isLeft = B[i]->previousNodePtr->isLeft;
					B[i]->isRight = B[i]->previousNodePtr->isRight;
				}

				A.insert(A.end(), B.begin(), B.end());
			}
		}
	}

	A.clear();
	A.push_back(&this->graph.nodes[this->inlet]);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextNodesPtr;

		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				//Process
				for (int j = 0; j < this->skel->lines.size(); j++) {
					if ((B[i]->index == this->skel->lines[j][0] - 1) || (B[i]->index == this->skel->lines[j][1] - 1)) {
						if ((B[i]->previousNodePtr->index == this->skel->lines[j][0] - 1) || (B[i]->previousNodePtr->index == this->skel->lines[j][1] - 1)) {
							B[i]->line = j;
							B[i]->diameter = LocalDiameterPerSkeletonEdge[j];
							B[i]->meshVerticesIndices = SkeletonEdge2MeshVertices[j];
							for (int k = 0; k < B[i]->meshVerticesIndices.size(); k++) {
								int v = B[i]->meshVerticesIndices[k];
								B[i]->meshVerticesPositions.push_back(Vector3f(this->model->vertices[v][0], this->model->vertices[v][1], this->model->vertices[v][2]));
							}
						}
					}
				}
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}

	//this->graph.nodes[this->inlet].nextBifurcationPtr[0]->generation = 0;

	A.clear();

	this->graph.nodes[this->inlet].generation = 0;

	this->graph.nodes[this->inlet].segmentIndex = 0; //<=

	A.push_back(&this->graph.nodes[this->inlet]);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextBifurcationPtr;

		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				//Process
				B[i]->generation = 1 + B[i]->previousBifurcationPtr->generation;

				B[i]->segmentIndex = 2 * (B[i]->previousBifurcationPtr->segmentIndex + powf(2, (B[i]->previousBifurcationPtr->generation - 1))) + i - powf(2, (B[i]->generation - 1)); //<=
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}

	A.clear();
	A.push_back(&this->graph.nodes[this->inlet]);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextNodesPtr;
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				//Process
				if (!B[i]->isBifurcation) {
					//B[i]->generation = B[i]->previousNodePtr->generation;

					//B[i]->segmentIndex = B[i]->previousNodePtr->segmentIndex; //<=

					B[i]->generation = A.back()->nextBifurcationPtr[i]->generation;

					B[i]->segmentIndex = A.back()->nextBifurcationPtr[i]->segmentIndex; //<=
				}
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}

	return;
}

void analysis::generateGraphFeaturesLRDiscrimination(void) {


	this->graph.generateGraphFeaturesLRDiscrimination();

	return;
}

void analysis::generateGraphFeaturesGeneration(void) {
	std::vector<gnode*> A;

	//this->graph.nodes[this->inlet].nextBifurcationPtr[0]->generation = 0;

	this->graph.nodes[this->inlet].generation = 0;
	this->graph.nodes[this->inlet].segmentIndex = 0;
	this->graph.nodes[this->inlet].HorsfieldOrder = -1;
	this->graph.nodes[this->inlet].StrahlerOrderStage_1 = -1;
	A.clear();
	A.push_back(&this->graph.nodes[this->inlet]);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextBifurcationPtr;
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				//Process
				B[i]->generation = 1 + B[i]->previousBifurcationPtr->generation;
				B[i]->segmentIndex = 2 * (B[i]->previousBifurcationPtr->segmentIndex + powf(2, (B[i]->previousBifurcationPtr->generation - 1))) + i - powf(2, (B[i]->generation - 1)); //<=
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}

	//Extracting Horsfield Order
	A.clear();
	A.push_back(&this->graph.nodes[this->inlet]);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextBifurcationPtr;
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				B[i]->HorsfieldOrder = -1;
				B[i]->StrahlerOrderStage_1 = -1;
				if (B[i]->isTerminal) {
					B[i]->HorsfieldOrder = 1;
					B[i]->StrahlerOrderStage_1 = 1;
				}
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}

	while (this->graph.nodes[this->inlet].HorsfieldOrder < 1) {
		A.clear();
		A.push_back(&this->graph.nodes[this->inlet]);
		while (!A.empty()) {
			std::vector<gnode*> B = A.back()->nextBifurcationPtr;
			if (A.back()->HorsfieldOrder < 1) {
				int minHorsfieldVal = 10000;
				int maxHorsfieldVal = -10000;

				int minStrahlerOrder = 10000;
				int maxStrahlerOrder = -10000;

				for (int i = 0; i < B.size(); i++) {
					if (!B.empty()) {
						if (B[i]->HorsfieldOrder < minHorsfieldVal) {
							minHorsfieldVal = B[i]->HorsfieldOrder;
						}
						if (B[i]->HorsfieldOrder > maxHorsfieldVal) {
							maxHorsfieldVal = B[i]->HorsfieldOrder;
						}

						if (B[i]->StrahlerOrderStage_1 < minStrahlerOrder) {
							minStrahlerOrder = B[i]->StrahlerOrderStage_1;
						}
						if (B[i]->StrahlerOrderStage_1 > maxStrahlerOrder) {
							maxStrahlerOrder = B[i]->StrahlerOrderStage_1;
						}
					}
				}
				if (minHorsfieldVal < 1) {
				}
				else {
					A.back()->HorsfieldOrder = maxHorsfieldVal + 1;
					A.back()->StrahlerOrderStage_1 = minStrahlerOrder + 1;
				}
			}
			A.pop_back();
			if (!B.empty()) {
				A.insert(A.end(), B.begin(), B.end());
			}
		}
	}

	//---------------------------------------------------------------
	A.clear();
	A.push_back(&this->graph.nodes[this->inlet]);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextNodesPtr;
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				//Process
				if (!B[i]->isBifurcation) {
					//B[i]->generation = B[i]->previousNodePtr->generation;

					//B[i]->segmentIndex = B[i]->previousNodePtr->segmentIndex; //<=

					B[i]->generation = A.back()->nextBifurcationPtr[i]->generation;
					B[i]->HorsfieldOrder = A.back()->nextBifurcationPtr[i]->HorsfieldOrder;
					B[i]->StrahlerOrderStage_1 = A.back()->nextBifurcationPtr[i]->StrahlerOrderStage_1;
					B[i]->segmentIndex = A.back()->nextBifurcationPtr[i]->segmentIndex;
				}
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}

	return;
}

//TODO, Incomplete
/*void analysis::generateGraphFeaturesDiameter(dotObj * geom) {
	for (int i = 0; i < geom->vertices.size(); i++) {
	}

	return;
}*/

void analysis::getSkeletonVertex2Generation(dotObj & skel) {
	bool found;
	int fullPathIndex;
	int fullPathIndexPos;
	this->SkeletonVertex2Generation.clear();
	this->SkeletonVertex2Generation.resize(skel.vertices.size());
	for (int i = 0; i < skel.vertices.size(); i++) {
		found = false;
		for (int j = 0; j < this->fullpathsSFInlet.size(); j++) {
			for (int k = 0; k < this->fullpathsSFInlet[j].size(); k++) {
				if (!found) {
					if (this->fullpathsSFInlet[j][k] == i) {
						found = true;
						fullPathIndex = j;
						fullPathIndexPos = k;
					}
				}
			}
		}
		int numOfThreeNeighb = 0;;
		for (int k = 0; k < fullPathIndexPos; k++) {
			if (this->neighboursPerVertex[fullpathsSFInlet[fullPathIndex][k]].size() > 2) {
				numOfThreeNeighb++;
			}
		}
		this->SkeletonVertex2Generation[i] = numOfThreeNeighb;
	}
	return;
}

distance analysis::getMaxPath(void) {
	distance d;
	d.index = 0;
	d.value = this->paths[0].nodes.size();
	for (int i = 0; i < this->paths.size(); i++) {
		if (this->paths[i].nodes.size() > d.value) {
			d.value = this->paths[i].nodes.size();
			d.index = i;
		}
	}
	return d;
}

//status constructor
status::status(void) {
	this->isComplete = false;
	this->error = false;
	progress = 0;
	action = "";
}

//Retrieve SDF values for a dotObj object
void dotObj::getSDF(double cone, size_t number_of_rays, bool postprocess, bool clustering, size_t number_of_clusters, double smoothing_lambda, bool exportClusters) {
	while (this->sdf_property_map.size() > 0) { this->sdf_property_map.pop_back(); }
	while (this->segment_property_map_per_vertex.size() > 0) { this->segment_property_map_per_vertex.pop_back(); }
	while (this->segment_property_map.size() > 0) { this->segment_property_map.pop_back(); }
	// create and read Polyhedron----------------------------------------------------------------------------------
	Polyhedron mesh;
	std::stringstream input = this->toOFF();
	//toOFF("input.off");
	//std::ifstream input("input.off");
	if (!input || !(input >> mesh) || mesh.empty()) {
		std::cerr << "Not a valid off file." << std::endl;
		return;
	}
	// assign id field for each facet--------------------------------------------------------------------------------
	std::size_t facet_id = 0;
	for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin();
		facet_it != mesh.facets_end(); ++facet_it, ++facet_id) {
		facet_it->id() = facet_id;
	}

	// create a property-map for SDF values-----------------------------------------------------------------------
	std::vector<double> sdf_values(mesh.size_of_facets());
	Facet_with_id_pmap<double> sdf_property_map(sdf_values);

	//Caution
	std::cout << "Get SDF values" << std::endl;
	CGAL::sdf_values(mesh, sdf_property_map, cone, number_of_rays, postprocess);//<============

	// access SDF values (with constant-complexity)
	std::cout << "Access SDF values" << std::endl;
	for (Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it) {
		this->sdf_property_map.push_back(sdf_property_map[facet_it]);
	}

	std::cout << "Some exports" << std::endl;
	std::string SDFExportDataCSVFile;
	SDFExportDataCSVFile = "sdf_property_map_per_face_starting_from_zero_for_cone_";
	SDFExportDataCSVFile += std::to_string(cone);
	SDFExportDataCSVFile += ".csv";
	this->functions.printVectorToFile(this->sdf_property_map, SDFExportDataCSVFile);

	if (clustering) {
		// create a property-map for segment-ids----------------------------------------------------------------------
		std::vector<std::size_t> segment_ids(mesh.size_of_facets());
		Facet_with_id_pmap<std::size_t> segment_property_map(segment_ids);
		//const std::size_t number_of_clusters = 4;       // use 4 clusters in soft clustering
		//const double smoothing_lambda = 0.0001;  // importance of surface features, suggested to be in-between [0,1]
		CGAL::segmentation_from_sdf_values(mesh, sdf_property_map, segment_property_map, number_of_clusters, smoothing_lambda, exportClusters); //<==========
		// access segment-ids (with constant-complexity)
		for (Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin();
			facet_it != mesh.facets_end(); ++facet_it) {
			this->segment_property_map.push_back(segment_property_map[facet_it]);
		}
		this->functions.printVectorToFile(this->segment_property_map, "segment_property_map_per_face_starting_from_zero.txt");

		//WritePerVertex
		for (int k = 0; k < this->vertices.size(); k++) {
			segment_property_map_per_vertex.push_back(0);
		}
		for (int k = 0; k < this->faces.size(); k++) {
			int v0, v1, v2, group;
			group = this->segment_property_map.at(k);
			v0 = this->faces.at(k).at(0) - 1;
			v1 = this->faces.at(k).at(3) - 1;
			v2 = this->faces.at(k).at(6) - 1;
			segment_property_map_per_vertex.at(v0) = group;
			segment_property_map_per_vertex.at(v1) = group;
			segment_property_map_per_vertex.at(v2) = group;
		}
	}

	return;
}

//Refine SDF values so that none is negative
void dotObj::refineSDF(void) {
	for (int i = 0; i < this->sdf_property_map.size(); i++) {
		if ((this->sdf_property_map.at(i) < 0.1) || (isnan(this->sdf_property_map.at(i)))) {
			this->sdf_property_map.at(i) = 0.1;
		}
	}
	return;
}

std::vector<float> dotObj::getSDFPropertyMap(double cone, int number_of_rays, bool postprocess)
{
	std::vector<float> sdf_property_map_temp;

	std::stringstream ifd = toOFF();
	//toOFF("input.off");

	// create and read Polyhedron----------------------------------------------------------------------------------
	Polyhedron mesh;

	ifd >> mesh;

	// assign id field for each facet--------------------------------------------------------------------------------
	std::size_t facet_id = 0;
	for (Polyhedron::Facet_iterator facet_it = mesh.facets_begin();
		facet_it != mesh.facets_end(); ++facet_it, ++facet_id) {
		facet_it->id() = facet_id;
	}

	// create a property-map for SDF values-----------------------------------------------------------------------
	std::vector<double> sdf_values(mesh.size_of_facets());
	Facet_with_id_pmap<double> sdf_property_map(sdf_values);

	//Caution
	CGAL::sdf_values(mesh, sdf_property_map, cone, number_of_rays, postprocess);//<============

	// access SDF values (with constant-complexity)
	for (Polyhedron::Facet_const_iterator facet_it = mesh.facets_begin(); facet_it != mesh.facets_end(); ++facet_it) {
		sdf_property_map_temp.push_back(sdf_property_map[facet_it]);
	}

	this->functions.printVectorToFile(sdf_property_map_temp, "sdf_property_map_per_face_starting_from_zero.txt");
	return sdf_property_map_temp;
}

void simulation::writeFile() {
	if (lungmodel.at(lungmodel.size() - 1).vertices.size() > 0) {
		if (lungmodel.size() > 2) {
			lungmodel.erase(lungmodel.begin());
		}
		lungmodel.push_back(*tempModel);
		lungmodel.at(lungmodel.size() - 1).exportToFile("update");
	}

	return;
}

void simulation::loadFile(std::string modelName) {
	dotObj buffer;
	buffer.initializeFromFile(modelName);
	lungmodel.at(lungmodel.size() - 1) = buffer;

	return;
}

void simulation::exportSegmented() {
	lungmodel.at(lungmodel.size() - 1).exportToFileSegmented("segmented.obj");

	return;
}

void simulation::exportVerticesToText(void) {
	lungmodel.at(lungmodel.size() - 1).exportSelectedIndicesToFile("exportedVertices.txt");
	return;
}

void simulation::exportFile(void) {
	lungmodel.at(lungmodel.size() - 1).exportToFile("model.obj");

	return;
}

void simulation::sdf(std::string segmentOrCluster, double coneAngleSDF, int numberOfRaysSDF, bool postProcessingSDF, int numberOfClustersSDF, double lamdaSDF) {
	if (lungmodel.at(lungmodel.size() - 1).vertices.size() > 0) {
		bool exportClusters = false;
		if (segmentOrCluster == "Per cluster") { exportClusters = true; }
		lungmodel.at(lungmodel.size() - 1).getSDF(coneAngleSDF, numberOfRaysSDF, postProcessingSDF, true, numberOfClustersSDF, lamdaSDF, exportClusters);  //getSDF(0.005,25,true,4,0.001);
	}
	return;
}

void simulation::narrow(double contractionPercentageBox, int contractionStrengthBox, bool useCustomFunction, double frequencyBox, bool useInterpolation, double interpolationPercentage, bool useExtendedSmoothing, int  numberOfRaysSDF, int numberOfClustersSDF, double lamdaSDF) {
	std::vector<float> diff, initialSDF, finalSDF;
	std::vector<int> partIndices;
	double initialMeanSDF = 1.0, finalMeanSDF = 1.0;
	bool doSDF = true;
	narrowingRatio = 0.0;

	//INITIALIZE============================================
	dotObj *inputModel, *airway, buffer;
	airway = new dotObj;
	inputModel = new dotObj;
	*inputModel = lungmodel.at(lungmodel.size() - 1);
	//Get selected vertices--------------------------------------------------------------------------------------
	inputModel->selectedVertices = lungmodel.at(lungmodel.size() - 1).selectedVertices;
	if (inputModel->selectedVertices.size() == 0) { return; }
	if (inputModel->selectedVertices.size() > 300000) { return; }

	//Segment---------------------------------------------------------------------------------------
	airway->registerPart(*inputModel);

	//SDF narrowing estimation=================================================
	if (doSDF) {
		airway->smoothingSparse("taubin", "weighted", 5, -0.7, 0.5);
		airway->getSDF(0.02, numberOfRaysSDF, false, false, numberOfClustersSDF, lamdaSDF, false);  //(4, 25, true, 4, 0.0001);
		airway->refineSDF();
		airway->functions.printVectorToFile(airway->sdf_property_map, "_before.txt");
		initialSDF = airway->sdf_property_map;
		initialMeanSDF = airway->functions.getSpecificMean(initialSDF);
		std::cout << "initialMeanSDF: " << initialMeanSDF << std::endl;
		this->localDIameterBeforeNarrowing = initialMeanSDF;
	}

	//=========================================================================
	//Get values from UI-------------------------------------------------------------------------------------------
	airway->customFunctionAmplitude = 1 - contractionPercentageBox;
	airway->customFunctionFrequency = frequencyBox;
	if (useCustomFunction) { airway->uniformContraction = false; }
	else { airway->uniformContraction = true; }
	airway->skeletonizationIterations = contractionStrengthBox;
	//Get weighting---------------------------------------------------------------------------------------------
	airway->geoDistPerVertexFromFurthestVertex = airway->geodesicBasedWeighting();
	//Apply laplacian narrowing ----------------------------------------------------------------------------------------

	if (!(useInterpolation)) {
		//airway->skeletonize("laplacian.horizontal", airway->skeletonizationIterations, 0);  //Weights!!
		airway->skeletonize(airway->skeletonizationIterations);
	}
	//Apply interpolation-----------------------------------------------------------------------------------------------
	if (useInterpolation) {
		buffer = *airway;
		//airway->skeletonize("laplacian.horizontal", 17, 0);
		airway->skeletonize(17);
		airway->interp(buffer, interpolationPercentage);  //FOR INTERPOLATION
	}
	//SDF narrowing estimation==================
	//===================================
	if (doSDF) {
		airway->smoothingSparse("taubin", "weighted", 5, -0.7, 0.5);
		std::cout << "0" << std::endl;
		airway->getSDF(0.02, numberOfRaysSDF, false, false, numberOfClustersSDF, lamdaSDF, false); //(4, 25, true, 4, 0.0001);
		std::cout << "1" << std::endl;
		airway->refineSDF();
		std::cout << "2" << std::endl;
		airway->functions.printVectorToFile(airway->sdf_property_map, "_after.txt");
		std::cout << "3" << std::endl;
		//Estimate SDF after narrowing
		finalSDF = airway->sdf_property_map;
		finalMeanSDF = airway->functions.getSpecificMean(finalSDF);
		this->localDIameterAfterNarrowing = finalMeanSDF;
		//
		if (finalMeanSDF < 0) { finalMeanSDF = 0; }
		std::cout << "finalMeanSDF: " << finalMeanSDF << std::endl;
		initialMeanSDF = abs(initialMeanSDF);
		finalMeanSDF = abs(finalMeanSDF);
		narrowingRatio = initialMeanSDF - finalMeanSDF;
		if (initialMeanSDF > 0) {
			narrowingRatio = narrowingRatio / initialMeanSDF;
		}
		else { narrowingRatio = 1.0; }
		if (narrowingRatio > 1.0) { narrowingRatio = 1.0; }
		std::cout << "narrowingRatio: " << narrowingRatio << std::endl;
	}

	//==================================
	//==================================

	//After processing reconnect segment--------------------------------------------------------------------------
	inputModel->getPositionsFromChild(*airway);
	delete airway;

	//Smooth the processed area
	if (useExtendedSmoothing) {
		inputModel->selectedVertices;
		for (int smoothingIters = 0; smoothingIters < 5; smoothingIters++) {
			inputModel->extendSelection();
		}
		inputModel->smoothingSparse("taubin", "weighted", 50, -0.7, 0.5, inputModel->selectedVertices);			//FOR EXTENDED SMOOTHING
	}

	/*if (false){
	dotObj tester;
	tester = *inputModel;
	inputModel->selectedVertices;
	for (int smoothingIters = 0; smoothingIters < 5; smoothingIters++){
	tester.extendSelection();
	}
	airway = new dotObj();
	airway->registerPart(tester);
	//airway->recalculateNormals();
	//airway->exportToFile("tester");
	airway->smoothingSparse("taubin", "weighted", 300, -0.7, 0.5);			//FOR EXTENDED SMOOTHING
	inputModel->getPositionsFromChild(*airway);
	//inputModel->smoothingSparse("taubin", "weighted", 10, -0.52, 0.5);
	//delete airway;
	}*/

	//Export to new obj
	inputModel->exportToFile("modelbuffer");	//object oriented implementation needed

	//Delete initial object
	//lungmodel.push_back(*inputModel);
	tempModel = new dotObj;
	*tempModel = *inputModel;
	delete inputModel;

	return;
}

void simulation::generateVolume(int numOfGenerations, int density, volume &thevol, dotObj boundary, std::string id) {
	int totalNumberOfPoints;
	int atLeastPerEdge;

	totalNumberOfPoints = (int)(powf(2.0, (float)numOfGenerations + 1));
	////SET BOUNDARY
	std::cout << std::endl << "Action :: Setting boundary" << std::endl;
	dotObj * theParentObj;
	theParentObj = new dotObj();
	*theParentObj = boundary;
	//theParentObj->exportToFile("volumeBoundary" + id, "obj"); //->Exports outsite the create function

	////SET VOLUME
	volume * vol;
	vol = new volume();

	if (thevol.vertices.size() == 0) {
		std::cout << std::endl << "Action :: Building volume from obj boundary" << std::endl;
		atLeastPerEdge = 2 * (int)cbrt((float)totalNumberOfPoints);
		std::cout << std::endl << "Per edge points : " << atLeastPerEdge << std::endl;
		vol->createVolumeFromObj(*theParentObj, atLeastPerEdge);
		thevol = *vol;
		delete theParentObj;
	}
	else {
		std::cout << std::endl << "Action :: Getting volume from file" << std::endl;
		*vol = thevol;
	}
	delete vol;
	return;
}

void simulation::generate1DTree(int numOfGenerations, volume &thevol, dotObj &boundary, dotObj &hostMesh, dotObj &oneDim)
{
	int totalNumberOfPoints = (int)(powf(2.0, (float)numOfGenerations + 1));
	std::cout << std::endl << "Action :: Setting boundary" << std::endl;
	dotObj * theParentObj;
	theParentObj = new dotObj();
	*theParentObj = boundary;
	std::cout << std::endl << "Action :: Getting volume from file" << std::endl;
	volume * vol;
	vol = new volume();
	*vol = thevol;
	////GENERATE 1D REPRESENTATION
	treeGeneration *treeOne;
	treeOne = new treeGeneration();
	treeGeneration *treeTwo;
	treeTwo = new treeGeneration();
	std::cout << std::endl << "Action :: Building 1-dimensional representation" << std::endl;
	treeOne->volumeFillingSupplyDemand(*vol, hostMesh, totalNumberOfPoints); //->Exports inside the create function
	//treeOne->volumeFillingSupplyDemandWithHull(*vol,boundary, hostMesh, totalNumberOfPoints); //->Exports inside the create function
	treeOne->assignDiameters();

	oneDim = treeOne->oneDimension;
	treeTwo->oneDimension = treeOne->oneDimension;

	delete treeOne;
	delete vol;

	return;
}

void simulation::generate3DPointCloudFromTree(int numOfGenerations, int density, float radiusMultiplier, dotObj oneDRep, dotObj hostMesh, dotObj &treePCexp, std::string id, bool useHostMeshForPointCloud)
{
	treeGeneration *treeTwo;
	treeTwo = new treeGeneration();
	treeTwo->oneDimension = oneDRep;
	std::cout << std::endl << "Action :: Building point cloud" << std::endl;
	treeTwo->assignDiameters();

	treeTwo->diameterPerLineSegment[0] = 2 * hostMesh.mAnalysis.graph.init->diameter;
	//treeTwo->diameterPerLineSegment[0] = 2 * treeTwo->oneDimension.mAnalysis.graph.init->diameter;
	//treeTwo->diameterPerLineSegment[0] = 2 * 14.0;

	treeTwo->get3DPointCloudv3(useHostMeshForPointCloud, hostMesh, numOfGenerations, density, radiusMultiplier); //->Exports inside the create function
	std::cout << std::endl << "Action :: Clean up point cloud" << std::endl;
	//treeTwo->cleanUpv2();
	//treeTwo->cleanUpv3();

	treeTwo->cleanUpv3(radiusMultiplier);

	/*
	for (int i = 0; i < treeTwo->treePC.vertices.size();i++) {
		treeTwo->treePC.vertices[i][0]= treeTwo->treePC.vertices[i][0]+this->functions.normalRandom();
		treeTwo->treePC.vertices[i][1] = treeTwo->treePC.vertices[i][1]*this->functions.normalRandom();
		treeTwo->treePC.vertices[i][2] = treeTwo->treePC.vertices[i][2]*this->functions.normalRandom();
	}
	*/

	treePCexp = treeTwo->treePC;

	//treeTwo->treePC.exportToFile("treePointCloud" + id, "obj");
	delete treeTwo;
	return;
}

void simulation::generate3DPointCloudFromTree(int density, float radiusMultiplier, dotObj oneDRep, dotObj &treePCexp, std::vector<float> diameters, std::string id) {
	treeGeneration *treeTwo;
	treeTwo = new treeGeneration();
	treeTwo->oneDimension = oneDRep;
	std::cout << std::endl << "Action :: Building point cloud" << std::endl;
	treeTwo->assignDiameters(diameters);
	treeTwo->get3DPointCloudv3(density, radiusMultiplier); //->Exports inside the create function
	std::cout << std::endl << "Action :: Clean up point cloud" << std::endl;
	treeTwo->cleanUpv3();
	treePCexp = treeTwo->treePC;
	delete treeTwo;
	return;
}

void simulation::generate3DPointCloudFromGraph(
	int numOfGenerations,
	int density,
	float radiusMultiplier,
	ggraph * theGraph,
	dotObj * solid) {
	std::vector<gnode*> A;
	bool doClean = true;

	//Stoping criteria for sampling 
	A.clear();
	A.push_back(theGraph->init);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextNodesPtr;
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				if (B[i]->generation > numOfGenerations) {
					B[i]->doStop = true;
				}
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}




	//Start sampling 
	A.clear();
	A.push_back(theGraph->init->nextBifurcationPtr[0]->nextNodesPtr[0]);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextNodesPtr;
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				if (B[i]->meshVerticesPositions.empty()) {
					dotObj * geom = new dotObj();
					treeGeneration*treeTwo = new treeGeneration();
					treeTwo->generateCloudForLineSegment(
						geom,
						A.back()->position,
						B[i]->position,
						density,
						radiusMultiplier,
						0.27,
						0,
						true,
						radiusMultiplier*A.back()->diameter,
						radiusMultiplier*B[i]->diameter
					);
					solid->vertices.insert(std::end(solid->vertices), std::begin(geom->vertices), std::end(geom->vertices));
					delete geom;
					delete treeTwo;

					if (B[i]->isTerminal || B[i]->doStop) {
						dotObj * geom = new dotObj();
						treeGeneration*treeTwo = new treeGeneration();
						treeTwo->generateCloudForLineTip(
							geom,
							A.back()->position,
							B[i]->position,
							density,
							radiusMultiplier,
							0.27,
							0,
							true,
							B[i]->diameter
						);
						solid->vertices.insert(std::end(solid->vertices), std::begin(geom->vertices), std::end(geom->vertices));
						delete geom;
						delete treeTwo;
					}
				}
			}
		}
		A.pop_back();
		if (!B.empty()) {
			//A.insert(A.end(), B.begin(), B.end());
			for (int i = 0; i < B.size(); i++) {
				if (!B[i]->doStop) {
					A.push_back(B[i]);
				}
			}
		}
	}

	A.clear();
	A.push_back(theGraph->init->nextBifurcationPtr[0]->nextNodesPtr[1]);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextNodesPtr;
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				if (B[i]->meshVerticesPositions.empty()) {
					dotObj * geom = new dotObj();
					treeGeneration*treeTwo = new treeGeneration();
					treeTwo->generateCloudForLineSegment(
						geom,
						A.back()->position,
						B[i]->position,
						density,
						radiusMultiplier,
						0.27,
						0,
						true,
						radiusMultiplier*A.back()->diameter,
						radiusMultiplier*B[i]->diameter
					);
					solid->vertices.insert(std::end(solid->vertices), std::begin(geom->vertices), std::end(geom->vertices));
					delete geom;
					delete treeTwo;

					if (B[i]->isTerminal || B[i]->doStop) {
						dotObj * geom = new dotObj();
						treeGeneration*treeTwo = new treeGeneration();
						treeTwo->generateCloudForLineTip(
							geom,
							A.back()->position,
							B[i]->position,
							density,
							radiusMultiplier,
							0.27,
							0,
							true,
							B[i]->diameter
						);
						solid->vertices.insert(std::end(solid->vertices), std::begin(geom->vertices), std::end(geom->vertices));
						delete geom;
						delete treeTwo;
					}
				}
			}
		}
		A.pop_back();
		if (!B.empty()) {
			//A.insert(A.end(), B.begin(), B.end());
			for (int i = 0; i < B.size(); i++) {
				if (!B[i]->doStop) {
					A.push_back(B[i]);
				}
			}
		}
	}

	A.clear();
	A.push_back(theGraph->init);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextNodesPtr;
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				dotObj * geom = new dotObj();
				geom->functions.eigentostdvec(B[i]->meshVerticesPositions, geom->vertices);
				solid->vertices.insert(std::end(solid->vertices), std::begin(geom->vertices), std::end(geom->vertices));
				delete geom;
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}



	//Clean up
	if (doClean) {
		dotObj * solidClear = new dotObj;
		for (int j = 0; j < solid->vertices.size(); j++) {
			//int lastval = 0;
			//int progress = (int)(j / (solid->vertices.size() / 100));
			//if (progress!=lastval) {
			//	lastval = progress;
			//	std::cout << progress << "%" << std::endl;
			//}
			Vector3f pointInQuestion = solid->functions.stdtoeigenvec(solid->vertices[j]);
			bool doKeep = true;
			A.clear();
			A.push_back(theGraph->init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr);
			while (!A.empty()) {
				std::vector<gnode*> B = A.back()->nextNodesPtr;
				for (int i = 0; i < B.size(); i++) {
					if (!B.empty()) {
						Vector3f a = A.back()->position;
						Vector3f b = B[i]->position;
						Vector3f x = solid->functions.projPointOnVector(pointInQuestion, a, b);
						float lim0 = radiusMultiplier * A.back()->diameter;
						float lim1 = radiusMultiplier * B[i]->diameter;
						float limCutOff = 1.05;
						float dist = ((a - x).norm() / (a - b).norm());
						float pointToProj = (x - pointInQuestion).norm();
						if ((a - x).norm() + (x - b).norm() <= limCutOff * (a - b).norm()) {
							if ((pointToProj < lim0 + dist * (lim1 - lim0))) {
								doKeep = false;
							}
						}
						if ((B[i]->isTerminal || B[i]->doStop) && (0.1*(a - b).norm() > (x - b).norm())) {
							doKeep = true;
						}
					}
				}
				A.pop_back();
				if (!B.empty()) {
					//A.insert(A.end(), B.begin(), B.end());
					for (int i = 0; i < B.size(); i++) {
						if (!B[i]->doStop) {
							A.push_back(B[i]);
						}
					}
				}
			}

			A.clear();
			A.push_back(theGraph->init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr);
			while (!A.empty()) {
				std::vector<gnode*> B = A.back()->nextNodesPtr;
				for (int i = 0; i < B.size(); i++) {
					if (!B.empty()) {
						Vector3f a = A.back()->position;
						Vector3f b = B[i]->position;
						Vector3f x = solid->functions.projPointOnVector(pointInQuestion, a, b);
						float lim0 = radiusMultiplier * A.back()->diameter;
						float lim1 = radiusMultiplier * B[i]->diameter;
						float limCutOff = 1.05;
						float dist = ((a - x).norm() / (a - b).norm());
						float pointToProj = (x - pointInQuestion).norm();
						if ((a - x).norm() + (x - b).norm() <= limCutOff * (a - b).norm()) {
							if ((pointToProj < lim0 + dist * (lim1 - lim0))) {
								doKeep = false;
							}
						}
						if ((B[i]->isTerminal || B[i]->doStop) && (0.1*(a - b).norm() > (x - b).norm())) {
							doKeep = true;
						}
					}
				}
				A.pop_back();
				if (!B.empty()) {
					//A.insert(A.end(), B.begin(), B.end());
					for (int i = 0; i < B.size(); i++) {
						if (!B[i]->doStop) {
							A.push_back(B[i]);
						}
					}
				}
			}

			if (doKeep) {
				solidClear->vertices.push_back(solid->vertices[j]);
			}
		}

		solid->vertices = solidClear->vertices;
		delete solidClear;
	}


	//Remove stopping criteria
	A.clear();
	A.push_back(theGraph->init);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextNodesPtr;
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				B[i]->doStop = false;
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}



	return;
}

void simulation::segmentTree2(dotObj &oneDL, dotObj &oneDR, dotObj &recModel, std::string result, std::string keyword) {
	treeGeneration *treeTwo;
	treeTwo = new treeGeneration();

	treeTwo->reconstructedModel = recModel;

	treeTwo->reconstructedModel.segment_property_map_per_vertex.resize(treeTwo->reconstructedModel.vertices.size());
	for (int i = 0; i < treeTwo->reconstructedModel.segment_property_map_per_vertex.size(); i++) {
		treeTwo->reconstructedModel.segment_property_map_per_vertex[i] = -5;
	}
	treeTwo->reconstructedModel.segment_property_map.resize(treeTwo->reconstructedModel.faces.size());
	for (int i = 0; i < treeTwo->reconstructedModel.faces.size(); i++) {
		treeTwo->reconstructedModel.segment_property_map[i] = -5;
	}
	treeTwo->outletCounter = 0;

	//treeTwo->oneDimension = oneDL;
	//treeTwo->flatten(result, "L");

	//treeTwo->oneDimension = oneDR;
	//treeTwo->flatten(result, "R");

	//treeTwo->reconstructedModel.smoothingSparse("taubin", "weighted", 2, -0.53, 0.60);

	treeTwo->oneDimension = oneDL;
	treeTwo->naming2(result, "L");

	treeTwo->oneDimension = oneDR;
	treeTwo->naming2(result, "R");

	for (int i = 0; i < treeTwo->reconstructedModel.faces.size(); i++) {
		int vm1 = treeTwo->reconstructedModel.faces.at(i).at(0) - 1;
		int vm2 = treeTwo->reconstructedModel.faces.at(i).at(3) - 1;
		int vm3 = treeTwo->reconstructedModel.faces.at(i).at(6) - 1;

		int group1 = treeTwo->reconstructedModel.segment_property_map_per_vertex.at(vm1);
		int group2 = treeTwo->reconstructedModel.segment_property_map_per_vertex.at(vm2);
		int group3 = treeTwo->reconstructedModel.segment_property_map_per_vertex.at(vm3);
		if ((group1 == group2) && (group3 == group2)) { treeTwo->reconstructedModel.segment_property_map.at(i) = group1; }
		else {
			if ((group1 != 99) && (group3 != 99) && (group2 != 99)) {
				if (group1 == group2) { treeTwo->reconstructedModel.segment_property_map.at(i) = group1; }

				if (group2 == group3) { treeTwo->reconstructedModel.segment_property_map.at(i) = group2; }

				if (group3 == group1) { treeTwo->reconstructedModel.segment_property_map.at(i) = group3; }

				if ((group1 != group2) && (group2 != group3)) { treeTwo->reconstructedModel.segment_property_map.at(i) = group1; }
			}
		}
	}

	/*
	std::cout << "Categorizing non segmented" << std::endl;
	int uncateg = 0;
	for (int i = 0; i < treeTwo->reconstructedModel.segment_property_map.size(); i++){
	if (treeTwo->reconstructedModel.segment_property_map[i] ==-5){ uncateg++; }
	}
	for (int k = 0; k < uncateg; k++){
	std::cout << "Iteration" << k << std::endl;
	for (int i = 0; i < treeTwo->reconstructedModel.faces.size(); i++){
	if (treeTwo->reconstructedModel.segment_property_map.at(i) == -5){
	std::cout << "Found one" << std::endl;

	int vi1 = treeTwo->reconstructedModel.faces.at(i).at(0);
	int vi2 = treeTwo->reconstructedModel.faces.at(i).at(3);
	int vi3 = treeTwo->reconstructedModel.faces.at(i).at(6);

	bool c1 = false, c2 = false, c3 = false;
	for (int j = 0; j < treeTwo->reconstructedModel.faces.size(); j++){
	if (treeTwo->reconstructedModel.segment_property_map.at(j) != -5){
	int vj1 = treeTwo->reconstructedModel.faces.at(j).at(0);
	int vj2 = treeTwo->reconstructedModel.faces.at(j).at(3);
	int vj3 = treeTwo->reconstructedModel.faces.at(j).at(6);
	//float fj1 = roundatdecimal(ref->vertices.at(vj1-1).at(0), presicion);
	//float fj2 = roundatdecimal(ref->vertices.at(vj2-1).at(1), presicion);
	//float fj3 = roundatdecimal(ref->vertices.at(vj3-1).at(2), presicion);
	if ((vi1 == vj1) || (vi1 == vj2) || (vi1 == vj3)){ c1 = true; }
	if ((vi2 == vj1) || (vi2 == vj2) || (vi2 == vj3)){ c2 = true; }
	if ((vi3 == vj1) || (vi3 == vj2) || (vi3 == vj3)){ c3 = true; }
	//if (k == 0){
	if (((!c1) && (c2) && (c3)) || ((c1) && (!c2) && (c3)) || ((c1) && (c2) && (!c3))){
	treeTwo->reconstructedModel.segment_property_map.at(i) = treeTwo->reconstructedModel.segment_property_map.at(j);
	}
	//}
	//if (c1 || c2 || c3){
	//	treeTwo->reconstructedModel.segment_property_map.at(i) = treeTwo->reconstructedModel.segment_property_map.at(j);
	//}
	}
	}
	}
	}
	}
	*/

	treeTwo->reconstructedModel.exportToFileSegmented(result, keyword);

	return;
}

void simulation::flatten(dotObj &oneDL, dotObj &oneDR, dotObj &recModel, std::string result, std::string keyword) {
	treeGeneration *treeTwo;
	treeTwo = new treeGeneration();

	treeTwo->reconstructedModel = recModel;

	treeTwo->reconstructedModel.segment_property_map_per_vertex.resize(treeTwo->reconstructedModel.vertices.size());
	for (int i = 0; i < treeTwo->reconstructedModel.segment_property_map_per_vertex.size(); i++) {
		treeTwo->reconstructedModel.segment_property_map_per_vertex[i] = -5;
	}
	treeTwo->reconstructedModel.segment_property_map.resize(treeTwo->reconstructedModel.faces.size());
	for (int i = 0; i < treeTwo->reconstructedModel.faces.size(); i++) {
		treeTwo->reconstructedModel.segment_property_map[i] = -5;
	}
	treeTwo->outletCounter = 0;

	//treeTwo->oneDimension = oneDL;
	//treeTwo->flatten(result, "L");

	//treeTwo->oneDimension = oneDR;
	//treeTwo->flatten(result, "R");

	//treeTwo->reconstructedModel.smoothingSparse("taubin", "weighted", 2, -0.53, 0.60);

	treeTwo->oneDimension = oneDL;
	treeTwo->flatten(result, "L");

	treeTwo->oneDimension = oneDR;
	treeTwo->flatten(result, "R");

	/*
	for (int i = 0; i < treeTwo->reconstructedModel.faces.size(); i++){
	int vm1 = treeTwo->reconstructedModel.faces.at(i).at(0) - 1;
	int vm2 = treeTwo->reconstructedModel.faces.at(i).at(3) - 1;
	int vm3 = treeTwo->reconstructedModel.faces.at(i).at(6) - 1;

	int group1 = treeTwo->reconstructedModel.segment_property_map_per_vertex.at(vm1);
	int group2 = treeTwo->reconstructedModel.segment_property_map_per_vertex.at(vm2);
	int group3 = treeTwo->reconstructedModel.segment_property_map_per_vertex.at(vm3);
	if ((group1 == group2) && (group3 == group2)){ treeTwo->reconstructedModel.segment_property_map.at(i) = group1; }
	else{
	if ((group1 != 99) && (group3 != 99) && (group2 != 99)){
	if (group1 == group2){ treeTwo->reconstructedModel.segment_property_map.at(i) = group1; }

	if (group2 == group3){ treeTwo->reconstructedModel.segment_property_map.at(i) = group2; }

	if (group3 == group1){ treeTwo->reconstructedModel.segment_property_map.at(i) = group3; }

	if ((group1 != group2) && (group2 != group3)){ treeTwo->reconstructedModel.segment_property_map.at(i) = group1; }
	}
	}
	}
	*/

	/*
	std::cout << "Categorizing non segmented" << std::endl;
	int uncateg = 0;
	for (int i = 0; i < treeTwo->reconstructedModel.segment_property_map.size(); i++){
	if (treeTwo->reconstructedModel.segment_property_map[i] ==-5){ uncateg++; }
	}
	for (int k = 0; k < uncateg; k++){
	std::cout << "Iteration" << k << std::endl;
	for (int i = 0; i < treeTwo->reconstructedModel.faces.size(); i++){
	if (treeTwo->reconstructedModel.segment_property_map.at(i) == -5){
	std::cout << "Found one" << std::endl;

	int vi1 = treeTwo->reconstructedModel.faces.at(i).at(0);
	int vi2 = treeTwo->reconstructedModel.faces.at(i).at(3);
	int vi3 = treeTwo->reconstructedModel.faces.at(i).at(6);

	bool c1 = false, c2 = false, c3 = false;
	for (int j = 0; j < treeTwo->reconstructedModel.faces.size(); j++){
	if (treeTwo->reconstructedModel.segment_property_map.at(j) != -5){
	int vj1 = treeTwo->reconstructedModel.faces.at(j).at(0);
	int vj2 = treeTwo->reconstructedModel.faces.at(j).at(3);
	int vj3 = treeTwo->reconstructedModel.faces.at(j).at(6);
	//float fj1 = roundatdecimal(ref->vertices.at(vj1-1).at(0), presicion);
	//float fj2 = roundatdecimal(ref->vertices.at(vj2-1).at(1), presicion);
	//float fj3 = roundatdecimal(ref->vertices.at(vj3-1).at(2), presicion);
	if ((vi1 == vj1) || (vi1 == vj2) || (vi1 == vj3)){ c1 = true; }
	if ((vi2 == vj1) || (vi2 == vj2) || (vi2 == vj3)){ c2 = true; }
	if ((vi3 == vj1) || (vi3 == vj2) || (vi3 == vj3)){ c3 = true; }
	//if (k == 0){
	if (((!c1) && (c2) && (c3)) || ((c1) && (!c2) && (c3)) || ((c1) && (c2) && (!c3))){
	treeTwo->reconstructedModel.segment_property_map.at(i) = treeTwo->reconstructedModel.segment_property_map.at(j);
	}
	//}
	//if (c1 || c2 || c3){
	//	treeTwo->reconstructedModel.segment_property_map.at(i) = treeTwo->reconstructedModel.segment_property_map.at(j);
	//}
	}
	}
	}
	}
	}
	*/

	treeTwo->reconstructedModel.exportToFile(result, keyword);

	return;
}

void simulation::extendBronchialTree(void) {
	bool buildVolumes = true;
	bool buildCenterline = true;
	bool build1DModel = true;
	bool surfaceSampling = true;
	bool buildNormals = true;
	bool build3DModel = true;
	bool refinements = false;
	bool segmentation = false;
	int depth = 8;
	int volumeDepth = 10;
	int density = 500;
	int poissonDepth = 11;										// poissonDepth->value()
	std::string path = "models\\";  //= workspacePath->text().toStdString() + "/";
	std::string boundaryMeshR = path + "L2.obj";					//= rightLungGeomPath->text().toStdString();
	std::string boundaryMeshL = path + "L1.obj";				//= leftLungGeomPath->text().toStdString();
	std::string existingModelMesh = path + "lungPartsFull.obj";	//= existingGeomPath->text().toStdString();
	dotObj * boundaryR = new dotObj();
	boundaryR->initializeFromFile(boundaryMeshR);
	dotObj * boundaryL = new dotObj();
	boundaryL->initializeFromFile(boundaryMeshL);
	std::cout << std::endl << "Model :: Existing geometry " << std::endl;
	dotObj *existingModel = new dotObj();
	existingModel->initializeFromFile(existingModelMesh);
	dotObj * trachea = new dotObj();
	trachea->initializeFromFile(path + "tracheav8PC.obj");
	if (strlen(path.c_str()) > 2) {
	}
	else {
		std::cout << "Error: Enter valid workspace path" << std::endl;
		return;
	}

	status * st = new status();

	this->extendBronchialTree(
		buildVolumes,
		buildCenterline,
		build1DModel,
		surfaceSampling,
		buildNormals,
		build3DModel,
		refinements,
		segmentation,
		depth,
		volumeDepth,
		density,
		poissonDepth,
		path,
		boundaryR,
		boundaryL,
		trachea,
		existingModel,
		st,
		false);
	return;
}

void simulation::extendBronchialTree(
	bool buildVolumes,
	bool buildCenterline,
	bool build1DModel,
	bool surfaceSampling,
	bool buildNormals,
	bool build3DModel,
	bool refinements,
	bool segmentation,
	int depth,
	int volumeDepth,
	int density,
	int poissonDepth,
	std::string path,
	dotObj * boundaryR,
	dotObj * boundaryL,
	dotObj * trachea,
	dotObj * existingModel,
	status * stat,
	bool verbose) {
	dotObj * hostR = new dotObj();
	dotObj * hostL = new dotObj();

	bool initvols = buildVolumes;
	bool initskeletonization = buildCenterline;
	bool init1D = build1DModel;
	bool initPC = surfaceSampling;
	bool initNormals = buildNormals;
	bool init3D = build3DModel;
	bool initRefine = refinements;
	bool initSegmentation = segmentation;
	bool preSmoothing = initskeletonization;
	bool buildSkel = initskeletonization;
	//1D
	bool doMatch = init1D;
	bool build1D = init1D;
	//PC
	bool buildPC = initPC;
	bool buildPLY = initPC;
	//Normals & Smoothing
	//Reconstruction
	bool build3D = init3D;
	//Refinements
	bool refine3D = initRefine;
	//Segmentation
	bool segmentModel = initSegmentation;
	//Depth
	int tree3DDepth = depth;
	int treeDepth = depth;
	if (existingModel->vertices.size() > 0 && existingModel->faces.size() && existingModel->normals.size() > 0 &&
		boundaryL->vertices.size() > 0 && boundaryL->faces.size() && boundaryL->normals.size() > 0 &&
		boundaryR->vertices.size() > 0 && boundaryR->faces.size() && boundaryR->normals.size() > 0) {
	}
	else {
		std::cout << "Error: Enter valid 3D models" << std::endl;

		return;
	}

	stat->action = "Simulation :: Commencing";

	std::cout << std::endl << "Simulation :: Commencing" << std::endl;
	simulation *sim;
	volume * volLeft;

	volume * volRight;
	dotObj * skelLeftGen;
	dotObj * skelLeftRight;
	dotObj cachedModel, skel;

	//BUILD VOLUMES
	if (buildVolumes) {
		stat->action = "Volume generation";
		skelLeftGen = new dotObj();
		skelLeftRight = new dotObj();
		volLeft = new volume();
		volRight = new volume();

		sim = new simulation();
		sim->generateVolume(volumeDepth, density, *volLeft, *boundaryL, "Left");
		delete sim;

		sim = new simulation();
		sim->generateVolume(volumeDepth, density, *volRight, *boundaryR, "Right");
		delete sim;
	}

	dotObj * skeletonPart1 = new dotObj();
	dotObj * skeletonPart2 = new dotObj();

	//BUILD SKELETON
	if (buildSkel) {
		stat->action = "Build Skeleton";

		cachedModel = *existingModel;
		//cachedModel.smoothingSparse("taubin", "onoff", 1000, -0.5, 0.45);
		skel = cachedModel.mcfskel();

		skel.mcfskelRefineStepOne(8);
		skel.mcfskelRefineStepTwo();

		skel.mcfskelRefine();
		skel.mcfskelRefine();

		skel.buildInitialTrees(*skeletonPart1, *skeletonPart2, 2);

		skeletonPart1->mcfskelRefine();
		skeletonPart2->mcfskelRefine();
	}

	//MATCH VOLUMES TO LINES======================================
	if (doMatch) {
		Vector3f M1 = skeletonPart1->findCenterOfMass();
		Vector3f M2 = skeletonPart2->findCenterOfMass();
		Vector3f B1 = boundaryL->findCenterOfMass();
		Vector3f B2 = boundaryR->findCenterOfMass();

		std::cout << "Cross checking" << std::endl;

		bool found1 = false;
		for (int i = 0; i < skeletonPart1->vertices.size(); i++) {
			std::cout << "skeletonPart1 not found inside boundaryL:" << found1 << std::endl;
			if (!found1) {
				if (!found1) {
					found1 = boundaryL->checkIfPointisInsidev2(Vector3f(skeletonPart1->vertices[i][0], skeletonPart1->vertices[i][1], skeletonPart1->vertices[i][2]));
				}
				else
				{
					std::cout << "skeletonPart1 found inside boundaryL:" << "found1 is true" << std::endl;
				}
			}
		}

		bool found2 = false;
		for (int i = 0; i < skeletonPart1->vertices.size(); i++) {
			std::cout << "skeletonPart1 not found inside boundaryR:" << found1 << std::endl;
			if (!found2) {
				if (!found2) {
					found2 = boundaryR->checkIfPointisInsidev2(Vector3f(skeletonPart1->vertices[i][0], skeletonPart1->vertices[i][1], skeletonPart1->vertices[i][2]));
				}
				else {
					std::cout << "skeletonPart1 found inside boundaryR:" << "found2 is true" << std::endl;
				}
			}
		}

		if (found1) {
			*hostL = *skeletonPart1;
			*hostR = *skeletonPart2;
		}
		else if (found2) {
			*hostL = *skeletonPart2;
			*hostR = *skeletonPart1;
		}
	}

	dotObj * skeletonLeftGen = new dotObj();
	dotObj * skeletonRightGen = new dotObj();

	//BUILD 1D
	if (build1D) {
		stat->action = "Extending geometry";
		std::cout << "Extending geometry" << std::endl;
		sim = new simulation();
		sim->generate1DTree(treeDepth, *volLeft, *boundaryL, *hostL, *skeletonLeftGen);
		if (verbose) skeletonLeftGen->exportToFile(path + "skeletonLeftGenerated");
		delete sim;
		sim = new simulation();
		sim->generate1DTree(treeDepth, *volRight, *boundaryR, *hostR, *skeletonRightGen);
		if (verbose) skeletonRightGen->exportToFile(path + "skeletonRightGenerated");
		delete sim;
	}

	dotObj * treePCExpL = new dotObj();
	dotObj * treePCExpR = new dotObj();

	//BUILD 3D
	if (buildPC) {
		stat->action = "Build 3D Point cloud";
		std::cout << "Build 3D Point cloud" << std::endl;
		sim = new simulation();
		sim->generate3DPointCloudFromTree(tree3DDepth, density, 0.5, *skeletonLeftGen, *hostL, *treePCExpL, "Left", false);
		delete sim;
		sim = new simulation();
		sim->generate3DPointCloudFromTree(tree3DDepth, density, 0.5, *skeletonRightGen, *hostR, *treePCExpR, "Right", false);
		delete sim;
	}

	//BUILD PLY FILE
	if (buildPLY) {
		treePCExpR->exportToXYZ(path + "lungPCRec", true);
		treePCExpL->exportToXYZ(path + "lungPCRec", false);
		dotObj * trachea = new dotObj();
		trachea->initializeFromFile(path + "tracheav8PC.obj");
		trachea->exportToXYZ(path + "lungPCRec", false);
		delete trachea;
	}

	//Calculation of normals and smoothing
	if (false) {
		std::cout << "Calculating normals" << std::endl;
		std::string strname = path + "lungPCRec.xyz";
		const char* fname = strname.c_str();
		// Reads a .xyz point set file in points[].
		std::vector<PointEx> pointsL;
		std::ifstream stream(fname);
		if (!stream ||
			!CGAL::read_xyz_points(stream, std::back_inserter(pointsL),
				CGAL::Identity_property_map<PointEx>()))
		{
			std::cerr << "Error: cannot read file " << fname << std::endl;
			return;
		}

		std::vector<PointEx> output;
		//parameters
		const double retain_percentage = 80.0;   // percentage of points to retain.
		const double neighbor_radius = 0.5;   // neighbors size.
		const bool require_uniform_sampling = true;
		const int max_iter_number = 100;
		CGAL::wlop_simplify_and_regularize_point_set
			<Concurrency_tag>
			(pointsL.begin(),
				pointsL.end(),
				std::back_inserter(output),
				retain_percentage,
				neighbor_radius,
				max_iter_number,
				require_uniform_sampling
				);

		std::ofstream out(path + "WLOP_lungPCRec.xyz");
		if (!out || !CGAL::write_xyz_points(
			out, output.begin(), output.end()))
		{
			//return EXIT_FAILURE;
		}
	}

	if (buildNormals) {
		stat->action = "Calculating normals";
		std::cout << "Calculating normals" << std::endl;
		std::string strname = path + "lungPCRec.xyz";
		const char* fname = strname.c_str();
		// Reads a .xyz point set file in points[].
		std::vector<PointEx> pointsL;
		std::ifstream stream(fname);
		if (!stream ||
			!CGAL::read_xyz_points(stream, std::back_inserter(pointsL),
				CGAL::Identity_property_map<PointEx>()))
		{
			std::cerr << "Error: cannot read file " << fname << std::endl;

			return;
		}

		/*
		// simplification by clustering using erase-remove idiom
		double cell_size = 0.001;
		pointsL.erase(CGAL::grid_simplify_point_set(pointsL.begin(), pointsL.end(), cell_size),
		pointsL.end());
		// Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
		std::vector<PointEx>(pointsL).swap(pointsL);*/

		/*CGAL::Timer task_timer; task_timer.start();
		// simplification by clustering using erase-remove idiom
		pointsL.erase(CGAL::hierarchy_simplify_point_set(pointsL.begin(), pointsL.end(),
		5, // Max cluster size
		0.0001), // Max surface variation
		pointsL.end());
		std::size_t memory = CGAL::Memory_sizer().virtual_size();

		std::cout << pointsL.size() << " point(s) kept, computed in "
		<< task_timer.time() << " seconds, "
		<< (memory >> 20) << " Mib allocated." << std::endl;
		*/

		std::list<PointVectorPair> points;

		for (int i = 0; i < pointsL.size(); i++) {
			PointVectorPair temp;
			temp.first = pointsL[i];
			points.push_back(temp);
		}

		/*std::ifstream stream(fname);
		if (!stream ||
		!CGAL::read_xyz_points(stream,
		std::back_inserter(points),
		CGAL::First_of_pair_property_map<PointVectorPair>()))
		{
		std::cerr << "Error: cannot read file " << fname << std::endl;
		return ;
		}*/
		stat->action = "Estimate normals direction";
		std::cout << "Estimate normals direction" << std::endl;
		// Note: pca_estimate_normals() requires an iterator over points
		// as well as property maps to access each point's position and normal.
		const int nb_neighbors = 40; // K-nearest neighbors = 3 rings
		CGAL::pca_estimate_normals<Concurrency_tag>(points.begin(), points.end(),
			CGAL::First_of_pair_property_map<PointVectorPair>(),
			CGAL::Second_of_pair_property_map<PointVectorPair>(),
			nb_neighbors);
		// Orients normals.
		// Note: mst_orient_normals() requires an iterator over points
		// as well as property maps to access each point's position and normal.
		std::list<PointVectorPair>::iterator unoriented_points_begin =
			CGAL::mst_orient_normals(points.begin(), points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>(),
				nb_neighbors);
		// Optional: delete points with an unoriented normal
		// if you plan to call a reconstruction algorithm that expects oriented normals.
		points.erase(unoriented_points_begin, points.end());

		std::ofstream out(path + "lungPC.ply");
		if (!out ||
			!CGAL::write_ply_points_and_normals(
				out, points.begin(), points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>()))
		{
			return;
		}
		stat->action = "Smoothing normals";
		std::cout << "Smoothing normals" << std::endl;
		// Algorithm parameters
		int k = 40;                 // size of neighborhood. The bigger the smoother the result will be.
		// This value should bigger than 1.
		double sharpness_angle = 25; // control sharpness of the result.
		// The bigger the smoother the result will be
		int iter_number = 4;         // number of times the projection is applied

		for (int i = 0; i < iter_number; ++i)
		{
			CGAL::bilateral_smooth_point_set <Concurrency_tag>(
				points.begin(),
				points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>(),
				k,
				sharpness_angle);
		}

		std::ofstream outSmoothedply(path + "lungPCSmoothed.ply");
		std::ofstream outSmoothedxyz(path + "lungPCSmoothed.xyz");
		if (!outSmoothedply ||
			!CGAL::write_ply_points_and_normals(
				outSmoothedply, points.begin(), points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>()))
		{
			return;
		}

		if (!outSmoothedxyz ||
			!CGAL::write_xyz_points_and_normals(
				outSmoothedxyz, points.begin(), points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>()))
		{
			return;
		}
	}

	if (build3D)
	{
		stat->action = "Poisson Recontruction Commencing";
		dotObj *mm;
		std::cout << "Poisson Recontruction Commencing" << std::endl;
		//PoissonRecon(path + "lungPCSmoothed.ply", path + "lungPC.screened.ply",poissonDepth->value());
		SSDRecon(path + "lungPCSmoothed.ply", path + "lungPC.screened.ply", poissonDepth, mm);
		//poissonRec(path + "lungPCSmoothed.xyz", path + "lungPCscreened.off");
		//isotropicRemesh(path + "lungPCscreened.off", path + "lungPCscreenedIso.off");
		std::cout << "Poisson Recontruction Complete" << std::endl;
	}

	/*

	//REFINE
	if (refine3D){
	sim = new simulation();

	std::cout << "Segmentation" << std::endl;
	//Reconstructed
	dotObj * Reconstructed = new dotObj();
	Reconstructed->initializeFromPLY(path + "lungPC.screened.ply");
	std::cout << "Recalculate normals" << std::endl;
	Reconstructed->recalculateNormals();
	Reconstructed->exportToFile(path + "lungReconstructed");
	Reconstructed->smoothingSparse("taubin", "weighted", 100, -0.7, 0.55);
	//Reconstructed->smoothingSparse("taubin", "weighted", 200, -0.7, 0.5);
	std::cout << "Recalculate normals" << std::endl;
	Reconstructed->recalculateNormals();
	Reconstructed->exportToFile(path + "lungReconstructedSmooth");

	//Left 1D
	std::string oneDR = path + "skeletonRightGenerated.obj";
	dotObj *oneDModelR = new dotObj();
	oneDModelR->initializeFromFile(oneDR);

	//Right 1D
	std::string oneDL = path + "skeletonLeftGenerated.obj";
	dotObj *oneDModelL = new dotObj();
	oneDModelL->initializeFromFile(oneDL);

	sim->flatten(*oneDModelL, *oneDModelR, *Reconstructed, path + "lungReconstructedFlattened", "obj");

	delete oneDModelR;
	delete oneDModelL;

	delete sim;
	}

	if (segmentModel){
	sim = new simulation();

	std::cout << "Segmentation" << std::endl;
	//Reconstructed
	dotObj * Reconstructed = new dotObj();
	//Reconstructed->initializeFromFile(path + "lungReconstructedFlattened.obj");
	Reconstructed->initializeFromFile(path + "lungReconstructedSmooth.obj");
	//Left 1D
	std::string oneDR = path + "skeletonRightGenerated.obj";
	dotObj *oneDModelR = new dotObj();
	oneDModelR->initializeFromFile(oneDR);

	//Right 1D
	std::string oneDL = path + "skeletonLeftGenerated.obj";
	dotObj *oneDModelL = new dotObj();
	oneDModelL->initializeFromFile(oneDL);

	sim->segmentTree2(*oneDModelL, *oneDModelR, *Reconstructed, path + "lungReconstructedSegmented", "full");
	delete oneDModelR;
	delete oneDModelL;
	//delete recModel;
	delete sim;
	}*/

	return;
}

void simulation::extendBronchialTreeV2(
	bool buildVolumes,
	bool buildCenterline,
	bool build1DModel,
	bool surfaceSampling,
	bool buildNormals,
	bool build3DModel,
	bool refinements,
	bool segmentation,
	int depth,
	int volumeDepth,
	int density,
	int poissonDepth,
	std::string path,
	dotObj * boundaryR,
	dotObj * boundaryL,
	dotObj * trachea,
	dotObj * existingModel,
	status * stat,
	bool verbose) {
	bool initvols = buildVolumes;
	bool initskeletonization = buildCenterline;
	bool init1D = build1DModel;
	bool initPC = surfaceSampling;
	bool initNormals = buildNormals;
	bool init3D = build3DModel;
	bool initRefine = refinements;
	bool initSegmentation = segmentation;
	bool preSmoothing = initskeletonization;
	bool buildSkel = initskeletonization;
	//1D
	bool doMatch = init1D;
	bool build1D = init1D;
	//PC
	bool buildPC = initPC;
	bool buildPLY = initPC;
	//Normals & Smoothing
	//Reconstruction
	bool build3D = init3D;
	//Refinements
	bool refine3D = initRefine;
	//Segmentation
	bool segmentModel = initSegmentation;
	//Depth
	int tree3DDepth = depth;
	int treeDepth = depth;
	if (existingModel->vertices.size() > 0 && existingModel->faces.size() && existingModel->normals.size() > 0 &&
		boundaryL->vertices.size() > 0 && boundaryL->faces.size() && boundaryL->normals.size() > 0 &&
		boundaryR->vertices.size() > 0 && boundaryR->faces.size() && boundaryR->normals.size() > 0) {
	}
	else {
		std::cout << "Error: Enter valid 3D models" << std::endl;
		return;
	}
	stat->action = "Simulation :: Commencing";
	std::cout << std::endl << "Simulation :: Commencing" << std::endl;
	simulation *sim;
	volume * volLeft;
	volume * volRight;
	dotObj * skelLeftGen;
	dotObj * skelLeftRight;
	dotObj cachedModel, skel;
	dotObj * hostR = new dotObj();
	dotObj * hostL = new dotObj();
	dotObj * recon = new dotObj();
	dotObj * fullTree = new dotObj();
	trachea = new dotObj;
	dotObj * skeletonLeftGen = new dotObj();
	dotObj * skeletonRightGen = new dotObj();

	//BUILD VOLUMES
	if (buildVolumes) {
		stat->action = "Volume generation";
		skelLeftGen = new dotObj();
		skelLeftRight = new dotObj();
		volLeft = new volume();
		volRight = new volume();

		sim = new simulation();
		sim->generateVolume(volumeDepth, density, *volLeft, *boundaryL, "Left");
		delete sim;

		sim = new simulation();
		sim->generateVolume(volumeDepth, density, *volRight, *boundaryR, "Right");
		delete sim;
	}

	//BUILD SKELETON
	if (buildSkel) {
		stat->action = "Build Skeleton";
		cachedModel = *existingModel;
		//cachedModel.smoothingSparse("taubin", "onoff", 2000, -0.95, 0.35);
		//cachedModel.exportToFile(path + "smoothed");

		skel = cachedModel.mcfskel();

		skel.exportToFile(path + "_1_skeletonRaw");

		skel.mcfskelRefineStepOne(8);
		skel.mcfskelRefineStepOne(8);
		skel.mcfskelRefineStepOne(8);
		skel.mcfskelRefineStepTwo();
		skel.exportToFile(path + "_2_skeletonRefined");

		//skel.mcfskelRefine();
		//skel.mcfskelRefine();
		/*skel.graphBasedPruning();
		ggraph skk;
		skk = skel.mAnalysis.graph;
		dotObj * sk = new dotObj();
		skk.ggraph2Model(sk);
		sk->exportToFile(path + "skOut2");
		skel = *sk;
		skel.exportToFile(path + "skeletonOut2");*/
		//=============================================================
		/*
		dotObj * treePCExp;
		simulation * sim;
		treePCExp = new dotObj();
		sim = new simulation();
		sim->generate3DPointCloudFromTree(120, 0.7, skel, *treePCExp, this->mAnalysis.LocalDiameterPerSkeletonEdge, "ok");
		treePCExp->exportToFile("../../../data/tree-0.5");
		delete sim;*/
		//=============================================================

		ggraph rebuilt;
		existingModel->mAnalysis.parseTerminals = true;
		existingModel->graphBasedAnalysis(skel);

		rebuilt = existingModel->mAnalysis.graph;
		recon = new dotObj();
		rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0]->doStop = true;
		rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0]->doStop = true;
		rebuilt.ggraph2Model(recon);
		recon->exportToFile(path + "_3_trachea_centerline");
		//===========================================================
		rebuilt = existingModel->mAnalysis.graph;
		trachea = new dotObj();
		rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0]->doStop = true;
		rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0]->doStop = true;
		rebuilt.ggraph2ModelV(trachea);
		trachea->exportToFile(path + "_4_trachea");
		//=============================================================
		/*std::vector<gnode*> A;
		A.clear();
		A.push_back(rebuilt->init);
		while (!A.empty()){
		std::vector<gnode*> B = A.back()->nextNodesPtr;
		A.pop_back();
		for (int i = 0; i < B.size(); i++){
		if (!B.empty()){
		if (B[i]->generation > 3){ B[i]->doStop = true; }
		}
		}
		}*/
		gnode * g;
		//=============================================================
		existingModel->graphBasedAnalysis(skel);
		for (int i = 0; i < existingModel->mAnalysis.graph.nodes.size(); i++) {
			if (existingModel->mAnalysis.graph.nodes[i].generation > 6)existingModel->mAnalysis.graph.nodes[i].doStop = true;
		}
		rebuilt = existingModel->mAnalysis.graph;
		g = new gnode();
		hostL = new dotObj();
		g = rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0];
		if (g->isLeft) {
			g = rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr;
		}
		else {
			g = rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr;
		}
		rebuilt.init = g;
		rebuilt.ggraph2Model(hostL, true);
		hostL->mAnalysis.generateGraph(*hostL, 0);
		hostL->mAnalysis.graph.init->diameter = g->diameter;
		hostL->exportToFile(path + "_5_host_centerline_left");
		//=============================================================
		existingModel->graphBasedAnalysis(skel);
		for (int i = 0; i < existingModel->mAnalysis.graph.nodes.size(); i++) {
			if (existingModel->mAnalysis.graph.nodes[i].generation > 3)existingModel->mAnalysis.graph.nodes[i].doStop = true;
		}
		rebuilt = existingModel->mAnalysis.graph;
		g = new gnode();
		hostR = new dotObj();
		g = rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0];
		if (g->isRight) {
			g = rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr;
		}
		else {
			g = rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr;
		}

		rebuilt.init = g;
		rebuilt.ggraph2Model(hostR, true);
		hostR->mAnalysis.generateGraph(*hostR, 0);
		hostR->mAnalysis.graph.init->diameter = g->diameter;
		hostR->exportToFile(path + "_5_host_centerline_right");
	}
	//BUILD 1D

	dotObj * fullTreeExtendedModel = new dotObj();
	ggraph fullTreeExtendedGraph;

	if (build1D) {
		stat->action = "Extending geometry";
		std::cout << "Extending geometry" << std::endl;
		sim = new simulation();
		sim->generate1DTree(treeDepth, *volLeft, *boundaryL, *hostL, *skeletonLeftGen);
		if (verbose) skeletonLeftGen->exportToFile(path + "_6_skeletonLeftGenerated");
		delete sim;
		sim = new simulation();
		sim->generate1DTree(treeDepth, *volRight, *boundaryR, *hostR, *skeletonRightGen);
		if (verbose) skeletonRightGen->exportToFile(path + "_6_skeletonRightGenerated");
		delete sim;
		skeletonLeftGen->mAnalysis.generateGraph(*skeletonLeftGen, 0);
		skeletonLeftGen->mAnalysis.graph.init->diameter = hostL->mAnalysis.graph.init->diameter;

		skeletonRightGen->mAnalysis.generateGraph(*skeletonRightGen, 0);
		skeletonRightGen->mAnalysis.graph.init->diameter = hostR->mAnalysis.graph.init->diameter;

		//Full tree rebuilt
		ggraph fullTreeGraph;
		existingModel->mAnalysis.parseTerminals = true;
		existingModel->graphBasedAnalysis(skel);

		fullTreeGraph = existingModel->mAnalysis.graph;
		dotObj * fullTreeModel = new dotObj();
		if (fullTreeGraph.init->nextBifurcationPtr[0]->nextNodesPtr[0]->isLeft) {
			fullTreeGraph.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0] = skeletonLeftGen->mAnalysis.graph.init;
			fullTreeGraph.init->nextBifurcationPtr[0]->nextBifurcationPtr[0] = skeletonLeftGen->mAnalysis.graph.init->nextBifurcationPtr[0];

			fullTreeGraph.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0] = skeletonRightGen->mAnalysis.graph.init;
			fullTreeGraph.init->nextBifurcationPtr[0]->nextBifurcationPtr[1] = skeletonRightGen->mAnalysis.graph.init->nextBifurcationPtr[0];
		}
		else {
			fullTreeGraph.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0] = skeletonRightGen->mAnalysis.graph.init;
			fullTreeGraph.init->nextBifurcationPtr[0]->nextBifurcationPtr[0] = skeletonRightGen->mAnalysis.graph.init->nextBifurcationPtr[0];

			fullTreeGraph.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0] = skeletonLeftGen->mAnalysis.graph.init;
			fullTreeGraph.init->nextBifurcationPtr[0]->nextBifurcationPtr[1] = skeletonLeftGen->mAnalysis.graph.init->nextBifurcationPtr[0];
		}

		fullTreeGraph.ggraph2Model(fullTreeModel);
		fullTreeModel->exportToFile(path + "_7_fullTree");
		fullTreeModel->splitSkeletonEdges(0.8);
		fullTreeModel->mAnalysis.generateGraph(*fullTreeModel, 0);
		fullTreeExtendedGraph = fullTreeModel->mAnalysis.graph;

		fullTreeExtendedGraph.ggraph2Model(fullTreeExtendedModel);
		fullTreeExtendedModel->exportToFile(path + "_7_fullTreeOversampled");
	}

	dotObj * treePCExpL = new dotObj();
	dotObj * treePCExpR = new dotObj();

	//BUILD 3D
	if (buildPC) {
		stat->action = "Build 3D Point cloud";
		std::cout << "Build 3D Point cloud" << std::endl;
		sim = new simulation();
		sim->generate3DPointCloudFromTree(tree3DDepth, density, 0.5, *skeletonLeftGen, *hostL, *treePCExpL, "Left", false);
		delete sim;
		sim = new simulation();
		sim->generate3DPointCloudFromTree(tree3DDepth, density, 0.5, *skeletonRightGen, *hostR, *treePCExpR, "Right", false);
		delete sim;
	}

	//BUILD PLY FILE
	if (buildPLY) {
		ggraph rebuilt;
		rebuilt = existingModel->mAnalysis.graph;
		if (rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0]->isLeft) {
			rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0] = skeletonLeftGen->mAnalysis.graph.init;
			rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0] = skeletonRightGen->mAnalysis.graph.init;
		}
		else {
			rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0] = skeletonRightGen->mAnalysis.graph.init;
			rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0] = skeletonLeftGen->mAnalysis.graph.init;
		}

		/*
		dotObj * fullModelCenteline = new dotObj();
		rebuilt.ggraph2Model(fullModelCenteline, false);
		fullModelCenteline->mAnalysis.generateGraph(*fullModelCenteline, 0);
		fullModelCenteline->exportToFile(path + "fullModelCenteline");
		*/

		treePCExpR->exportToXYZ("lungPCRec", true);
		treePCExpL->exportToXYZ("lungPCRec", false);
		trachea->exportToXYZ("lungPCRec", false);
	}

	//Calculation of normals and smoothing
	/*if (false){
	std::cout << "Calculating normals" << std::endl;
	std::string strname = path + "lungPCRec.xyz";
	const char* fname = strname.c_str();
	// Reads a .xyz point set file in points[].
	std::vector<PointEx> pointsL;
	std::ifstream stream(fname);
	if (!stream ||
	!CGAL::read_xyz_points(stream, std::back_inserter(pointsL),
	CGAL::Identity_property_map<PointEx>()))
	{
	std::cerr << "Error: cannot read file " << fname << std::endl;
	return;
	}
	std::vector<PointEx> output;
	//parameters
	const double retain_percentage = 80.0;   // percentage of points to retain.
	const double neighbor_radius = 0.5;   // neighbors size.
	const bool require_uniform_sampling = true;
	const int max_iter_number = 100;
	CGAL::wlop_simplify_and_regularize_point_set
	<Concurrency_tag>
	(pointsL.begin(),
	pointsL.end(),
	std::back_inserter(output),
	retain_percentage,
	neighbor_radius,
	max_iter_number,
	require_uniform_sampling
	);
	std::ofstream out(path + "WLOP_lungPCRec.xyz");
	if (!out || !CGAL::write_xyz_points(
	out, output.begin(), output.end()))
	{
	//return EXIT_FAILURE;
	}
	}*/

	if (buildNormals) {
		stat->action = "Calculating normals";
		std::cout << "Calculating normals" << std::endl;
		std::string strname = "lungPCRec.xyz";
		const char* fname = strname.c_str();
		std::vector<PointEx> pointsL;
		std::ifstream stream(fname);
		if (!stream ||
			!CGAL::read_xyz_points(stream, std::back_inserter(pointsL),
				CGAL::Identity_property_map<PointEx>()))
		{
			std::cerr << "Error: cannot read file " << fname << std::endl;

			return;
		}

		/*
		// simplification by clustering using erase-remove idiom
		double cell_size = 0.001;
		pointsL.erase(CGAL::grid_simplify_point_set(pointsL.begin(), pointsL.end(), cell_size),
		pointsL.end());
		// Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
		std::vector<PointEx>(pointsL).swap(pointsL);*/

		/*CGAL::Timer task_timer; task_timer.start();
		// simplification by clustering using erase-remove idiom
		pointsL.erase(CGAL::hierarchy_simplify_point_set(pointsL.begin(), pointsL.end(),
		5, // Max cluster size
		0.0001), // Max surface variation
		pointsL.end());
		std::size_t memory = CGAL::Memory_sizer().virtual_size();

		std::cout << pointsL.size() << " point(s) kept, computed in "
		<< task_timer.time() << " seconds, "
		<< (memory >> 20) << " Mib allocated." << std::endl;
		*/

		std::list<PointVectorPair> points;

		for (int i = 0; i < pointsL.size(); i++) {
			PointVectorPair temp;
			temp.first = pointsL[i];
			points.push_back(temp);
		}

		/*std::ifstream stream(fname);
		if (!stream ||
		!CGAL::read_xyz_points(stream,
		std::back_inserter(points),
		CGAL::First_of_pair_property_map<PointVectorPair>()))
		{
		std::cerr << "Error: cannot read file " << fname << std::endl;
		return ;
		}*/
		stat->action = "Estimate normals direction";
		std::cout << "Estimate normals direction" << std::endl;
		// Note: pca_estimate_normals() requires an iterator over points
		// as well as property maps to access each point's position and normal.
		const int nb_neighbors = 30; // K-nearest neighbors = 3 rings
		CGAL::pca_estimate_normals<Concurrency_tag>(points.begin(), points.end(),
			CGAL::First_of_pair_property_map<PointVectorPair>(),
			CGAL::Second_of_pair_property_map<PointVectorPair>(),
			nb_neighbors);
		// Orients normals.
		// Note: mst_orient_normals() requires an iterator over points
		// as well as property maps to access each point's position and normal.
		std::list<PointVectorPair>::iterator unoriented_points_begin =
			CGAL::mst_orient_normals(points.begin(), points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>(),
				nb_neighbors);
		// Optional: delete points with an unoriented normal
		// if you plan to call a reconstruction algorithm that expects oriented normals.
		points.erase(unoriented_points_begin, points.end());

		std::ofstream out("lungPC.ply");
		if (!out ||
			!CGAL::write_ply_points_and_normals(
				out, points.begin(), points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>()))
		{
			return;
		}
		stat->action = "Smoothing normals";
		std::cout << "Smoothing normals" << std::endl;
		// Algorithm parameters
		int k = 50;                 // size of neighborhood. The bigger the smoother the result will be.
		// This value should bigger than 1.
		double sharpness_angle = 70; // control sharpness of the result.
		// The bigger the smoother the result will be
		int iter_number = 5;         // number of times the projection is applied

		for (int i = 0; i < iter_number; ++i)
		{
			CGAL::bilateral_smooth_point_set <Concurrency_tag>(
				points.begin(),
				points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>(),
				k,
				sharpness_angle);
		}

		std::ofstream outSmoothedply("lungPCSmoothed.ply");
		std::ofstream outSmoothedxyz("lungPCSmoothed.xyz");
		if (!outSmoothedply ||
			!CGAL::write_ply_points_and_normals(
				outSmoothedply, points.begin(), points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>()))
		{
			return;
		}

		if (!outSmoothedxyz ||
			!CGAL::write_xyz_points_and_normals(
				outSmoothedxyz, points.begin(), points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>()))
		{
			return;
		}
	}

	dotObj * mMOdel = new dotObj;

	if (build3D)
	{
		stat->action = "Poisson Recontruction Commencing";
		std::cout << "Poisson Recontruction Commencing" << std::endl;
		//PoissonRecon(path + "lungPCSmoothed.ply", path + "lungPC.screened.ply",poissonDepth->value());
		SSDRecon("lungPCSmoothed.ply", "lungPC.screened.ply", poissonDepth, mMOdel);
		//poissonRec(path + "lungPCSmoothed.xyz", path + "lungPCscreened.off");
		//isotropicRemesh(path + "lungPCscreened.off", path + "lungPCscreenedIso.off");
		std::cout << "Poisson Recontruction Complete" << std::endl;
		mMOdel->recalculateNormals();

		mMOdel->exportToFile(path + "_9_fullReconstructedModel");

		mMOdel->simplificationEdgeCollapse(1600000);

		mMOdel->recalculateNormals();

		mMOdel->exportToFile(path + "_9_simplifiedModel");

		//if (lungmodel.size() > 0)lungmodel.back() = *mMOdel;

		*existingModel = *mMOdel;
	}

	namespace PMP = CGAL::Polygon_mesh_processing;
	namespace params = PMP::parameters;

	dotObj* cubef = new dotObj();
	cubef->initializeFromFile(path + "cubef.obj");

	std::stringstream modelmeshStream;

	if (refine3D) {
		mMOdel->graphBasedAnalysis(*fullTreeExtendedModel);

		Mesh out, modelmesh;
		bool valid_union = true;
		for (int i = 0; i < mMOdel->mAnalysis.graph.nodes.size(); i++) {
			if (mMOdel->mAnalysis.graph.nodes[i].isTerminal) {
				dotObj * newcube = new dotObj();
				*newcube = *cubef;
				Vector3f direction = mMOdel->mAnalysis.graph.nodes[i].previousNodePtr->position - mMOdel->mAnalysis.graph.nodes[i].position;
				float d = 2 * mMOdel->mAnalysis.graph.nodes[i].diameter;
				if (d > 0) {
					newcube->scale(d, d, d);
					newcube->rotate(direction);
					direction.normalize();
					Vector3f translatePosition = (mMOdel->mAnalysis.graph.nodes[i].position - 0.7*d*direction);
					newcube->translate(translatePosition);
					std::stringstream input = newcube->toOFF();
					if (input) {
						if (!input || !(input >> out))
						{
							std::cerr << "First mesh is not a valid off file." << std::endl;
							return;
						}
					}
				}
			}
			if (mMOdel->mAnalysis.graph.nodes[i].isInlet) {
				dotObj * newcube = new dotObj();
				*newcube = *cubef;
				Vector3f direction = mMOdel->mAnalysis.graph.nodes[i].nextNodesPtr[0]->position - mMOdel->mAnalysis.graph.nodes[i].position;
				float d = 2 * mMOdel->mAnalysis.graph.nodes[i].nextNodesPtr[0]->diameter;
				if (d > 0) {
					newcube->scale(d, d, d);
					newcube->rotate(direction);
					direction.normalize();
					Vector3f translatePosition = (mMOdel->mAnalysis.graph.nodes[i].position - 0.7*d*direction);
					newcube->translate(translatePosition);
					std::stringstream input = newcube->toOFF();
					if (input) {
						if (!input || !(input >> out))
						{
							std::cerr << "First mesh is not a valid off file." << std::endl;
							return;
						}
					}
				}
			}
		}
		if (valid_union)
		{
			std::cout << "Union was successfully computed\n";
			std::ofstream output(path + "union.off");
			output << out;

			std::stringstream model = mMOdel->toOFF();
			if (!model || !(model >> modelmesh))
			{
				std::cerr << "First mesh is not a valid off file." << std::endl;
				return;
			}

			//create a property on edges to indicate whether they are constrained
			Mesh::Property_map<edge_descriptor, bool> is_constrained_map =
				modelmesh.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;
			// update mesh1 to contain the mesh bounding the difference
			// of the two input volumes.
			bool valid_difference =
				PMP::corefine_and_compute_difference(modelmesh,
					out,
					modelmesh,
					params::all_default(), // default parameters for mesh1
					params::all_default(), // default parameters for mesh2
					params::edge_is_constrained_map(is_constrained_map));
			if (valid_difference)
			{
				std::cout << "Difference was successfully computed\n";
				std::ofstream output1(path + "difference.off");
				output1 << modelmesh;
			}
			else {
				std::cout << "Difference could not be computed\n";
				return;
			}

			// collect faces incident to a constrained edge
			std::vector<face_descriptor> selected_faces;
			std::vector<bool> is_selected(num_faces(modelmesh), false);
			BOOST_FOREACH(edge_descriptor e, edges(modelmesh)) {
				if (is_constrained_map[e])
				{
					// insert all faces incident to the target vertex
					BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(halfedge(e, modelmesh), modelmesh))
					{
						if (!is_border(h, modelmesh))
						{
							face_descriptor f = face(h, modelmesh);
							if (!is_selected[f])
							{
								selected_faces.push_back(f);
								is_selected[f] = true;
							}
						}
					}
				}
			}
			// increase the face selection
			CGAL::expand_face_selection(selected_faces, modelmesh, 2, Vector_pmap_wrapper(is_selected), std::back_inserter(selected_faces));
			std::cout << selected_faces.size()
				<< " faces were selected for the remeshing step\n";
			// remesh the region around the intersection polylines
			PMP::isotropic_remeshing(
				selected_faces,
				0.4,
				modelmesh,
				params::edge_is_constrained_map(is_constrained_map));
			//std::ofstream output2("D:/LungModeller/models/difference_remeshed.off");
			//output2 << modelmesh;
			modelmeshStream << modelmesh;
		}
	}

	if (segmentModel) {
		std::string type = "OFF";
		dotObj * Reconstructed = new dotObj();
		std::cout << "Segment Model" << std::endl;
		Reconstructed->initializeFromFile(modelmeshStream, type);
		Reconstructed->recalculateNormals();
		Reconstructed->segmentationMode = "NORMALS";
		Reconstructed->segmentByGeneration(*fullTreeExtendedModel);
		Reconstructed->exportToFileSegmentedPerFace(path + "final_segmented_model_2.obj");
		*existingModel = *Reconstructed;
	}

	stat->isComplete = true;

	stat->error = false;

	/*
	//REFINE
	if (refine3D){
	sim = new simulation();

	std::cout << "Segmentation" << std::endl;
	//Reconstructed
	dotObj * Reconstructed = new dotObj();
	Reconstructed->initializeFromPLY(path + "lungPC.screened.ply");
	std::cout << "Recalculate normals" << std::endl;
	Reconstructed->recalculateNormals();
	Reconstructed->exportToFile(path + "lungReconstructed");
	Reconstructed->smoothingSparse("taubin", "weighted", 100, -0.7, 0.55);
	//Reconstructed->smoothingSparse("taubin", "weighted", 200, -0.7, 0.5);
	std::cout << "Recalculate normals" << std::endl;
	Reconstructed->recalculateNormals();
	Reconstructed->exportToFile(path + "lungReconstructedSmooth");

	//Left 1D
	std::string oneDR = path + "skeletonRightGenerated.obj";
	dotObj *oneDModelR = new dotObj();
	oneDModelR->initializeFromFile(oneDR);

	//Right 1D
	std::string oneDL = path + "skeletonLeftGenerated.obj";
	dotObj *oneDModelL = new dotObj();
	oneDModelL->initializeFromFile(oneDL);

	sim->flatten(*oneDModelL, *oneDModelR, *Reconstructed, path + "lungReconstructedFlattened", "obj");

	delete oneDModelR;
	delete oneDModelL;

	delete sim;
	}

	if (segmentModel){
	sim = new simulation();

	std::cout << "Segmentation" << std::endl;
	//Reconstructed
	dotObj * Reconstructed = new dotObj();
	//Reconstructed->initializeFromFile(path + "lungReconstructedFlattened.obj");
	Reconstructed->initializeFromFile(path + "lungReconstructedSmooth.obj");
	//Left 1D
	std::string oneDR = path + "skeletonRightGenerated.obj";
	dotObj *oneDModelR = new dotObj();
	oneDModelR->initializeFromFile(oneDR);

	//Right 1D
	std::string oneDL = path + "skeletonLeftGenerated.obj";
	dotObj *oneDModelL = new dotObj();
	oneDModelL->initializeFromFile(oneDL);

	sim->segmentTree2(*oneDModelL, *oneDModelR, *Reconstructed, path + "lungReconstructedSegmented", "full");
	delete oneDModelR;
	delete oneDModelL;
	//delete recModel;
	delete sim;
	}*/

	return;
}

//===================================================
// CT based version test
//===================================================
void simulation::extendBronchialTreeV2(
	bool buildVolumes,
	bool buildCenterline,
	bool useCline,
	bool build1DModel,
	bool surfaceSampling,
	bool buildNormals,
	bool build3DModel,
	bool refinements,
	bool segmentation,
	int depth,
	int volumeDepth,
	int density,
	int poissonDepth,
	std::string path,
	dotObj * boundaryR,
	dotObj * boundaryL,
	dotObj * trachea,
	dotObj * existingModel,
	dotObj * existingCenterline,
	status * stat,
	bool verbose,
	volume * volLeft = nullptr,
	volume * volRight = nullptr
) {
	bool initvols = buildVolumes;
	bool initskeletonization = buildCenterline;
	bool init1D = build1DModel;
	bool initPC = surfaceSampling;
	bool initNormals = buildNormals;
	bool init3D = build3DModel;
	bool initRefine = refinements;
	bool initSegmentation = segmentation;

	bool useCenterLineInstead = useCline;

	bool preSmoothing = initskeletonization;
	bool buildSkel = initskeletonization;
	//1D
	bool doMatch = init1D;
	bool build1D = init1D;
	//PC
	bool buildPC = initPC;
	bool buildPLY = initPC;
	//Normals & Smoothing
	//Reconstruction
	bool build3D = init3D;
	//Refinements
	bool refine3D = initRefine;
	//Segmentation
	bool segmentModel = initSegmentation;
	//Depth
	int tree3DDepth = 6; //Caution

	int treeDepth = depth;
	if (existingModel->vertices.size() > 0 && existingModel->faces.size() && existingModel->normals.size() > 0 &&
		boundaryL->vertices.size() > 0 && boundaryL->faces.size() && boundaryL->normals.size() > 0 &&
		boundaryR->vertices.size() > 0 && boundaryR->faces.size() && boundaryR->normals.size() > 0) {
	}
	else {
		std::cout << "Error: Enter valid 3D models" << std::endl;
		//return;
	}
	stat->action = "Simulation :: Commencing";
	std::cout << std::endl << "Simulation :: Commencing" << std::endl;
	simulation *sim;
	dotObj * skelLeftGen;
	dotObj * skelLeftRight;
	dotObj cachedModel, skel;
	dotObj * hostR = new dotObj();
	dotObj * hostL = new dotObj();
	dotObj * recon = new dotObj();
	dotObj * fullTree = new dotObj();
	trachea = new dotObj;
	dotObj * skeletonLeftGen = new dotObj();
	dotObj * skeletonRightGen = new dotObj();

	//BUILD VOLUMES
	if (buildVolumes) {
		stat->action = "Volume generation";
		skelLeftGen = new dotObj();
		skelLeftRight = new dotObj();

		if ((volLeft == nullptr) || (volRight == nullptr)) {
			volLeft = new volume();
			volRight = new volume();
		}
		sim = new simulation();
		sim->generateVolume(volumeDepth, density, *volLeft, *boundaryL, "Left");
		delete sim;

		sim = new simulation();
		sim->generateVolume(volumeDepth, density, *volRight, *boundaryR, "Right");
		delete sim;
	}

	//BUILD SKELETON
	if (buildSkel) {
		stat->action = "Build Skeleton";

		if (!useCenterLineInstead) {
			cachedModel = *existingModel;
			//cachedModel.smoothingSparse("taubin", "onoff", 2000, -0.95, 0.35);
			//cachedModel.exportToFile(path + "smoothed");
			skel = cachedModel.mcfskel();
			skel.exportToFile(path + "_1_skeletonRaw");
			skel.mcfskelRefineStepOne(8);
			skel.mcfskelRefineStepOne(8);
			skel.mcfskelRefineStepOne(8);
			skel.mcfskelRefineStepTwo();
			skel.exportToFile(path + "_2_skeletonRefined");
		}
		else {
			skel = *existingCenterline;
			skel.exportToFile(path + "_1_skeletonRaw");
			skel.mcfskelRefineStepOne(8);
			skel.mcfskelRefineStepOne(8);
			skel.mcfskelRefineStepOne(8);
			skel.mcfskelRefineStepTwo();
			skel.mcfskelRefineStepTwo();
			skel.exportToFile(path + "_2_skeletonRefined");
		}

		/*
		ggraph rebuilt;
		//existingModel->mAnalysis.parseTerminals = true;
		//existingModel->graphBasedAnalysis(skel);
		skel.graphBased1DModelAnalysis();
		rebuilt = skel.mAnalysis.graph;
		existingModel->mAnalysis.graph = rebuilt;
		*/

		//dotObj skel;

		recon = new dotObj();
		trachea = new dotObj();
		existingModel->mAnalysis.parseTerminals = true;
		existingModel->graphBasedAnalysis(skel);
		std::cout << "Initialization index:" << skel.mAnalysis.inlet << std::endl;

		/*
		existingModel->mAnalysis.graph.init->nextBifurcationPtr[0]->
			nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0]->doStop = true;
		existingModel->mAnalysis.graph.init->nextBifurcationPtr[0]->
			nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0]->doStop = true;
		*/

		/*existingModel->mAnalysis.graph.init->nextBifurcationPtr[0]->
			nextBifurcationPtr[0]->
			nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0]->doStop = true;

		existingModel->mAnalysis.graph.init->nextBifurcationPtr[0]->nextBifurcationPtr[1]->
			nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0]->doStop = true;

		existingModel->mAnalysis.graph.init->nextBifurcationPtr[0]-> nextBifurcationPtr[0]->
			nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0]->doStop = true;

		existingModel->mAnalysis.graph.init->nextBifurcationPtr[0]->nextBifurcationPtr[1]->
			nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0]->doStop = true;*/

		existingModel->mAnalysis.graph.init->nextBifurcationPtr[0]->
			nextBifurcationPtr[0]->
			nextBifurcationPtr[0]->doStop = true;
		existingModel->mAnalysis.graph.init->nextBifurcationPtr[0]->
			nextBifurcationPtr[0]->
			nextBifurcationPtr[1]->doStop = true;
		existingModel->mAnalysis.graph.init->nextBifurcationPtr[0]->
			nextBifurcationPtr[1]->
			nextBifurcationPtr[0]->doStop = true;
		existingModel->mAnalysis.graph.init->nextBifurcationPtr[0]->
			nextBifurcationPtr[1]->
			nextBifurcationPtr[1]->doStop = true;

		existingModel->mAnalysis.graph.ggraph2Model(recon);
		existingModel->mAnalysis.graph.ggraph2ModelV(trachea);
		recon->exportToFile(path + "_3_trachea_centerline");
		trachea->exportToFile(path + "_4_trachea");
		delete(recon);

		//existingModel->mAnalysis.parseTerminals = true;
		//existingModel->graphBasedAnalysis(skel);
		/*
		rebuilt = skel.mAnalysis.graph;
		recon = new dotObj();
		std::cout << "Initialization index:"<< skel.mAnalysis.inlet << std::endl;
		rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0]->doStop = true;
		rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0]->doStop = true;
		rebuilt.ggraph2Model(recon);
		recon->exportToFile(path + "_3_trachea_centerline");
		delete(recon);
		*/

		//==========================================================
		ggraph rebuilt;
		//skel.graphBased1DModelAnalysis();
		//rebuilt = skel.mAnalysis.graph;

		existingModel->graphBasedAnalysis(skel);
		rebuilt = existingModel->mAnalysis.graph;

		dotObj skel2;
		std::vector<gnode*> passed = rebuilt.init->nextBifurcationPtr;
		for (int i = 0; i < passed.size(); i++) {
			std::vector<gnode*> passed2 = passed[i]->nextBifurcationPtr;
			for (int j = 0; j < passed2.size(); j++) {
				std::vector<gnode*> passed3 = passed2[j]->nextBifurcationPtr;
				for (int k = 0; k < passed3.size(); k++) {
					passed3[k]->doStop = true;
					//std::vector<gnode*> passed4 = passed3[j]->nextBifurcationPtr;
					//for (int k = 0; k < passed4.size(); k++) {
						//passed4[k]->doStop = true;
						//std::vector<gnode*> passed5 = passed4[j]->nextBifurcationPtr;
						//for (int k = 0; k < passed5.size(); k++) {
						//	passed5[k]->doStop = true;
						//}
					//}
				}
			}
		}
		rebuilt.ggraph2Model(&skel2);
		skel2.exportToFile(path + "_3b_nMod");
		existingModel->graphBasedAnalysis(skel2);
		rebuilt = existingModel->mAnalysis.graph;
		//skel2.graphBased1DModelAnalysis();
		//rebuilt = skel2.mAnalysis.graph;

		//===========================================================

		int genLim = 5;

		gnode * g;
		//=============================================================

		for (int i = 0; i < rebuilt.nodes.size(); i++) {
			if (rebuilt.nodes[i].generation > genLim)rebuilt.nodes[i].doStop = true;
		}
		g = new gnode();
		hostL = new dotObj();
		g = rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0];
		if (g->isLeft) {
			g = rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr;
		}
		else {
			g = rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr;
		}
		rebuilt.init = g;
		rebuilt.ggraph2Model(hostL, true);
		hostL->mAnalysis.generateGraph(*hostL, 0);
		hostL->mAnalysis.graph.init->diameter = g->diameter;
		hostL->exportToFile(path + "_5_host_centerline_left");
		//hostL->splitSkeletonEdges(0.8);
		//hostL->exportToFile(path + "_5_host_centerline_left_oversampled");

		//=============================================================
		//rebuilt = skel2.mAnalysis.graph;
		rebuilt = existingModel->mAnalysis.graph;
		for (int i = 0; i < rebuilt.nodes.size(); i++) {
			if (rebuilt.nodes[i].generation > genLim)rebuilt.nodes[i].doStop = true;
		}
		g = new gnode();
		hostR = new dotObj();
		g = rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0];
		if (g->isRight) {
			g = rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr;
		}
		else {
			g = rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr;
		}

		rebuilt.init = g;
		rebuilt.ggraph2Model(hostR, true);
		hostR->mAnalysis.generateGraph(*hostR, 0);
		hostR->mAnalysis.graph.init->diameter = g->diameter;
		hostR->exportToFile(path + "_5_host_centerline_right");
		//hostR->splitSkeletonEdges(0.8);
		//hostR->exportToFile(path + "_5_host_centerline_right_oversampled");
	}

	//BUILD 1D

	dotObj * fullTreeExtendedModel = new dotObj();
	ggraph fullTreeExtendedGraph;
	ggraph fullTreeGraph;

	if (build1D) {
		stat->action = "Extending geometry";
		std::cout << "Extending geometry" << std::endl;
		sim = new simulation();
		sim->generate1DTree(treeDepth, *volLeft, *boundaryL, *hostL, *skeletonLeftGen);
		//skeletonLeftGen->splitSkeletonEdges(0.5);
		if (verbose) skeletonLeftGen->exportToFile(path + "_6_skeletonLeftGenerated");
		delete sim;
		sim = new simulation();
		sim->generate1DTree(treeDepth, *volRight, *boundaryR, *hostR, *skeletonRightGen);
		//skeletonRightGen->splitSkeletonEdges(0.5);
		if (verbose) skeletonRightGen->exportToFile(path + "_6_skeletonRightGenerated");
		delete sim;


		//skeletonLeftGen->mAnalysis.generateGraph(*skeletonLeftGen, 0);

		std::cout << "skeletonLeftGen->graphBased1DModelAnalysis();" << std::endl;
		skeletonLeftGen->graphBased1DModelAnalysis(true, true,0);
		skeletonLeftGen->mAnalysis.graph.init->diameter = hostL->mAnalysis.graph.init->diameter;
		

		std::cout << "skeletonRightGen->graphBased1DModelAnalysis();" << std::endl;
		skeletonRightGen->graphBased1DModelAnalysis(true, true,0);
		skeletonRightGen->mAnalysis.graph.init->diameter = hostR->mAnalysis.graph.init->diameter;

		//Full tree rebuilt

		existingModel->mAnalysis.parseTerminals = true;
		existingModel->graphBasedAnalysis(skel);

		//skel.graphBased1DModelAnalysis(false);
		//skel.mAnalysis.graph.exportGraphFeatures(path + "graph.csv");

		fullTreeGraph = existingModel->mAnalysis.graph;
		fullTreeGraph.initializeDiameters();
		fullTreeGraph.generateGraphFeaturesLRDiscrimination();
		fullTreeGraph.exportGraphFeatures(path + "existingModel.json");

		if (fullTreeGraph.init->nextBifurcationPtr[0]->nextNodesPtr[0]->isLeft) {
			fullTreeGraph.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0] = skeletonLeftGen->mAnalysis.graph.init;
			fullTreeGraph.init->nextBifurcationPtr[0]->nextBifurcationPtr[0] = skeletonLeftGen->mAnalysis.graph.init->nextBifurcationPtr[0];
			skeletonLeftGen->mAnalysis.graph.init->nextBifurcationPtr[0]->previousBifurcationPtr = fullTreeGraph.init->nextBifurcationPtr[0];

			fullTreeGraph.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0] = skeletonRightGen->mAnalysis.graph.init;
			fullTreeGraph.init->nextBifurcationPtr[0]->nextBifurcationPtr[1] = skeletonRightGen->mAnalysis.graph.init->nextBifurcationPtr[0];
			skeletonRightGen->mAnalysis.graph.init->nextBifurcationPtr[0]->previousBifurcationPtr = fullTreeGraph.init->nextBifurcationPtr[0];
		}
		else {
			fullTreeGraph.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0] = skeletonRightGen->mAnalysis.graph.init;
			fullTreeGraph.init->nextBifurcationPtr[0]->nextBifurcationPtr[0] = skeletonRightGen->mAnalysis.graph.init->nextBifurcationPtr[0];
			skeletonRightGen->mAnalysis.graph.init->nextBifurcationPtr[0]->previousBifurcationPtr = fullTreeGraph.init->nextBifurcationPtr[0];

			fullTreeGraph.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0] = skeletonLeftGen->mAnalysis.graph.init;
			fullTreeGraph.init->nextBifurcationPtr[0]->nextBifurcationPtr[1] = skeletonLeftGen->mAnalysis.graph.init->nextBifurcationPtr[0];
			skeletonLeftGen->mAnalysis.graph.init->nextBifurcationPtr[0]->previousBifurcationPtr = fullTreeGraph.init->nextBifurcationPtr[0];
		}






		std::cout << "---------------------fullTreeGraph:Init analysis--------------------" << std::endl;
		fullTreeGraph.initializeDiameters();
		fullTreeGraph.generateGraphFeaturesGeneration();
		fullTreeGraph.generateGraphFeaturesDiameter();
		//std::cout << "fullTreeGraph-->LR" << std::endl;
		//fullTreeGraph.generateGraphFeaturesLRDiscrimination();
		//fullTreeGraph.enrichGraph(0.5);
		fullTreeGraph.generateGraphFeaturesLength("branch");
		fullTreeGraph.generateGraphFeaturesLength("node");
		fullTreeGraph.generateGraphFeaturesRatios();
		fullTreeGraph.propagateGraphFeatures();
		fullTreeGraph.exportGraphFeatures(path + "_7_fullTree.json");
		
		/*gnode * RU, RM, RL, LU, LL;
		if (fullTreeGraph.init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->isLeft) {
		}
		if (fullTreeGraph.init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->isRight) {
			if((fullTreeGraph.init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->nextBifurcationPtr[0]->position- fullTreeGraph.init->position).norm()
				<)
		}else {
		}

		fullTreeGraph.exportGraphFeatures(path + "_7_fullTree_RU.json", fullTreeGraph.init->nextBifurcationPtr[0]);
		fullTreeGraph.exportGraphFeatures(path + "_7_fullTree_RM.json");
		fullTreeGraph.exportGraphFeatures(path + "_7_fullTree_RL.json");
		fullTreeGraph.exportGraphFeatures(path + "_7_fullTree_LU.json");
		fullTreeGraph.exportGraphFeatures(path + "_7_fullTree_LL.json");*/

		//skeletonRightGen->mAnalysis.graph.exportGraphFeatures(path + "skeletonRightGen.csv");
		//skeletonLeftGen->mAnalysis.graph.exportGraphFeatures(path + "skeletonLeftGen.csv");

		dotObj * fullTreeModel = new dotObj();
		fullTreeGraph.ggraph2Model(fullTreeModel);
		fullTreeModel->exportToFile(path + "_7_fulltree");

		/*fullTreeModel->splitSkeletonEdges(0.5);
		fullTreeModel->graphBased1DModelAnalysis(); //<==Horsfield & Strahler
		fullTreeExtendedGraph = fullTreeModel->mAnalysis.graph;
		fullTreeExtendedGraph.exportGraphFeatures(path + "_7_fullTreeOversampled.csv");
		fullTreeExtendedGraph.ggraph2Model(fullTreeExtendedModel);
		fullTreeExtendedModel->exportToFile(path + "_7_fullTreeOversampled");
		dotObj * fullTreeModelColored = new dotObj();
		fullTreeExtendedGraph.ggraph2ModelColorize(fullTreeModelColored);
		fullTreeModelColored->exportToFileSegmentedPerLine("obj", path + "_7_fulltreeSegmented");*/
	}

	dotObj * treePCExpL = new dotObj();
	dotObj * treePCExpR = new dotObj();

	//BUILD 3D
	if (buildPC) {
		/*stat->action = "Build 3D Point cloud";
		std::cout << "Build 3D Point cloud" << std::endl;
		sim = new simulation();
		sim->generate3DPointCloudFromTree(tree3DDepth, density, 0.5, *skeletonLeftGen, *hostL, *treePCExpL, "Left", true);
		delete sim;
		sim = new simulation();
		sim->generate3DPointCloudFromTree(tree3DDepth, density, 0.5, *skeletonRightGen, *hostR, *treePCExpR, "Right", true);
		delete sim;*/

		dotObj * solid = new dotObj();
		sim = new simulation();
		sim->generate3DPointCloudFromGraph(tree3DDepth, density, 0.5, &fullTreeGraph, solid);
		delete sim;
		solid->exportToFile(path + "_8_fulltree");
		solid->exportToXYZ(path + "lungPCRec");
	}

	//BUILD PLY FILE
	if (buildPLY) {
		/*ggraph rebuilt;
		rebuilt = existingModel->mAnalysis.graph;
		if (rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0]->isLeft) {
			rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0] = skeletonLeftGen->mAnalysis.graph.init;
			rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0] = skeletonRightGen->mAnalysis.graph.init;
		}
		else {
			rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->nextNodesPtr[0] = skeletonRightGen->mAnalysis.graph.init;
			rebuilt.init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->nextNodesPtr[0] = skeletonLeftGen->mAnalysis.graph.init;
		}
		treePCExpR->exportToXYZ(path + "lungPCRec", true);
		treePCExpL->exportToXYZ(path + "lungPCRec", false);
		trachea->exportToXYZ(path + "lungPCRec", false);*/
	}

	//Calculation of normals and smoothing
	/*if (false){
	std::cout << "Calculating normals" << std::endl;
	std::string strname = path + "lungPCRec.xyz";
	const char* fname = strname.c_str();
	// Reads a .xyz point set file in points[].
	std::vector<PointEx> pointsL;
	std::ifstream stream(fname);
	if (!stream ||
	!CGAL::read_xyz_points(stream, std::back_inserter(pointsL),
	CGAL::Identity_property_map<PointEx>()))
	{
	std::cerr << "Error: cannot read file " << fname << std::endl;
	return;
	}
	std::vector<PointEx> output;
	//parameters
	const double retain_percentage = 80.0;   // percentage of points to retain.
	const double neighbor_radius = 0.5;   // neighbors size.
	const bool require_uniform_sampling = true;
	const int max_iter_number = 100;
	CGAL::wlop_simplify_and_regularize_point_set
	<Concurrency_tag>
	(pointsL.begin(),
	pointsL.end(),
	std::back_inserter(output),
	retain_percentage,
	neighbor_radius,
	max_iter_number,
	require_uniform_sampling
	);
	std::ofstream out(path + "WLOP_lungPCRec.xyz");
	if (!out || !CGAL::write_xyz_points(
	out, output.begin(), output.end()))
	{
	//return EXIT_FAILURE;
	}
	}*/

	if (buildNormals) {
		stat->action = "Calculating normals";
		std::cout << "Calculating normals" << std::endl;
		std::string strname = path + "lungPCRec.xyz";
		const char* fname = strname.c_str();
		std::vector<PointEx> pointsL;
		std::ifstream stream(fname);
		if (!stream ||
			!CGAL::read_xyz_points(stream, std::back_inserter(pointsL),
				CGAL::Identity_property_map<PointEx>()))
		{
			std::cerr << "Error: cannot read file " << fname << std::endl;

			return;
		}

		/*
		// simplification by clustering using erase-remove idiom
		double cell_size = 0.001;
		pointsL.erase(CGAL::grid_simplify_point_set(pointsL.begin(), pointsL.end(), cell_size),
		pointsL.end());
		// Optional: after erase(), use Scott Meyer's "swap trick" to trim excess capacity
		std::vector<PointEx>(pointsL).swap(pointsL);*/

		/*CGAL::Timer task_timer; task_timer.start();
		// simplification by clustering using erase-remove idiom
		pointsL.erase(CGAL::hierarchy_simplify_point_set(pointsL.begin(), pointsL.end(),
		5, // Max cluster size
		0.0001), // Max surface variation
		pointsL.end());
		std::size_t memory = CGAL::Memory_sizer().virtual_size();

		std::cout << pointsL.size() << " point(s) kept, computed in "
		<< task_timer.time() << " seconds, "
		<< (memory >> 20) << " Mib allocated." << std::endl;
		*/

		std::list<PointVectorPair> points;

		for (int i = 0; i < pointsL.size(); i++) {
			PointVectorPair temp;
			temp.first = pointsL[i];
			points.push_back(temp);
		}

		/*std::ifstream stream(fname);
		if (!stream ||
		!CGAL::read_xyz_points(stream,
		std::back_inserter(points),
		CGAL::First_of_pair_property_map<PointVectorPair>()))
		{
		std::cerr << "Error: cannot read file " << fname << std::endl;
		return ;
		}*/
		stat->action = "Estimate normals direction";
		std::cout << "Estimate normals direction" << std::endl;
		// Note: pca_estimate_normals() requires an iterator over points
		// as well as property maps to access each point's position and normal.
		const int nb_neighbors = 30; // K-nearest neighbors = 3 rings
		CGAL::pca_estimate_normals<Concurrency_tag>(points.begin(), points.end(),
			CGAL::First_of_pair_property_map<PointVectorPair>(),
			CGAL::Second_of_pair_property_map<PointVectorPair>(),
			nb_neighbors);
		// Orients normals.
		// Note: mst_orient_normals() requires an iterator over points
		// as well as property maps to access each point's position and normal.
		std::list<PointVectorPair>::iterator unoriented_points_begin =
			CGAL::mst_orient_normals(points.begin(), points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>(),
				nb_neighbors);
		// Optional: delete points with an unoriented normal
		// if you plan to call a reconstruction algorithm that expects oriented normals.
		points.erase(unoriented_points_begin, points.end());

		std::ofstream out("lungPC.ply");
		if (!out ||
			!CGAL::write_ply_points_and_normals(
				out, points.begin(), points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>()))
		{
			return;
		}
		stat->action = "Smoothing normals";
		std::cout << "Smoothing normals" << std::endl;
		// Algorithm parameters
		int k = 50;                 // size of neighborhood. The bigger the smoother the result will be.
		// This value should bigger than 1.
		double sharpness_angle = 70; // control sharpness of the result.
		// The bigger the smoother the result will be
		int iter_number = 5;         // number of times the projection is applied

		for (int i = 0; i < iter_number; ++i)
		{
			CGAL::bilateral_smooth_point_set <Concurrency_tag>(
				points.begin(),
				points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>(),
				k,
				sharpness_angle);
		}

		std::ofstream outSmoothedply("lungPCSmoothed.ply");
		std::ofstream outSmoothedxyz("lungPCSmoothed.xyz");
		if (!outSmoothedply ||
			!CGAL::write_ply_points_and_normals(
				outSmoothedply, points.begin(), points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>()))
		{
			return;
		}

		if (!outSmoothedxyz ||
			!CGAL::write_xyz_points_and_normals(
				outSmoothedxyz, points.begin(), points.end(),
				CGAL::First_of_pair_property_map<PointVectorPair>(),
				CGAL::Second_of_pair_property_map<PointVectorPair>()))
		{
			return;
		}
	}

	dotObj * mMOdel = new dotObj;

	if (build3D)
	{
		stat->action = "Poisson Recontruction Commencing";
		std::cout << "Poisson Recontruction Commencing" << std::endl;
		//PoissonRecon(path + "lungPCSmoothed.ply", path + "lungPC.screened.ply",poissonDepth->value());
		SSDRecon("lungPCSmoothed.ply", "lungPC.screened.ply", poissonDepth, mMOdel);
		//poissonRec(path + "lungPCSmoothed.xyz", path + "lungPCscreened.off");
		//isotropicRemesh(path + "lungPCscreened.off", path + "lungPCscreenedIso.off");
		std::cout << "Poisson Recontruction Complete" << std::endl;
		mMOdel->recalculateNormals();

		mMOdel->exportToFile(path + "_9_fullReconstructedModel");

		mMOdel->simplificationEdgeCollapse(1600000);

		mMOdel->recalculateNormals();

		mMOdel->exportToFile(path + "_9_simplifiedModel");

		//if (lungmodel.size() > 0)lungmodel.back() = *mMOdel;

		*existingModel = *mMOdel;
	}

	namespace PMP = CGAL::Polygon_mesh_processing;
	namespace params = PMP::parameters;

	dotObj* cubef = new dotObj();
	cubef->initializeFromFile(path + "cubef.obj");

	std::stringstream modelmeshStream;

	if (refine3D) {
		mMOdel->graphBasedAnalysis(*fullTreeExtendedModel);

		Mesh out, modelmesh;
		bool valid_union = true;
		for (int i = 0; i < mMOdel->mAnalysis.graph.nodes.size(); i++) {
			if (mMOdel->mAnalysis.graph.nodes[i].isTerminal) {
				dotObj * newcube = new dotObj();
				*newcube = *cubef;
				Vector3f direction = mMOdel->mAnalysis.graph.nodes[i].previousNodePtr->position - mMOdel->mAnalysis.graph.nodes[i].position;
				float d = 2 * mMOdel->mAnalysis.graph.nodes[i].diameter;
				if (d > 0) {
					newcube->scale(d, d, d);
					newcube->rotate(direction);
					direction.normalize();
					Vector3f translatePosition = (mMOdel->mAnalysis.graph.nodes[i].position - 0.7*d*direction);
					newcube->translate(translatePosition);
					std::stringstream input = newcube->toOFF();
					if (input) {
						if (!input || !(input >> out))
						{
							std::cerr << "First mesh is not a valid off file." << std::endl;
							return;
						}
					}
				}
			}
			if (mMOdel->mAnalysis.graph.nodes[i].isInlet) {
				dotObj * newcube = new dotObj();
				*newcube = *cubef;
				Vector3f direction = mMOdel->mAnalysis.graph.nodes[i].nextNodesPtr[0]->position - mMOdel->mAnalysis.graph.nodes[i].position;
				float d = 2 * mMOdel->mAnalysis.graph.nodes[i].nextNodesPtr[0]->diameter;
				if (d > 0) {
					newcube->scale(d, d, d);
					newcube->rotate(direction);
					direction.normalize();
					Vector3f translatePosition = (mMOdel->mAnalysis.graph.nodes[i].position - 0.7*d*direction);
					newcube->translate(translatePosition);
					std::stringstream input = newcube->toOFF();
					if (input) {
						if (!input || !(input >> out))
						{
							std::cerr << "First mesh is not a valid off file." << std::endl;
							return;
						}
					}
				}
			}
		}
		if (valid_union)
		{
			std::cout << "Union was successfully computed\n";
			std::ofstream output(path + "union.off");
			output << out;

			std::cout << "mMOdel->toOFF()" << std::endl;
			std::stringstream model = mMOdel->toOFF();
			if (!model || !(model >> modelmesh))
			{
				std::cerr << "First mesh is not a valid off file." << std::endl;
				return;
			}

			std::cout << "Create a property on edges to indicate whether they are constrained" << std::endl;

			Mesh::Property_map<edge_descriptor, bool> is_constrained_map =
				modelmesh.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;
			// update mesh1 to contain the mesh bounding the difference
			// of the two input volumes.
			bool valid_difference =
				PMP::corefine_and_compute_difference(modelmesh,
					out,
					modelmesh,
					params::all_default(), // default parameters for mesh1
					params::all_default(), // default parameters for mesh2
					params::edge_is_constrained_map(is_constrained_map));
			if (valid_difference)
			{
				std::cout << "Difference was successfully computed\n";
				std::ofstream output1(path + "difference.off");
				output1 << modelmesh;
			}
			else {
				std::cout << "Difference could not be computed\n";
				return;
			}

			// collect faces incident to a constrained edge
			std::vector<face_descriptor> selected_faces;
			std::vector<bool> is_selected(num_faces(modelmesh), false);
			BOOST_FOREACH(edge_descriptor e, edges(modelmesh)) {
				if (is_constrained_map[e])
				{
					// insert all faces incident to the target vertex
					BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(halfedge(e, modelmesh), modelmesh))
					{
						if (!is_border(h, modelmesh))
						{
							face_descriptor f = face(h, modelmesh);
							if (!is_selected[f])
							{
								selected_faces.push_back(f);
								is_selected[f] = true;
							}
						}
					}
				}
			}
			// increase the face selection
			CGAL::expand_face_selection(selected_faces, modelmesh, 2, Vector_pmap_wrapper(is_selected), std::back_inserter(selected_faces));
			std::cout << selected_faces.size()
				<< " faces were selected for the remeshing step\n";
			// remesh the region around the intersection polylines
			PMP::isotropic_remeshing(
				selected_faces,
				0.4,
				modelmesh,
				params::edge_is_constrained_map(is_constrained_map));
			//std::ofstream output2("D:/LungModeller/models/difference_remeshed.off");
			//output2 << modelmesh;
			modelmeshStream << modelmesh;
		}
	}

	if (segmentModel) {
		std::string type = "OFF";
		dotObj * Reconstructed = new dotObj();

		Reconstructed->initializeFromFile(modelmeshStream, type);
		Reconstructed->recalculateNormals();
		Reconstructed->segmentationMode = "NORMALS";
		Reconstructed->segmentByGeneration(*fullTreeExtendedModel);
		Reconstructed->exportToFileSegmentedPerFace(path + "final_segmented_model_2.obj");
		*existingModel = *Reconstructed;
	}

	stat->isComplete = true;

	stat->error = false;

	/*
	//REFINE
	if (refine3D){
	sim = new simulation();

	std::cout << "Segmentation" << std::endl;
	//Reconstructed
	dotObj * Reconstructed = new dotObj();
	Reconstructed->initializeFromPLY(path + "lungPC.screened.ply");
	std::cout << "Recalculate normals" << std::endl;
	Reconstructed->recalculateNormals();
	Reconstructed->exportToFile(path + "lungReconstructed");
	Reconstructed->smoothingSparse("taubin", "weighted", 100, -0.7, 0.55);
	//Reconstructed->smoothingSparse("taubin", "weighted", 200, -0.7, 0.5);
	std::cout << "Recalculate normals" << std::endl;
	Reconstructed->recalculateNormals();
	Reconstructed->exportToFile(path + "lungReconstructedSmooth");

	//Left 1D
	std::string oneDR = path + "skeletonRightGenerated.obj";
	dotObj *oneDModelR = new dotObj();
	oneDModelR->initializeFromFile(oneDR);

	//Right 1D
	std::string oneDL = path + "skeletonLeftGenerated.obj";
	dotObj *oneDModelL = new dotObj();
	oneDModelL->initializeFromFile(oneDL);

	sim->flatten(*oneDModelL, *oneDModelR, *Reconstructed, path + "lungReconstructedFlattened", "obj");

	delete oneDModelR;
	delete oneDModelL;

	delete sim;
	}

	if (segmentModel){
	sim = new simulation();

	std::cout << "Segmentation" << std::endl;
	//Reconstructed
	dotObj * Reconstructed = new dotObj();
	//Reconstructed->initializeFromFile(path + "lungReconstructedFlattened.obj");
	Reconstructed->initializeFromFile(path + "lungReconstructedSmooth.obj");
	//Left 1D
	std::string oneDR = path + "skeletonRightGenerated.obj";
	dotObj *oneDModelR = new dotObj();
	oneDModelR->initializeFromFile(oneDR);

	//Right 1D
	std::string oneDL = path + "skeletonLeftGenerated.obj";
	dotObj *oneDModelL = new dotObj();
	oneDModelL->initializeFromFile(oneDL);

	sim->segmentTree2(*oneDModelL, *oneDModelR, *Reconstructed, path + "lungReconstructedSegmented", "full");
	delete oneDModelR;
	delete oneDModelL;
	//delete recModel;
	delete sim;
	}*/

	return;
}