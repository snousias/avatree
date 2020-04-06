#include "lungModelling.h"
#include <iostream>

#include "itkImage.h"
#include <itkImageFileWriter.h>
#include <itkImageFileReader.h>
#include <itkMetaImageIO.h>
#include <itkImageRegionIterator.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/selection.h>
#include <fstream>


/*
namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

int generate(void) {
	bool buildVolumes = false;
	bool buildCenterline = true;
	bool build1DModel = true;
	bool surfaceSampling = false;
	bool buildNormals = false;
	bool build3DModel = false;

	bool refinements = false;
	bool segmentation = false;
	int depth = 7;
	int volumeDepth = 16;
	int density = 500;
	int poissonDepth = 11;

	//std::string path = "D:/LungModeller/models/Vessel02/";
	std::string path = "E:/_Groundwork/Lung/LungVis/Data/VESSEL12_02/";

	//std::string boundaryMeshL = path + "L1.obj";
	//std::string boundaryMeshR = path + "L2.obj";
	//std::string existingModelMesh = path + "lungPartsFullB.obj";
	//std::string existingModelMesh = path + "trachea-down-final.obj";

	//std::string boundaryMeshL = path + "V2_LL_simplified_01.obj";
	//std::string boundaryMeshR = path + "V2_LR_simplified_01.obj";
	//std::string existingModelMesh = path + "V2_mainAirways.obj";

	std::string boundaryMeshL = path + "VESSEL12_02_Lungs_L_3.obj";
	std::string boundaryMeshR = path + "VESSEL12_02_Lungs_R_3.obj";
	std::string existingModelMesh = path + "VESSEL12_02_AIRWAYS_g.obj";

	std::string extractedCenterline = path + "c2.vtk";
	dotObj * extractedCenterlineModel = new dotObj();
	extractedCenterlineModel->initializeFromVTK(extractedCenterline);
	extractedCenterlineModel->exportToFile(path + "cccline");

	std::cout << std::endl << "Model :: Right host mesh " << std::endl;
	dotObj * hostR = new dotObj();
	std::cout << std::endl << "Model :: Right boundary " << std::endl;
	dotObj * boundaryR = new dotObj();
	boundaryR->initializeFromFile(boundaryMeshR);
	std::cout << std::endl << "Model :: Left host mesh " << std::endl;
	dotObj * hostL = new dotObj();
	std::cout << std::endl << "Model :: Left boundary " << std::endl;
	dotObj * boundaryL = new dotObj();
	boundaryL->initializeFromFile(boundaryMeshL);
	std::cout << std::endl << "Model :: Existing geometry " << std::endl;
	dotObj *existingModel = new dotObj();
	existingModel->initializeFromFile(existingModelMesh);
	dotObj * trachea = new dotObj();
	//trachea->initializeFromFile(path + "tracheav8PC.obj");
	if (strlen(path.c_str()) > 2) {
	}
	else {
		std::cout << "Error: Enter valid workspace path" << std::endl;
		return 0;
	}
	simulation *s = new simulation();

	status * st = new status();

	try {
		/*s->extendBronchialTreeV2(
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
			true);

		s->extendBronchialTreeV2(
			buildVolumes,
			true,
			buildCenterline,
			true,
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
			extractedCenterlineModel,
			st,
			true);*/

		/*



	}
	catch (const std::exception& e) { // caught by reference to base
		std::cout << " a standard exception was caught, with message '" << e.what() << "'\n";
	}

	return 1;
}

template<typename TImageType> static void ReadFile(std::string filename, typename TImageType::Pointer image);


template<typename TImageType> void ReadFile(std::string filename, typename TImageType::Pointer image)
{
	typedef itk::ImageFileReader<TImageType> ReaderType;
	typename ReaderType::Pointer reader = ReaderType::New();

	reader->SetFileName(filename);
	reader->Update();

	image->Graft(reader->GetOutput());
}

void flatten(void) {
	dotObj* Reconstructed = new dotObj();
	Reconstructed = new dotObj();
	std::string path = "D:/LungModeller/models/";
	Reconstructed->initializeFromFile(path + "fullReconstructedModel.obj");
	Reconstructed->simplificationEdgeCollapse(1600000);
	Reconstructed->recalculateNormals();
	dotObj* cubef = new dotObj();
	cubef->initializeFromFile(path + "cubef.obj");
	dotObj * sk = new dotObj();
	sk->initializeFromFile(path + "fullTreeExtendedModel.obj");

	Reconstructed->graphBasedAnalysis(*sk);
	Mesh out, modelmesh;
	bool valid_union = true;;
	for (int i = 0; i < Reconstructed->mAnalysis.graph.nodes.size(); i++) {
		if (Reconstructed->mAnalysis.graph.nodes[i].isTerminal) {
			dotObj * newcube = new dotObj();
			*newcube = *cubef;
			Vector3f direction = Reconstructed->mAnalysis.graph.nodes[i].previousNodePtr->position - Reconstructed->mAnalysis.graph.nodes[i].position;
			float d = 2 * Reconstructed->mAnalysis.graph.nodes[i].diameter;
			newcube->scale(d, d, d);
			newcube->rotate(direction);
			direction.normalize();
			Vector3f translatePosition = (Reconstructed->mAnalysis.graph.nodes[i].position - 0.7*d*direction);
			newcube->translate(translatePosition);

			//newcube->exportToFile("oo" + std::to_string(i));
			Mesh mesh;
			std::stringstream input = newcube->toOFF();
			if (!input || !(input >> out))
			{
				std::cerr << "First mesh is not a valid off file." << std::endl;
				return;
			}
		}
		if (Reconstructed->mAnalysis.graph.nodes[i].isInlet) {
			dotObj * newcube = new dotObj();
			*newcube = *cubef;
			Vector3f direction = Reconstructed->mAnalysis.graph.nodes[i].nextNodesPtr[0]->position - Reconstructed->mAnalysis.graph.nodes[i].position;

			float d = 2 * Reconstructed->mAnalysis.graph.nodes[i].nextNodesPtr[0]->diameter;
			newcube->scale(d, d, d);
			newcube->rotate(direction);
			direction.normalize();
			Vector3f translatePosition = (Reconstructed->mAnalysis.graph.nodes[i].position - 0.7*d*direction);
			newcube->translate(translatePosition);
			Mesh mesh;
			std::stringstream input = newcube->toOFF();
			if (!input || !(input >> out))
			{
				std::cerr << "First mesh is not a valid off file." << std::endl;
				return;
			}
		}
	}
	if (valid_union)
	{
		std::cout << "Union was successfully computed\n";
		std::ofstream output("union.off");
		output << out;

		std::stringstream model = Reconstructed->toOFF();
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
			std::ofstream output1("difference.off");
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
		std::ofstream output2("D:/LungModeller/models/difference_remeshed.off");
		output2 << modelmesh;
	}

	return;
}

int segment(void) {
	//dotObj * Reconstructed = new dotObj();
	//Reconstructed->initializeFromFile("D:/LungModeller/models/lungPartsFull.obj");
	//Reconstructed->segmentByGeneration();
	//dotObj * sk = new dotObj();
	//sk->initializeFromFile("../../../../data/skeletonInitialRef.obj");
	//Reconstructed->segmentByGeneration(*sk);
	//Reconstructed->initializeFromFile("D:/LungModeller/models/simplifiedModel.obj");
	//Reconstructed->initializeFromFile("D:/LungModeller/models/fullReconstructedModel.obj");
	//Reconstructed->segmentByGeneration();

	dotObj * Reconstructed = new dotObj();
	dotObj * sk = new dotObj();

	//Reconstructed->initializeFromFile("D:/LungModeller/models/fullReconstructedModel.obj");
	//Reconstructed->simplificationEdgeCollapse(800000);
	//Reconstructed->recalculateNormals();
	//Reconstructed->exportToFile("D:/LungModeller/models/simplifiedModel");
	//Reconstructed = new dotObj();

	std::string type = "OFF";
	Reconstructed->initializeFromFile("D:/LungModeller/models/difference_remeshed.off", type);
	Reconstructed->recalculateNormals();
	//Reconstructed->initializeFromFile("D:/LungModeller/models/difference_remeshed.obj");
	//Reconstructed->smoothingSparse("taubin", "weighted", 4, -0.50, 0.50);
	sk->initializeFromFile("D:/LungModeller/models/fullTreeExtendedModel.obj");
	Reconstructed->segmentByGeneration(*sk);
	Reconstructed->exportToFileSegmentedPerFace("D:/LungModeller/models/final_segmented_model.obj");
	return 0;
}

void main_generate_flatten_segment(void) {
	generate();
	//flatten();
	//segment();
	return;
}

void main_OFF_reader_test(void) {
	dotObj * Reconstructed = new dotObj();
	std::string type = "OFF";
	Reconstructed->initializeFromFile("D:/LungModeller/models/difference_remeshed.off", type);
	Reconstructed->recalculateNormals();
	Reconstructed->exportToFile("testt");
	return;
}

void main_median_narrow(void) {
	dotObj * Reconstructed = new dotObj();
	Reconstructed->initializeFromFile("D:/LungModeller/models/lungPartsFull.obj");
	dotObj skel = Reconstructed->mcfskel();
	skel.exportToFile("skeletonInitial");
	skel.mcfskelRefineStepOne(1);
	skel.mcfskelRefineStepOne(2);
	skel.mcfskelRefineStepOne(3);
	skel.mcfskelRefineStepOne(8);
	skel.mcfskelRefineStepTwo();
	skel.exportToFile("skeletonInitialRef");
	Reconstructed->graphBasedAnalysis(skel);
	Reconstructed->mAnalysis.graph.init;

	std::vector<gnode*> A;
	A.clear();
	gnode* initPoint = &Reconstructed->mAnalysis.graph.nodes[Reconstructed->mAnalysis.inlet];
	A.push_back(initPoint);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextBifurcationPtr;
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				//std::cout << B[i]->nextNodesPtr.size() << std::endl;
				for (int j = 0; j < B[i]->nextNodesPtr.size(); j++) {
					Reconstructed->selectedVertices = B[i]->nextNodesPtr[j]->nextMedianPtr->meshVerticesIndices;
					Reconstructed->extendSelection();
					Reconstructed->extendSelection();
					Reconstructed->extendSelection();
					for (int k = 0; k < Reconstructed->selectedVertices.size(); k++) {
						Reconstructed->vertices[Reconstructed->selectedVertices[k]][0] = 2 * Reconstructed->vertices[Reconstructed->selectedVertices[k]][0];
						Reconstructed->vertices[Reconstructed->selectedVertices[k]][1] = 2 * Reconstructed->vertices[Reconstructed->selectedVertices[k]][1];
						Reconstructed->vertices[Reconstructed->selectedVertices[k]][2] = 2 * Reconstructed->vertices[Reconstructed->selectedVertices[k]][2];
					}
				}
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}

	Reconstructed->exportToFile("20180516");

	return;
}










*/




int main2(void) {
	/*
	dotObj * Reconstructed = new dotObj();
	Reconstructed->initializeFromFile("D:/LungModeller/models/lungPartsFullA.obj");

	dotObj * skel = new dotObj();
	*skel = Reconstructed->mcfskel();
	skel->exportToFile("kokok");
	*/

	//main_generate_flatten_segment();

	//generate();

	return 0;
}

/*

dotObj * Reconstructed = new dotObj();
std::string path = "D:/LungModeller/models/";

dotObj * fullTreeExtendedModel = new dotObj();
fullTreeExtendedModel->initializeFromFile(path + "_7_fullTreeOversampled.obj");

dotObj* cubef = new dotObj();
cubef->initializeFromFile(path + "cubef.obj");

dotObj * mMOdel = new dotObj();
mMOdel->initializeFromFile(path + "_9_simplifiedModel.obj");
mMOdel->graphBasedAnalysis(*fullTreeExtendedModel);
Mesh out, modelmesh;
bool valid_union = true;
for (int i = 0; i < mMOdel->mAnalysis.graph.nodes.size(); i++){
if (mMOdel->mAnalysis.graph.nodes[i].isTerminal){
dotObj * newcube = new dotObj();
*newcube = *cubef;
Vector3f direction = mMOdel->mAnalysis.graph.nodes[i].previousNodePtr->position - mMOdel->mAnalysis.graph.nodes[i].position;
float d = 2 * mMOdel->mAnalysis.graph.nodes[i].diameter;
if (d > 0){
newcube->scale(d, d, d);
newcube->rotate(direction);
direction.normalize();
Vector3f translatePosition = (mMOdel->mAnalysis.graph.nodes[i].position - 0.7*d*direction);
newcube->translate(translatePosition);
Mesh mesh;
std::stringstream input = newcube->toOFF();
if (input){
if (!input || !(input >> out))
{
std::cerr << "First mesh is not a valid off file." << std::endl;
return;
}
}
}
}

if (mMOdel->mAnalysis.graph.nodes[i].isInlet){
dotObj * newcube = new dotObj();
*newcube = *cubef;
Vector3f direction = mMOdel->mAnalysis.graph.nodes[i].nextNodesPtr[0]->position - mMOdel->mAnalysis.graph.nodes[i].position;
float d = 2 * mMOdel->mAnalysis.graph.nodes[i].nextNodesPtr[0]->diameter;
if (d > 0){
newcube->scale(d, d, d);
newcube->rotate(direction);
direction.normalize();
Vector3f translatePosition = (mMOdel->mAnalysis.graph.nodes[i].position - 0.7*d*direction);
newcube->translate(translatePosition);
Mesh mesh;
std::stringstream input = newcube->toOFF();
if (input){
if (!input || !(input >> out))
{
std::cerr << "First mesh is not a valid off file." << std::endl;
return;
}
}
}
}
}

*/





	/*for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			for (int k = 0; k < A[i][j].size(); k++) {
				if (A[i][j][k] == 200) {
					result3->vertices.push_back({ (float)k*dz,(float)j*dy,(float)i*dx });
				}
			}
		}
	}*/

	
	//mGrid grid;
	/*
	grid.boundingBox = rec->boundingBox;
	grid.boundingBox.maxX = float(imageIO->GetDimensions(0) - 1);
	grid.boundingBox.minX = 0.0;
	grid.boundingBox.maxY = float(imageIO->GetDimensions(1) - 1);
	grid.boundingBox.minY = 0.0;
	grid.boundingBox.maxZ = float(imageIO->GetDimensions(2) - 1);
	grid.boundingBox.minZ = 0.0;
	grid.discretization.nCX = imageIO->GetDimensions(0);
	grid.discretization.dx = (grid.boundingBox.maxX - grid.boundingBox.minX) / grid.discretization.nCX;
	grid.discretization.nCY = imageIO->GetDimensions(1);
	grid.discretization.dy = (grid.boundingBox.maxY - grid.boundingBox.minY) / grid.discretization.nCY;
	grid.discretization.nCZ = imageIO->GetDimensions(2);
	grid.discretization.dz = (grid.boundingBox.maxZ - grid.boundingBox.minZ) / grid.discretization.nCZ;

	grid.cells.resize(grid.discretization.nCX);

	for (int i = 0; i < grid.discretization.nCX; i++) {
		grid.cells[i].resize(grid.discretization.nCY);
	}

	for (int i = 0; i < grid.discretization.nCX; i++) {
		for (int j = 0; j < grid.discretization.nCY; j++) {
			grid.cells[i][j].resize(grid.discretization.nCZ);
		}
	}

	for (int i = 0; i < grid.discretization.nCX; i++) {
		for (int j = 0; j < grid.discretization.nCY; j++) {
			for (int k = 0; k < grid.discretization.nCZ; k++) {
				mCell * c;
				c = new mCell();
				float dx = grid.discretization.dx;
				float dy = grid.discretization.dy;
				float dz = grid.discretization.dz;
				float xx = grid.boundingBox.minX + i * dx;
				float yy = grid.boundingBox.minY + j * dy;
				float zz = grid.boundingBox.minZ + k * dz;

				c->position = Vector3f(xx + (dx / 2), yy + (dy / 2), zz + (dz / 2));
				c->boundingBox.minX = xx;
				c->boundingBox.maxX = xx + dx;
				c->boundingBox.minY = yy;
				c->boundingBox.maxY = yy + dy;
				c->boundingBox.minZ = zz;
				c->boundingBox.maxZ = zz + dz;
				c->index.x = i;
				c->index.y = j;
				c->index.z = k;

				grid.cells[i][j][k] = c;
			}
		}
	}
	*/

	/*
	std::cout << "Grid Generation complete" << std::endl;

	for (int i = 0; i < fullTreeExtendedModel->mAnalysis.graph.nodes.size(); i++) {
		Vector3f * p = &fullTreeExtendedModel->mAnalysis.graph.nodes[i].position;

		int ix = grid.discretization.nCX*((p->x() - grid.boundingBox.minX) / (grid.boundingBox.maxX - grid.boundingBox.minX + 0.1));
		int iy = grid.discretization.nCY*((p->y() - grid.boundingBox.minY) / (grid.boundingBox.maxY - grid.boundingBox.minY + 0.1));
		int iz = grid.discretization.nCZ*((p->z() - grid.boundingBox.minZ) / (grid.boundingBox.maxZ - grid.boundingBox.minZ + 0.1));

		grid.cells[ix][iy][iz]->nodes.push_back(&fullTreeExtendedModel->mAnalysis.graph.nodes[i]);
	}

	for (int ix = 0; ix < grid.cells.size(); ix++) {
		for (int iy = 0; iy < grid.cells[ix].size(); iy++) {
			for (int iz = 0; iz < grid.cells[ix][iy].size(); iz++) {
				int * g = &grid.cells[ix][iy][iz]->generation;
				*g = -5;
				gnode * n = new gnode();
				bool doBypass = true;;
				//if ((doBypass) || (LL->checkIfPointisInsidev2(grid.cells[ix][iy][iz]->position)) || (LR->checkIfPointisInsidev2(grid.cells[ix][iy][iz]->position))){
				if (true) {
					*g = -1;
					if (grid.cells[ix][iy][iz]->nodes.size() > 0) {
						n = grid.cells[ix][iy][iz]->nodes[0];
						*g = n->generation;
						for (int im = 0; im < grid.cells[ix][iy][iz]->nodes.size(); im++) {
							n = grid.cells[ix][iy][iz]->nodes[im];
							if (n->generation < *g) {
								*g = n->generation;
							}
						}
					}
				}
			}
		}
	}

	std::cout << "Distribution complete" << std::endl;
	*/



	/*
	for (int ix = 0; ix < grid.cells.size(); ix++){
		for (int iy = 0; iy < grid.cells[ix].size(); iy++){
			for (int iz = 0; iz < grid.cells[ix][iy].size(); iz++){
				mCell * iCell = grid.cells[ix][iy][iz];
				iCell->distanceToGeneration.resize(28);
				for (int z = 0; z < iCell->distanceToGeneration.size(); z++){
					iCell->distanceToGeneration[z].value = 1000000;
					iCell->distanceToGeneration[z].index3D[0] = -1;
					iCell->distanceToGeneration[z].index3D[1] = -1;
					iCell->distanceToGeneration[z].index3D[2] = -1;
				}
			}
		}
	}
	*/

	/*
	for (int ix = 0; ix < grid.cells.size(); ix++){
		for (int iy = 0; iy < grid.cells[ix].size(); iy++){
			for (int iz = 0; iz < grid.cells[ix][iy].size(); iz++){
				int in = grid.cells[ix][iy][iz]->generation;
				mCell * iCell = grid.cells[ix][iy][iz];
				if (in < 0){
					for (int jx = 0; jx < grid.cells.size(); jx++){
						for (int jy = 0; jy < grid.cells[jx].size(); jy++){
							for (int jz = 0; jz < grid.cells[jx][jy].size(); jz++){
								int jn = grid.cells[ix][iy][iz]->generation;
								mCell * jCell = grid.cells[ix][iy][iz];
								if (jn >= 0){
									float currentDistance = (iCell->position - jCell->position).norm();
									if (currentDistance < jCell->distanceToGeneration[jn].value){
										jCell->distanceToGeneration[jn].value = currentDistance;
										jCell->distanceToGeneration[jn].index3D[jx] = currentDistance;
										jCell->distanceToGeneration[jn].index3D[jy] = currentDistance;
										jCell->distanceToGeneration[jn].index3D[jz] = currentDistance;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	*/
	/*
	std::ofstream myfile;
	myfile.open("F:\\_Groundwork\\data_vessel_2_.csv");

	for (int ix = 0; ix < grid.cells.size(); ix++) {
		for (int iy = 0; iy < grid.cells[ix].size(); iy++) {
			for (int iz = 0; iz < grid.cells[ix][iy].size(); iz++) {
				float x = grid.cells[ix][iy][iz]->position.x();
				float y = grid.cells[ix][iy][iz]->position.y();
				float z = grid.cells[ix][iy][iz]->position.z();
				float g = grid.cells[ix][iy][iz]->generation;

				myfile << x << "," << y << "," << z << "," << ix << "," << iy << "," << iz << "," << g << "\n";
			}
		}
	}

	myfile.close();
	*/




	//Extract 1D representation
int main_1d_rep(void)
{
	dotObj * Reconstructed = new dotObj();
	Reconstructed->initializeFromFile("D:/LungModeller/models/lungPartsFullA.obj");
	dotObj * skel = new dotObj();
	*skel = Reconstructed->mcfskel();
	skel->exportToFile("D:/LungModeller/results/representation");
	return 0;
}



//Extend bronchial tree
int main_extend(void)
{
	std::string root = "F:\\_Datasets\\CT-MRI-DATABASE\\_VESSEL\\Scans\\VESSEL12_08\\";
	std::string mhdLungFileName = root + "VESSEL12_08.mhd";
	std::string rawLungFileName = root + "VESSEL12_08.raw";
	std::string mhdLungMaskFileName = root + "VESSEL12_08_Lungs.mhd";
	std::string rawLungMaskFileName = root + "VESSEL12_08_Lungs.raw";

	int seed, label, lim;

	lim = 500;
	std::vector<int> labelsArray = { 100,200 };
	bool buildVolumes = false;
	bool buildCenterline = true;
	bool build1DModel = true;
	bool surfaceSampling = false;
	bool buildNormals = false;
	bool build3DModel = false;
	bool refinements = false;
	bool segmentation = false;
	int depth = 16;
	int volumeDepth = 16;
	int density = 500;
	int poissonDepth = 11;

	voxelSpace * vox = new voxelSpace();
	typedef itk::ImageIOBase::IOComponentType ScalarPixelType;
	itk::ImageIOBase::Pointer imageIO = itk::ImageIOFactory::CreateImageIO(
		mhdLungMaskFileName.c_str(), itk::ImageIOFactory::ReadMode);
	if (!imageIO)
	{
		std::cerr << "Could not CreateImageIO for: " << mhdLungMaskFileName << std::endl;
		return EXIT_FAILURE;
	}
	imageIO->SetFileName(mhdLungMaskFileName.c_str());
	imageIO->ReadImageInformation();
	std::cout << "Pixel Type is " << imageIO->GetComponentType() << std::endl;
	const size_t numDimensions = imageIO->GetNumberOfDimensions();
	std::cout << "Number of Dimensions: " << numDimensions << std::endl;
	std::cout << "Component size: " << imageIO->GetComponentSize() << std::endl;

	std::cout << "Pixel type: " << imageIO->GetPixelTypeAsString(imageIO->GetPixelType()) << std::endl;

	std::cout << "Byte Order: " << imageIO->GetByteOrder() << std::endl;
	std::cout << "Dimension  " << imageIO->GetDimensions(0) << std::endl;
	std::cout << "Dimension  " << imageIO->GetDimensions(1) << std::endl;
	std::cout << "Dimension  " << imageIO->GetDimensions(2) << std::endl;

	std::cout << "Component type: " << imageIO->GetComponentTypeAsString(imageIO->GetComponentType()) << std::endl;

	std::cout << "Byte Order: " << imageIO->GetByteOrderAsString(imageIO->GetByteOrder()) << std::endl;

	float dx = imageIO->GetSpacing(0);
	float dy = imageIO->GetSpacing(1);
	float dz = imageIO->GetSpacing(2);
	int LX = imageIO->GetDimensions(0);
	int LY = imageIO->GetDimensions(1);
	int LZ = imageIO->GetDimensions(2);

	std::vector<std::vector<std::vector<int>>> A;
	vox->resize3DMat(A, mDiscretization(LZ, LY, LX));

	VectorXf  a2 = VectorXf::Zero(imageIO->GetDimensions(0)*imageIO->GetDimensions(1)*imageIO->GetDimensions(2));
	std::ifstream myData(rawLungMaskFileName, std::ios::binary);
	unsigned char value;
	int i = 0;
	char buf[sizeof(unsigned char)];
	while (myData.read(buf, sizeof(buf)))
	{
		memcpy(&value, buf, sizeof(value));
		//std::cout << value << " ";
		a2(i) = value;
		i++;
	}
	std::cout << std::endl << "Total count: " << i << std::endl;
	int idx = 0;
	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			for (int k = 0; k < A[i][j].size(); k++) {
				int val = a2(idx++);
				A[i][j][k] = val;
			}
		}
	}

	std::vector<mCellSimplified> Asparse;
	vox->Matrix2SparseRepresentationByValue(A, Asparse, 1.0);

	label = labelsArray[0];
	std::vector<mIndex> toInvestigate;
	seed = vox->findSeed(A, Asparse, 1);
	mIndex seedIndex = Asparse[seed].index;
	A[seedIndex.x][seedIndex.y][seedIndex.z] = label;

	toInvestigate.clear();
	toInvestigate.push_back(seedIndex);
	int regionGrowingFront = 0;
	for (int pp = 0; pp < toInvestigate.size(); pp++) {
		mIndex * ind = &toInvestigate[pp];
		if (ind->isInsideMatrix(A, 2)) {
			bool found = vox->regionHasNotLabel(ind, A, 2, 0);
			bool isNeighb = vox->regionHasOnlyLabel(ind, A, 1, label);
			if (found && isNeighb) {
				A[ind->x][ind->y][ind->z] = label;
			}
		}
		if (pp == toInvestigate.size() - 1) {
			regionGrowingFront++;
			int c = regionGrowingFront;
			int r1, r2, r3;
			toInvestigate.clear();
			pp = 0;
			for (r1 = seedIndex.x - c; r1 <= seedIndex.x + c; r1 = r1 + 2 * c) {
				for (r2 = seedIndex.y - c; r2 <= seedIndex.y + c; r2++) {
					for (r3 = seedIndex.z - c; r3 <= seedIndex.z + c; r3++) {
						toInvestigate.push_back(mIndex(r1, r2, r3));
					}
				}
			}
			for (r2 = seedIndex.y - c; r2 <= seedIndex.y + c; r2 = r2 + 2 * c) {
				for (r1 = seedIndex.x - c; r1 <= seedIndex.x + c; r1++) {
					for (r3 = seedIndex.z - c; r3 <= seedIndex.z + c; r3++) {
						toInvestigate.push_back(mIndex(r1, r2, r3));
					}
				}
			}
			for (r3 = seedIndex.z - c; r3 <= seedIndex.z + c; r3 = r3 + 2 * c) {
				for (r1 = seedIndex.x - c; r1 <= seedIndex.x + c; r1++) {
					for (r2 = seedIndex.y - c; r2 <= seedIndex.y + c; r2++) {
						toInvestigate.push_back(mIndex(r1, r2, r3));
					}
				}
			}
		}
		if (regionGrowingFront == lim) {
			regionGrowingFront = 0;
			toInvestigate.clear();
		}
	}

	label = labelsArray[1];

	seed = vox->findSeed(A, Asparse, 1);
	seedIndex = Asparse[seed].index;
	A[seedIndex.x][seedIndex.y][seedIndex.z] = label;
	toInvestigate.clear();
	toInvestigate.push_back(seedIndex);
	regionGrowingFront = 0;

	for (int pp = 0; pp < toInvestigate.size(); pp++) {
		mIndex * ind = &toInvestigate[pp];
		if (ind->isInsideMatrix(A, 2)) {
			if (vox->regionHasNotLabel(ind, A, 2, 0) && vox->regionHasNotLabel(ind, A, 2, label) && vox->regionHasOnlyLabel(ind, A, 1, label)) {
				A[ind->x][ind->y][ind->z] = label;
			}
		}
		if (pp == toInvestigate.size() - 1) {
			regionGrowingFront++;
			int c = regionGrowingFront;
			int r1, r2, r3;
			toInvestigate.clear();
			pp = 0;
			for (r1 = seedIndex.x - c; r1 <= seedIndex.x + c; r1 = r1 + 2 * c) {
				for (r2 = seedIndex.y - c; r2 <= seedIndex.y + c; r2++) {
					for (r3 = seedIndex.z - c; r3 <= seedIndex.z + c; r3++) {
						toInvestigate.push_back(mIndex(r1, r2, r3));
					}
				}
			}
			for (r2 = seedIndex.y - c; r2 <= seedIndex.y + c; r2 = r2 + 2 * c) {
				for (r1 = seedIndex.x - c; r1 <= seedIndex.x + c; r1++) {
					for (r3 = seedIndex.z - c; r3 <= seedIndex.z + c; r3++) {
						toInvestigate.push_back(mIndex(r1, r2, r3));
					}
				}
			}
			for (r3 = seedIndex.z - c; r3 <= seedIndex.z + c; r3 = r3 + 2 * c) {
				for (r1 = seedIndex.x - c; r1 <= seedIndex.x + c; r1++) {
					for (r2 = seedIndex.y - c; r2 <= seedIndex.y + c; r2++) {
						toInvestigate.push_back(mIndex(r1, r2, r3));
					}
				}
			}
		}
		if (regionGrowingFront == lim) {
			regionGrowingFront = 0;
			toInvestigate.clear();
		}
	}

	for (int i = 0; i < A.size(); i++) {
		for (int j = 0; j < A[i].size(); j++) {
			for (int k = 0; k < A[i][j].size(); k++) {
				if (A[i][j][k] == 1) {
					A[i][j][k] = 0;
				}
			}
		}
	}

	/*
	int x,y,z;
	bbox limits;
	limits.maxX = LX-1;
	limits.minX = 0;
	limits.maxY = LY-1;
	limits.minY = 0;
	limits.maxZ = LZ-1;
	limits.minZ = 0;
	*/

	volume * result2 = new volume();
	volume * result3 = new volume();

	//result2->initializeFromFile(path + "LL.obj");
	//result3->initializeFromFile(path + "RL.obj");

	int maxPoints = 30000;
	helper *fns = new helper();

	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A[i].size(); j++)
		{
			for (int k = 0; k < A[i][j].size(); k++)
			{
				int ii = (int)((A.size() - 1)*fns->uniformRandom());
				int jj = (int)((A[i].size() - 1)*fns->uniformRandom());
				int kk = (int)((A[i][j].size() - 1)*fns->uniformRandom());
				if (result2->vertices.size() <= maxPoints) {
					if (A[ii][jj][kk] == labelsArray[0]) {
						result2->vertices.push_back({ (float)kk*dz,(float)jj*dy,(float)ii*dx });
					}
				}
			}
		}
	}

	for (int i = 0; i < A.size(); i++)
	{
		for (int j = 0; j < A[i].size(); j++)
		{
			for (int k = 0; k < A[i][j].size(); k++)
			{
				int ii = (int)((A.size() - 1)*fns->uniformRandom());
				int jj = (int)((A[i].size() - 1)*fns->uniformRandom());
				int kk = (int)((A[i][j].size() - 1)*fns->uniformRandom());
				if (result3->vertices.size() <= maxPoints) {
					if (A[ii][jj][kk] == labelsArray[1]) {
						result3->vertices.push_back({ (float)kk*dz,(float)jj*dy,(float)ii*dx });
					}
				}
			}
		}
	}

	result2->exportToFile(root + "LL.obj");
	result3->exportToFile(root + "RL.obj");

	std::string path = root;
	std::string extractedCenterline = path + "VESSEL12_08_Centerline.vtk";
	std::string existingModelMesh = path + "existingModel.obj";
	dotObj * extractedCenterlineModel = new dotObj();
	extractedCenterlineModel->initializeFromVTK(extractedCenterline);
	std::cout << std::endl << "Model :: Right host mesh " << std::endl;
	dotObj * hostR = new dotObj();
	std::cout << std::endl << "Model :: Right boundary " << std::endl;
	dotObj * boundaryR = new dotObj();
	//boundaryR->initializeFromFile(boundaryMeshR);
	std::cout << std::endl << "Model :: Left host mesh " << std::endl;
	dotObj * hostL = new dotObj();
	std::cout << std::endl << "Model :: Left boundary " << std::endl;
	dotObj * boundaryL = new dotObj();
	//boundaryL->initializeFromFile(boundaryMeshL);
	std::cout << std::endl << "Model :: Existing geometry " << std::endl;
	dotObj *existingModel = new dotObj();
	existingModel->initializeFromFile(existingModelMesh);
	dotObj * trachea = new dotObj();

	if (strlen(path.c_str()) < 2) {
		std::cout << "Error: Enter valid workspace path" << std::endl;
		return 0;
	}

	simulation *s = new simulation();
	status * st = new status();

	try {
		/*s->extendBronchialTreeV2(
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
			true);*/

		s->extendBronchialTreeV2(
			buildVolumes,
			buildCenterline,
			true,
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
			extractedCenterlineModel,
			st,
			true,
			result2,
			result3
		);
	}
	catch (const std::exception& e) {
		std::cout << " a standard exception was caught, with message '" << e.what() << "'\n";
	}

	delete s;
}

