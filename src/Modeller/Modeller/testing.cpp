#include "lungModelling.h"

/*
void ggraph::ggraph2ObjB(dotObj * result){
//result = new dotObj();
result->vertices.push_back({
this->nodes[this->init->index].position.x(),
this->nodes[this->init->index].position.y(),
this->nodes[this->init->index].position.z()
});
std::vector<gnode*> A;
//A.push_back(&this->nodes[this->init->index]);
A.push_back(this->init);
//this->nodes[this->init->index].index = result.vertices.size() - 1;
this->init->index = result->vertices.size() - 1;
while (!A.empty()){
std::vector<gnode*> B = A.back()->nextBifurcationPtr;
A.pop_back();
for (int i = 0; i < B.size(); i++){
gnode * n = B[i];
result->vertices.push_back({ B[i]->position.x(), B[i]->position.y(), B[i]->position.z() });
B[i]->index = result->vertices.size() - 1;
result->lines.push_back({ B[i]->previousBifurcationPtr->index + 1, B[i]->index + 1 });
}
//
//A.erase(A.begin());
if (!B.empty()){
A.insert(A.end(), B.begin(), B.end());
}
}
return;
}
*/

int main_1(void){
	try {
		std::cout << "Throwing an integer exception...\n";
		throw 42;
	}
	catch (int i) {
		std::cout << " the integer exception was caught, with value: " << i << '\n';
	}

	try {
		std::cout << "Creating a vector of size 5... \n";
		std::vector<int> v(5);
		std::cout << "Accessing the 11th element of the vector...\n";
		std::cout << v.at(10); // vector::at() throws std::out_of_range
	}
	catch (const std::exception& e) { // caught by reference to base
		std::cout << " a standard exception was caught, with message '"
			<< e.what() << "'\n";
	}
	return 0;
}

int main_4(void){
	dotObj * Reconstructed = new dotObj();
	Reconstructed->initializeFromFile("D:/LungModeller/models/lungPartsFull.obj");
	Reconstructed->smoothingSparse("taubin", "weighted", 300, -0.7, 0.5);

	Reconstructed->exportToFile("SmoothedTau");

	return 0;
}

int main_FUCKED(void){
	dotObj * fullTreeModel = new dotObj();
	fullTreeModel->initializeFromFile("D:/LungModeller/models/fulltree.obj");
	fullTreeModel->splitSkeletonEdges(0.2);
	fullTreeModel->mAnalysis.generateGraph(*fullTreeModel, 0);
	ggraph fullTreeExtendedGraph = fullTreeModel->mAnalysis.graph;
	dotObj * fullTreeExtendedModel = new dotObj();
	fullTreeExtendedGraph.ggraph2Model(fullTreeExtendedModel);
	fullTreeExtendedModel->exportToFile("D:/LungModeller/models/fullTreeExtendedModel");

	return 0;
}




/*void graphBasedFlattening(void){
Mesh mesh1, mesh2;

std::stringstream input1;// = mod->toOFF();
std::stringstream input2;// = cube->toOFF();

if (!input1 || !(input1 >> mesh1))
{
std::cerr << "First mesh is not a valid off file." << std::endl;
return;
}

if (!input2 || !(input2 >> mesh2))
{
std::cerr << "Second mesh is not a valid off file." << std::endl;
return;
}

//create a property on edges to indicate whether they are constrained
Mesh::Property_map<edge_descriptor, bool> is_constrained_map =
mesh1.add_property_map<edge_descriptor, bool>("e:is_constrained", false).first;
// update mesh1 to contain the mesh bounding the difference
// of the two input volumes.
bool valid_difference =
PMP::corefine_and_compute_difference(mesh1,
mesh2,
mesh1,
params::all_default(), // default parameters for mesh1
params::all_default(), // default parameters for mesh2
params::edge_is_constrained_map(is_constrained_map));
if (valid_difference)
{
std::cout << "Difference was successfully computed\n";
std::ofstream output("difference.off");
output << mesh1;
}
else{
std::cout << "Difference could not be computed\n";
return;
}

// collect faces incident to a constrained edge
std::vector<face_descriptor> selected_faces;
std::vector<bool> is_selected(num_faces(mesh1), false);
BOOST_FOREACH(edge_descriptor e, edges(mesh1)){
if (is_constrained_map[e])
{
// insert all faces incident to the target vertex
BOOST_FOREACH(halfedge_descriptor h, halfedges_around_target(halfedge(e, mesh1), mesh1))
{
if (!is_border(h, mesh1))
{
face_descriptor f = face(h, mesh1);
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
CGAL::expand_face_selection(selected_faces, mesh1, 2, Vector_pmap_wrapper(is_selected), std::back_inserter(selected_faces));
std::cout << selected_faces.size()
<< " faces were selected for the remeshing step\n";
// remesh the region around the intersection polylines
PMP::isotropic_remeshing(
selected_faces,
0.02,
mesh1,
params::edge_is_constrained_map(is_constrained_map));
std::ofstream output("difference_remeshed.off");
output << mesh1;
return;
}*/