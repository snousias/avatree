#include "mainwindow.h"
#include "glwidget.h"
#include "modellerToViewer.h"
#include <QMenuBar>
#include <QMenu>
#include <QMessageBox>
#include <QFileDialog>
#include <QDialog>
#include <QTimer>
#include <QDebug>

namespace PMP = CGAL::Polygon_mesh_processing;



MainWindow::MainWindow() : m_timeWindowframe(0),m_lastFrame(10000)
{
	setupUi(this);
	connect(sdfButton, SIGNAL(pressed()), this, SLOT(sdf()));
	connect(reloadButton, SIGNAL(pressed()), this, SLOT(reload()));
	connect(exposrtSelectedVerticesToFile, SIGNAL(pressed()), this, SLOT(exportVerticesToText()));
	connect(writeButton, SIGNAL(pressed()), this, SLOT(writeObj()));
	connect(narrowingButton, SIGNAL(pressed()), this, SLOT(narrow()));
	connect(exportButton, SIGNAL(pressed()), this, SLOT(exportObj()));
	connect(exportObjectsAsOBJGroups, SIGNAL(pressed()), this, SLOT(exportSegmented()));
	connect(existingGeomPathButton, SIGNAL(pressed()), this, SLOT(existingGeomSetPath()));
	connect(leftLungGeomPathButton, SIGNAL(pressed()), this, SLOT(leftLungGeomSetPath()));
	connect(rightLungGeomPathButton, SIGNAL(pressed()), this, SLOT(rightLungGeomSetPath()));
	connect(workspacePathButton, SIGNAL(pressed()), this, SLOT(workspaceFolderPathToLineEdit()));
	connect(useBrush, SIGNAL(clicked(bool)), this, SLOT(updateBrushType()));
	connect(brushSizeBox, SIGNAL(valueChanged(int)), this, SLOT(updateBrushSize()));
	connect(brushDistance, SIGNAL(valueChanged(double)), this, SLOT(updateBrushDistance()));
	connect(distanceMode, SIGNAL(currentIndexChanged(const QString&)), this, SLOT(updateDistanceMode()));
	connect(loadModelButton, SIGNAL(pressed()), this, SLOT(loadFile()));
	connect(loadSegmentPropertyListButton, SIGNAL(pressed()), this, SLOT(loadSegmentPropertyList()));
	connect(extendSelectionButton, SIGNAL(pressed()), this, SLOT(extendSelection()));
	connect(selectWorkspaceButton, SIGNAL(pressed()), this, SLOT(initializeWorkspace()));
	connect(cfdSim, SIGNAL(pressed()), this, SLOT(doSimulate()));
	connect(segmentByGenerationButton, SIGNAL(pressed()), this, SLOT(segmentByGeneration()));
	connect(segmentBySkeletonButton, SIGNAL(pressed()), this, SLOT(surfaceMeshSkeletonBasedSegmentation()));
	//connect(extendLungModel, SIGNAL(pressed()), this, SLOT(extendGeometryV3()));










	QString myAirCoachLogoPath = ":/Images/mAClogo.png";
	QImage myAirCoachLogoImage(myAirCoachLogoPath);

	this->myAirCoachLogo->setPixmap(QPixmap::fromImage(myAirCoachLogoImage).scaledToHeight(80));
	this->myAirCoachLogo->setAlignment(Qt::AlignCenter | Qt::AlignBaseline);

	QString H2020LogoPath = ":/Images/H2020logo.png";
	QImage H2020LogoImage(H2020LogoPath);

	this->H2020Logo->setPixmap(QPixmap::fromImage(H2020LogoImage.scaledToHeight(90)));
	this->H2020Logo->setAlignment(Qt::AlignCenter | Qt::AlignBaseline);
	//Check if model exists to activate buttons
	if (lungmodel.size() > 0){
		if ((lungmodel.at(lungmodel.size() - 1).vertices.size() > 0) && (lungmodel.at(lungmodel.size() - 1).normals.size() > 0) && (lungmodel.at(lungmodel.size() - 1).faces.size() > 0)){
			sdfButton->setEnabled(true);
			narrowingButton->setEnabled(true);
		}
		else{
			sdfButton->setEnabled(false);
			narrowingButton->setEnabled(false);
		}
	}
	else{
		sdfButton->setEnabled(false);
		narrowingButton->setEnabled(false);
	}
}

MainWindow::~MainWindow()
{}

//GUI RELATED FUNCTIONALITIES

void MainWindow::reload(){
	openGLWidget->loadNewOBJ(&lungmodel.at(lungmodel.size() - 1));
	openGLWidget->m_camera.calibrate(openGLWidget->m_model->bbox());

	//Check if model exists to activate buttons
	if (lungmodel.size() > 0){
		if ((lungmodel.at(lungmodel.size() - 1).vertices.size() > 0) && (lungmodel.at(lungmodel.size() - 1).normals.size() > 0) && (lungmodel.at(lungmodel.size() - 1).faces.size() > 0)){
			sdfButton->setEnabled(true);
			narrowingButton->setEnabled(true);
		}
		else{
			sdfButton->setEnabled(false);
			narrowingButton->setEnabled(false);
		}
	}
	else{
		sdfButton->setEnabled(false);
		narrowingButton->setEnabled(false);
	}

	return;
}

void MainWindow::extendSelection(){
	if (lungmodel.at(lungmodel.size() - 1).selectedVertices.size() > 0){
		lungmodel.at(lungmodel.size() - 1).extendSelection();
	}

	this->openGLWidget->m_model->populateColorsPerVertexList(lungmodel.at(lungmodel.size() - 1).selectedVertices);

	return;
}

void MainWindow::updateBrushType(void)
{
	brushSelectionFunctionality = useBrush->isChecked();
	partSelectionFunctionality = !brushSelectionFunctionality;
	qDebug() << "Brush type changed : " << brushSelectionFunctionality;
	return;
}

void MainWindow::updateBrushSize(void)
{
	_brushSizeBox = brushSizeBox->value();

	qDebug() << "Brush size changed : " << _brushSizeBox;
	return;
}

void MainWindow::updateDistanceMode(void)
{
	_brushdistanceMode = distanceMode->currentText().toStdString();
	_brushdistanceModeIndex = distanceMode->currentIndex();
	qDebug() << "Distance Mode Changed : " << QString::fromStdString(_brushdistanceMode);
	qDebug() << "Distance Mode Index Changed : " << _brushdistanceModeIndex;
	return;
}

void MainWindow::updateBrushDistance(void)
{
	_brushDistance = brushDistance->value();

	qDebug() << "Brush size changed : " << _brushDistance;
	return;
}

void MainWindow::writeObj(void)
{
	if (lungmodel.at(lungmodel.size() - 1).vertices.size() > 0){
		if (lungmodel.size() > 2){
			lungmodel.erase(lungmodel.begin());
		}
		lungmodel.push_back(*tempModel);
		lungmodel.at(lungmodel.size() - 1).exportToFile("update");
	}

	return;
}

void MainWindow::existingGeomSetPath(void){
	QString dir = QString("../../../../models/");
	QString filename = QFileDialog::getOpenFileName(this, tr("Select surface file"), dir, tr("Wavefront obj (*.obj)"));
	std::string path = filename.toStdString();
	existingGeomPath->setText(QString::fromStdString(path));
	return;
}

void MainWindow::leftLungGeomSetPath(void){
	QString dir = QString("../../../../models/");
	QString filename = QFileDialog::getOpenFileName(this, tr("Select surface file"), dir, tr("Wavefront obj (*.obj)"));
	std::string path = filename.toStdString();
	leftLungGeomPath->setText(QString::fromStdString(path));
	return;
}

void MainWindow::rightLungGeomSetPath(void){
	QString dir = QString("../../../../models/");
	QString filename = QFileDialog::getOpenFileName(this, tr("Select surface file"), dir, tr("Wavefront obj (*.obj)"));
	std::string path = filename.toStdString();
	rightLungGeomPath->setText(QString::fromStdString(path));
	return;
}

void MainWindow::workspaceFolderPathToLineEdit(void){
	QString dir = QString("../../../../models/");
	QString directory = QFileDialog::getExistingDirectory(this, tr("Select Working Directory"), dir, QFileDialog::ShowDirsOnly
		| QFileDialog::DontResolveSymlinks);
	std::string path = directory.toStdString();
	workspacePath->setText(QString::fromStdString(path));
	return;
}

void MainWindow::loadFile(void){
	QString dir = QString("../../../../models/");
	QString filename = QFileDialog::getOpenFileName(this, tr("Select surface file"), dir, tr("Wavefront obj (*.obj)"));
	std::string modelstring = filename.toStdString();
	dotObj buffer;
	buffer.initializeFromFile(modelstring);

	//buffer.segment_property_map = lungmodel.at(lungmodel.size() - 1).segment_property_map;
	//buffer.segment_property_map_per_vertex = lungmodel.at(lungmodel.size() - 1).segment_property_map_per_vertex;

	lungmodel.at(lungmodel.size() - 1) = buffer;

	//Check if model exists to activate buttons
	if (lungmodel.size() > 0){
		if ((lungmodel.at(lungmodel.size() - 1).vertices.size() > 0) && (lungmodel.at(lungmodel.size() - 1).normals.size() > 0) && (lungmodel.at(lungmodel.size() - 1).faces.size() > 0)){
			sdfButton->setEnabled(true);
			narrowingButton->setEnabled(true);
		}
		else{
			sdfButton->setEnabled(false);
			narrowingButton->setEnabled(false);
		}
	}
	else{
		sdfButton->setEnabled(false);
		narrowingButton->setEnabled(false);
	}

	openGLWidget->loadNewOBJ(&lungmodel.at(lungmodel.size() - 1));
	openGLWidget->m_camera.calibrate(openGLWidget->m_model->bbox());
	return;
}

//SIMULATION RELATED FUNCTIONALITIES

void MainWindow::exportSegmented(){
	QString format = "obj";
	QString initialPath = QDir::currentPath() + tr("/untitled.") + format;
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Model"), initialPath,
		tr("%1 Files (*.%2);;All Files (*)").arg(format.toUpper()).arg(format));
	if (!fileName.isEmpty()){
		lungmodel.at(lungmodel.size() - 1).exportToFileSegmented(fileName.toStdString());
	}

	return;
}

void MainWindow::threadtester(){
	return;
}

void MainWindow::narrow()
{
	std::cout << "Narrowing commencing" << std::endl;

	if (lungmodel.at(lungmodel.size() - 1).selectedVertices.size() > 0){
		simulation *sim;
		sim = new simulation();
		sim->lungmodel = lungmodel;
		sim->narrow(contractionPercentageBox->value(), contractionStrengthBox->value(), useCuFnCheckBox->isChecked(), frequencyBox->value(), useInterpolation->isChecked(), interpolationPercentage->value(), useExtendedSmoothing->isChecked(), numberOfRaysSDF->value(), numberOfClustersSDF->value(), lamdaSDF->value());
		narrowingPercentage->setValue((int)(sim->narrowingRatio * 100));
		tempModel = sim->tempModel;
		//Reload new model
		openGLWidget->loadNewOBJ(tempModel);
		openGLWidget->m_camera.calibrate(openGLWidget->m_model->bbox());
	}
	else{
		std::cout << "\n";
		std::cout << "Error no area selected" << "\n";

		msgBox.setText("Error no area selected");
		msgBox.exec();
	}

	return;
		
}

void MainWindow::exportObj(void){
	QString format = "obj";
	QString initialPath = QDir::currentPath() + tr("/untitled.") + format;
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Model"), initialPath,
		tr("%1 Files (*.%2);;All Files (*)").arg(format.toUpper()).arg(format));
	if (!fileName.isEmpty())
		openGLWidget->saveOBJ(fileName);
	// Our camera changes in this example.
	openGLWidget->m_camera.calibrate(openGLWidget->m_model->bbox());
}

void MainWindow::exportVerticesToText(){
	QString format = "txt";
	QString initialPath = QDir::currentPath() + tr("/untitled.") + format;
	QString fileName = QFileDialog::getSaveFileName(this, tr("Export vertices"), initialPath,
		tr("%1 Files (*.%2);;All Files (*)").arg(format.toUpper()).arg(format));

	std::string newfile = fileName.toStdString();
	lungmodel.at(lungmodel.size() - 1).exportSelectedIndicesToFile(newfile);

	return;
}

void MainWindow::loadSegmentPropertyList(){
	QString dir = QString("../../../../models/");
	QString filename = QFileDialog::getOpenFileName(this, tr("Select selections file"), dir, tr("CSV file (*.csv)"));
	std::string fname = filename.toStdString();

	std::ifstream myfile(fname);
	if (!myfile.fail()){
		float output;
		std::string line;
		std::vector<int> v;

		if (myfile.is_open())
		{
			while (!myfile.eof())
			{
				getline(myfile, line);

				output = (float)atof(line.c_str());
				if (output > 0){ v.push_back(output); }
			}
		}

		std::cout << "File list size:" << v.size() << std::endl;
		std::cout << "Faces size:" << lungmodel.at(lungmodel.size() - 1).faces.size() << std::endl;

		if (v.size() == lungmodel.at(lungmodel.size() - 1).faces.size()){
			std::vector<int>::iterator resultMax, resultMin;
			resultMax = std::max_element(v.begin(), v.end());
			resultMin = std::min_element(v.begin(), v.end());
			int max = v.at(std::distance(v.begin(), resultMax));
			int min = v.at(std::distance(v.begin(), resultMin));
			std::vector<int> list;
			std::cout << "Min is  " << min << " Max is  " << max << std::endl;
			for (int j = 0; j < max; j++){
				if (list.size() > 0){ list.clear(); }
				for (int i = 0; i < v.size(); i++){
					if (v.at(i) == j){
						list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(0) - 1);
						list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(3) - 1);
						list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(6) - 1);
					}
				}
				openGLWidget->m_model->populateColorsPerVertexList(list, j*(1.0 / (double)max));
			}
		}
	}
	return;
}

void MainWindow::sdf(){
	simulation *sim;
	sim = new simulation();

	if (lungmodel.at(lungmodel.size() - 1).vertices.size() > 0){
		if (lungmodel.at(lungmodel.size() - 1).normals.size() > 0){
			if (lungmodel.at(lungmodel.size() - 1).faces.size() > 0){
				sim->lungmodel = lungmodel;
				QString des = segmentOrCluster->currentText();
				sim->sdf(des.toStdString(), coneAngleSDF->value(), numberOfRaysSDF->value(), postProcessingSDF->isChecked(), numberOfClustersSDF->value(), lamdaSDF->value());

				/*
				if (lungmodel.at(lungmodel.size() - 1).vertices.size() > 0){
				bool exportClusters = false;
				if (segmentOrCluster->currentText() == "Per cluster"){ exportClusters = true; }
				lungmodel.at(lungmodel.size() - 1).getSDF(coneAngleSDF->value(), numberOfRaysSDF->value(), postProcessingSDF->isChecked(), true, numberOfClustersSDF->value(), lamdaSDF->value(), exportClusters);  //getSDF(0.005,25,true,4,0.001);
				*/

				lungmodel = sim->lungmodel;

				std::vector<int> v = lungmodel.at(lungmodel.size() - 1).segment_property_map;

				openGLWidget->m_model->populateColorsPerVertexWithSegmentProperty(v);

				/*std::vector<int>::iterator resultMax, resultMin;
				resultMax = std::max_element(v.begin(), v.end());
				resultMin = std::min_element(v.begin(), v.end());
				int max = v.at(std::distance(v.begin(), resultMax));
				int min = v.at(std::distance(v.begin(), resultMin));
				std::vector<int> list;
				std::cout << "Min is  " << min << " Max is  " << max << std::endl;
				for (int j = 0; j < max; j++){
					if (list.size() > 0){ list.clear(); }
					for (int i = 0; i < v.size(); i++){
						if (v.at(i) == j){
							list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(0) - 1);
							list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(3) - 1);
							list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(6) - 1);
						}
					}
					openGLWidget->m_model->populateColorsPerVertexList(list, j*(1.0 / (double)max));
				}*/

			}
		}
	}
	return;
}





void MainWindow::initializeWorkspace(void){
	QString dir = QString("../../../../models/");
	QString directory = QFileDialog::getExistingDirectory(this, tr("Select Working Directory"), dir, QFileDialog::ShowDirsOnly
		| QFileDialog::DontResolveSymlinks);
	std::string path = directory.toStdString();
	workspacePath->setText(QString::fromStdString(path));
	leftLungGeomPath->setText(QString::fromStdString(path + "/L1.obj"));
	rightLungGeomPath->setText(QString::fromStdString(path + "/L2.obj"));
	existingGeomPath->setText(QString::fromStdString(path + "/lungPartsFull.obj"));
	
	return;
}

void MainWindow::doSimulate(void){
#ifdef _WIN64
	system("InitSim.exe");
#elif _WIN32
	system("InitSim.exe");
#elif __unix
	system("./InitSim");
#endif

	return;
}



void MainWindow::segmentByGeneration(void){
	//simulation *sim;
	//sim = new simulation();

	if (lungmodel.at(lungmodel.size() - 1).vertices.size() > 0){
		if (lungmodel.at(lungmodel.size() - 1).normals.size() > 0){
			if (lungmodel.at(lungmodel.size() - 1).faces.size() > 0){
				lungmodel.at(lungmodel.size() - 1).segmentByGeneration();

				std::vector<int> v = lungmodel.at(lungmodel.size() - 1).segment_property_map;
				std::vector<int>::iterator resultMax, resultMin;
				resultMax = std::max_element(v.begin(), v.end());
				resultMin = std::min_element(v.begin(), v.end());
				int max = v.at(std::distance(v.begin(), resultMax));
				int min = v.at(std::distance(v.begin(), resultMin));
				std::vector<int> list;
				std::cout << "Min is  " << min << " Max is  " << max << std::endl;

				for (int j = 0; j < max; j++){
					if (list.size() > 0){ list.clear(); }
					for (int i = 0; i < v.size(); i++){
						if (v.at(i) == j){
							list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(0) - 1);
							list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(3) - 1);
							list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(6) - 1);
						}
					}
					openGLWidget->m_model->populateColorsPerVertexList(list, j*(1.0 / (double)max));
				}
			}
		}
	}




	return;
}

void MainWindow::segmentBySkeleton(void){
	simulation *sim;
	sim = new simulation();

	if (lungmodel.at(lungmodel.size() - 1).vertices.size() > 0){
		if (lungmodel.at(lungmodel.size() - 1).normals.size() > 0){
			if (lungmodel.at(lungmodel.size() - 1).faces.size() > 0){
				lungmodel.at(lungmodel.size() - 1).segmentBySkeleton();

				std::vector<int> v = lungmodel.at(lungmodel.size() - 1).segment_property_map;
				std::vector<int>::iterator resultMax, resultMin;
				resultMax = std::max_element(v.begin(), v.end());
				resultMin = std::min_element(v.begin(), v.end());
				int max = v.at(std::distance(v.begin(), resultMax));
				int min = v.at(std::distance(v.begin(), resultMin));
				std::vector<int> list;
				std::cout << "Min is  " << min << " Max is  " << max << std::endl;

				for (int j = 0; j < max; j++){
					if (list.size() > 0){ list.clear(); }
					for (int i = 0; i < v.size(); i++){
						if (v.at(i) == j){
							list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(0) - 1);
							list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(3) - 1);
							list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(6) - 1);
						}
					}
					openGLWidget->m_model->populateColorsPerVertexList(list, j*(1.0 / (double)max));
				}
			}
		}
	}

	return;
}




/*


void MainWindow::shapeAnalysis(void){
simulation *sim;
sim = new simulation();

if (lungmodel.at(lungmodel.size() - 1).vertices.size() > 0){
if (lungmodel.at(lungmodel.size() - 1).normals.size() > 0){
if (lungmodel.at(lungmodel.size() - 1).faces.size() > 0){
dotObj skel = lungmodel.at(lungmodel.size() - 1).mcfskel();
skel.exportToFile("skeletonInitial");
skel.mcfskelRefineStepOne();
skel.exportToFile("skeletonInitialRef");
skel.mcfskelRefine();
skel.exportToFile("skeletonInitialRefined");

//lungmodel.at(lungmodel.size() - 1).surf_mesh_skel_based_segm();

std::vector<int> v = lungmodel.at(lungmodel.size() - 1).segment_property_map;
std::vector<int>::iterator resultMax, resultMin;
resultMax = std::max_element(v.begin(), v.end());
resultMin = std::min_element(v.begin(), v.end());
int max = v.at(std::distance(v.begin(), resultMax));
int min = v.at(std::distance(v.begin(), resultMin));
std::vector<int> list;
std::cout << "Min is  " << min << " Max is  " << max << std::endl;

for (int j = 0; j < max; j++){
if (list.size() > 0){ list.clear(); }
for (int i = 0; i < v.size(); i++){
if (v.at(i) == j){
list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(0) - 1);
list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(3) - 1);
list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(6) - 1);
}
}
openGLWidget->m_model->populateColorsPerVertexList(list, j*(1.0 / (double)max));
}
}
}
}

return;
}


*/


struct halfedge2edge
{
	halfedge2edge(const Mesh& m, std::vector<edge_descriptor>& edges)
		: m_mesh(m), m_edges(edges)
	{}
	void operator()(const halfedge_descriptor& h) const
	{
		m_edges.push_back(edge(h, m_mesh));
	}
	const Mesh& m_mesh;
	std::vector<edge_descriptor>& m_edges;
};


/*

void poissonRec(std::string in, std::string outresult){
	// Poisson options
	FT sm_angle = 20.0; // Min triangle angle in degrees.
	FT sm_radius = 30; // Max triangle size w.r.t. point set average spacing.
	FT sm_distance = 0.25; // Surface Approximation error w.r.t. point set average spacing.
	// Reads the point set file in points[].
	// Note: read_xyz_points_and_normals() requires an iterator over points
	// + property maps to access each point's position and normal.
	// The position property map can be omitted here as we use iterators over Point_3 elements.
	PointList points;
	std::ifstream stream(in);
	if (!stream ||
		!CGAL::read_xyz_points_and_normals(
		stream,
		std::back_inserter(points),
		CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type())))
	{
		std::cerr << "Error: cannot read file data" << std::endl;
		return;
	}
	// Creates implicit function from the read points using the default solver.
	// Note: this method requires an iterator over points
	// + property maps to access each point's position and normal.
	// The position property map can be omitted here as we use iterators over Point_3 elements.
	Poisson_reconstruction_function function(points.begin(), points.end(),
		CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type()));
	// Computes the Poisson indicator function f()
	// at each vertex of the triangulation.
	if (!function.compute_implicit_function())
		return;
	// Computes average spacing
	FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(points.begin(), points.end(), 6);  //knn = 1 ring
	// Gets one point inside the implicit surface
	// and computes implicit function bounding sphere radius.
	PointEx inner_point = function.get_inner_point();
	Sphere bsphere = function.bounding_sphere();
	FT radius = std::sqrt(bsphere.squared_radius());
	// Defines the implicit surface: requires defining a
	// conservative bounding sphere centered at inner point.
	FT sm_sphere_radius = 5.0 * radius;
	FT sm_dichotomy_error = sm_distance*average_spacing / 1000.0; // Dichotomy error must be << sm_distance
	Surface_3 surface(function,
		Sphere(inner_point, sm_sphere_radius*sm_sphere_radius),
		sm_dichotomy_error / sm_sphere_radius);
	// Defines surface mesh generation criteria
	CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
		sm_radius*average_spacing,  // Max triangle size
		sm_distance*average_spacing); // Approximation error
	// Generates surface mesh with manifold option
	STr tr; // 3D Delaunay triangulation for surface mesh generation
	C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
	CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
		surface,                              // implicit surface
		criteria,                             // meshing criteria
		CGAL::Manifold_with_boundary_tag());  // require manifold mesh
	if (tr.number_of_vertices() == 0)
		return;

	// saves reconstructed surface mesh
	std::ofstream out(outresult);
	PolyhedronEx output_mesh;
	CGAL::output_surface_facets_to_polyhedron(c2t3, output_mesh);
	out << output_mesh;
	return;
}


void isotropicRemesh(std::string in, std::string outfilename){
	std::ifstream input(in);
	Mesh mesh;
	if (!input || !(input >> mesh)) {
		std::cerr << "Not a valid off file." << std::endl;
		return;
	}

	double target_edge_length = 0.04;
	unsigned int nb_iter = 3;
	std::cout << "Split border...";
	std::vector<edge_descriptor> border;
	PMP::border_halfedges(faces(mesh),
		mesh,
		boost::make_function_output_iterator(halfedge2edge(mesh, border)));
	PMP::split_long_edges(border, target_edge_length, mesh);
	std::cout << "done." << std::endl;
	PMP::isotropic_remeshing(
		faces(mesh),
		target_edge_length,
		mesh,
		PMP::parameters::number_of_iterations(nb_iter)
		.protect_constraints(true)//i.e. protect border, here
		);
	std::cout << "Remeshing done." << std::endl;
	std::ofstream cube_off(outfilename);
	cube_off << mesh;
}

*/