#include "lungModelling.h"

#ifndef _DOTOBJ_
#define _DOTOBJ_

//Graph

class mColor {
public:
	mColor() {
		R = 1.0;
		G = 1.0;
		B = 1.0;
		H = 0.0;
		S = 0.0;
		V = 1.0;
	}

	float R;
	float G;
	float B;
	float H;
	float S;
	float V;

	void initializeRBG(float R, float G, float B) {
		this->R = R;
		this->G = G;
		this->B = B;
		//float delta = max(max(R, G),B)- min(min(R, G), B);
		return;
	}
	void initializeHSV(float H, float S, float V) {
		this->H = H;
		this->S = S;
		this->V = V;
		return;
	}
};




class dotObj {
public:

	helper functions;
	cellArray cells;
	bbox boundingBox;
	analysis mAnalysis;

	facePropertyMap fMap;
	vertexPropertyMap vMap;

	//VARIABLE MEMBERS-------------------------------------------------------------------------------------------------
	std::string partName;
	std::vector<std::vector<float>> vertices;
	std::vector<std::vector<float>> normals;
	std::vector<std::vector<int>> faces;
	std::vector<std::vector<int>> edges;
	std::vector<std::vector<int>> lines;
	std::vector<std::vector<float>> vertexToNormalMatch;

	//Submodule analysis : Current selection of vertices
	std::vector<int> selectedVertices;
	std::vector<int> selectedNormals;
	std::vector<int> selectedFaces;
	std::vector<int> selectedSubset;
	Vector3f centerOfMass;

	std::vector<std::vector<float>> deformableMeshPressurePoints;
	std::vector<std::vector<int>> assemblyToChildPartVerticesMatch;

	Matrix<double, Dynamic, Dynamic> LaplacianMatrix;
	Matrix<double, Dynamic, Dynamic> WL;
	Matrix<double, Dynamic, Dynamic> WH;
	Matrix<double, Dynamic, Dynamic> WH0;
	Matrix<double, Dynamic, Dynamic> V;
	Matrix<double, Dynamic, Dynamic> initialPerVertexRingAreas;
	Matrix<double, Dynamic, Dynamic> adjacencyMatrix;

	//Sparce logic
	std::vector<std::vector<int>> adjacencyMatrixSimplified;
	SparseMatrix<double> adjacencyMatrixSparse;
	SparseMatrix<double> VSparse;
	SparseMatrix<double> initialPerVertexRingAreasSparse;
	SparseMatrix<double> LaplacianMatrixSparse;
	SparseMatrix<double> WLSparse;
	SparseMatrix<double> WHSparse;
	SparseMatrix<double> WH0Sparse;
	double initialAverageFaceArea, SL;

	//SURGERY CONNECTIVITY
	std::vector<float> QEMshapeCostMatrix;
	std::vector<float> QEMsamplingCostMatrix;
	std::vector<float> QEMtotalCostMatrix;
	std::vector<MatrixXf> Qs;

	//Geodesic distance
	std::vector<double> geoDistPerVertexFromFurthestVertex;
	std::vector< std::vector<double>> pointsOfInterest;

	//TIME
	clock_t * ticTocTimerBegin;
	clock_t * ticTocTimerCheckPoint;
	double ticTocTimerElapsedSecs;

	//PCA
	PCA attachedPCA;

	//LAPLACIAN MESH CONTRACTION
	double customFunctionFrequency = Pi / 2;
	double customFunctionAmplitude = 0.8;
	bool uniformContraction = false;
	int skeletonizationIterations = 8;
	int skeletonizationSelection = 1000;

	//SDF
	std::vector<float> sdf_property_map;
	std::vector<float> sdf_property_map_per_vertex;
	std::vector<std::string> segment_property_name;
	std::vector<int> segment_property_map_per_vertex;
	std::vector<int> segment_property_map_per_line;
	std::vector<int> segment_property_map; //Per Face
	
	

	//ONE-DIMENSIONAL
	std::vector<std::vector<int>> neighbouringVerticesPerVertexOnSkeleton;
	std::vector<std::vector<int>> neighbouringLinesPerVertexOnSkeleton;
	std::vector<int> skeletonInlets;
	std::vector<int> skeletonOutlets;

	std::vector<int> edgesToBeDeleted;
	std::vector<int> facesToBeDeleted;

	//Version II

	std::vector<mVertex> mVertices;
	std::vector<mFace> mFaces;



public:

	//SKELETONIZATION=========================================================================================
	void getSkel(void);
	Matrix<double, Dynamic, Dynamic> getLaplacian(std::string type);
	void setWL(void);
	void setWH(void);
	void setSL(void);
	void setSL(double value);
	void updateWL(void);
	void updateWH(void);
	void initializeSkeletonizationProcess(void);
	double computeOmega(int firstPoint, int secondPoint);
	double computeTriangeArea(int firstIndex, int secondIndex, int thirdIndex);
	double computeAngle(Vector3d a, Vector3d b);
	//Skeletonization-Narrowing Approach #1
	void skelStep();
	void skeletonize(int iterations);
	Matrix<double, Dynamic, Dynamic> solve(Matrix<double, Dynamic, Dynamic> & L, Matrix<double, Dynamic, Dynamic> & WL, Matrix<double, Dynamic, Dynamic> & WH, Matrix<double, Dynamic, Dynamic> & V);
	Matrix<double, Dynamic, Dynamic> getRingAreas(void);
	double getRingArea(int vertex);
	void convToMat(Matrix<double, Dynamic, Dynamic> &V, std::vector<std::vector<float>>&vert, int direction);
	dotObj mcfskel(void);
	dotObj mcfskel(bool contract, bool convertToSkeleton);
	//SPARSE SKELETONIZATION
	void skeletonizeSparse(std::string method, int iterations, double lamda);
	void skelStepSparse(int iteration, std::string method);
	SparseMatrix<double> getLaplacianSparse(std::string type);
	SparseMatrix<double> getRingAreasSparse(void);
	void initializeSparse(void);
	void setWHSparse(void);
	void updateWHSparse(void);
	void setWLSparse(void);
	void updateWLSparse(void);
	SparseMatrix<double> solveSparse(SparseMatrix<double> &LSparse, SparseMatrix<double>  &WLSparse, SparseMatrix<double>  &WHSparse, SparseMatrix<double>  &VSparse);
	void convToMatSparse(SparseMatrix<double> &V, std::vector<std::vector<float>>&vert, int direction);
	double computeOmegaSparse(int firstIndex, int secondIndex);
	std::vector<int> findCommonNeighboursSparse(int p1, int p2);
	//CONNECTIVITY SURGERY
	MatrixXf getK(int i, int j);
	MatrixXf getQ(int i);
	void getQs(void);
	float getShapeCost(int i, int j);
	float getSamplingCost(int i, int j);
	float getErrorMetricOfVertexToPoint(int i, Vector3f p);
	void halfEdgeCollapse(int selectedVertex, int withVertex);
	void computeShapeCosts(void);
	void computeSamplingCosts(void);
	float computeTotalCostForEdge(int i, int j);
	void computeTotalCosts(void);

	////UNIVERSAL======================================================================
	//TIME
	void tic();
	void toc();

	////3D MESH======================================================================
	//IO related

	void dotObj::exportToFileSegmentedPerLine(std::string fileType, std::string fileName);

	std::string  exportToFile(std::string name, bool overwrite = true);
	std::string  exportToXYZ(std::string name, bool overwrite = true);
	std::string  exportToFile(std::string name, std::string type, bool overwrite = true);


	void initializeFromFile(std::string thefile, bool verbose = true, bool readvertices = true, bool readnormals = true, bool readfaces = true, bool readedges = false);
	void initializeFromFile(std::string thefile, std::string fileType, bool verbose = true, bool readvertices = true, bool readnormals = true, bool readfaces = true, bool readedges = false);
	void initializeFromFile(std::stringstream & filestream, std::string fileType, bool verbose = true, bool readvertices = true, bool readnormals = true, bool readfaces = true, bool readedges = false);
	std::ostringstream initializeFromFileE(std::string thefile, bool verbose = true, bool readvertices = true, bool readnormals = true, bool readfaces = true, bool readedges = false);
	void initializeFromPLY(std::string fnamein, std::string fnameout = "", bool verbose = true, bool readvertices = true, bool readnormals = true, bool readfaces = true, bool readedges = false);
	void initializeFromVTK(std::string thefile, bool verbose = true, bool readvertices = true, bool readnormals = true, bool readfaces = true, bool readedges = false);


	void initializeDotObj(std::string filename, bool verbose=false);


	void exportSelectedIndicesToFile(std::string name);
	void toOFF(std::string name);
	std::stringstream toOFF();
	void exportToFileSegmented(std::string output);
	void exportToFileSegmentedPerFace(std::string output);
	void exportToFileSegmented(std::string output, std::string keyword);
	bool faceHasVertex(int faceIndex, int vertexIndex);
	void mapNormals2Faces(void);
	void getEdges(void);
	//ANALYSIS
	Vector3f findCenterOfMass(void);
	//DELTA COORDINATES
	Vector3f getDeltaCoordinate(int vertex);
	void getDeltaCoordinates(Matrix<double, Dynamic, Dynamic> &M);
	void getDeltaCoordinates(std::vector<Vector3f> &M);
	//SMOOTHING
	void smoothing(std::string type, std::string laplacian, int iterations, double lamda, double mu);
	void smoothingSparse(std::string type, std::string laplacian, int iterations, double lamda, double mu);
	void smoothingSparse(std::string type, std::string laplacian, int iterations, double lamda, double mu, std::vector<int> &selection);
	void taubin(int iterations, double lamda, double mu); //Smoothing a mesh
	//DIFF & INTERPOLATION
	void interp(dotObj initial, double percentage); //Interp in diff meshes
	//SELECTION
	void selectorBasedOnCenterline(int startingPointIndex, double distance);
	void selectorBasedOnGeodesic(int startingPointIndex, double distance);

	void extendSelection(void);
	void selector(int startingPointIndex, int numberOfPoints);
	void getNeighboursByFaceTri(int faceIndex, int numOfNeighbours);


	//TRANSLATION
	void rotate(Vector3f & target);
	void scale(double scaleX, double scaleY, double scaleZ);
	void translate(Vector3f & target);
	//NORMALS
	Vector3f calculateNormal(int faceIndexUnderInvestigation);
	void recalculateNormals(void);
	//FACES
	double getAverageFaceArea(void);
	void getFaceIndicesFromVertexList(std::vector<int> &v, std::vector<int> &f);
	//BOUNDING BOX
	void getBoundingBox(void);
	//CONNECTIVITY & ADJACENCY
	std::vector<int> findAdjacentFaces(int faceIndexUnderInvestigation);
	std::vector<int> findAdjacentFacesPerVertex(int vertexIndexUnderInvestigation);
	std::vector<int> searchImmediateNeighboors(int startingpoint);
	std::vector<int> searchImmediateLineNeighboors(int startingpoint);
	void searchpoints(int startingpoint);
	void searchpoints(int startingpoint, int numberOfPoints);
	std::vector<int> findCommonNeighbours(int p1, int p2);
	std::vector<int> getFaceIndexContainingVertexIndex(int vertexIndex);
	std::vector<int> getLineIndexContainingVertexIndex(int vertexIndex);

	template <typename  T>
	distance findMinimumDistance(std::vector<T> & r, std::vector<std::vector<T>> & F) {
		distance result;
		Vector3f x = Vector3f(r[0], r[1], r[2]);
		result.value = 1000;
		result.index = 0;

		for (int j = 0; j < F.size(); j++) {
			Vector3f p = Vector3f(F[j][0], F[j][1], F[j][2]);
			if ((p - x).norm() < result.value) { result.index = j; result.value = (p - x).norm(); }
		}

		return result;
	}

	template <typename  T>
	distance findMinimumDistance(std::vector<T> & r, std::vector<std::vector<T>> & F, std::vector<int> & indices) {
		distance result;
		Vector3f x = Vector3f(r[0], r[1], r[2]);
		result.value = 1000;
		result.index = 0;

		for (int m = 0; m < indices.size(); m++) {
			int j = indices[m];
			Vector3f p = Vector3f(F[j][0], F[j][1], F[j][2]);
			if ((p - x).norm() < result.value) { result.index = j; result.value = (p - x).norm(); result.nestedIndex = m; }
		}

		return result;
	}

	//SHAPE DIAMETER FUNCTION
	std::vector<float> getSDFPropertyMap(double cone, int number_of_rays, bool postprocess);
	void getSDF(double cone, size_t number_of_rays, bool postprocess, bool clustering, size_t number_of_clusters, double smoothing_lambda, bool exportClusters);
	void refineSDF(void);
	//INSIDE OR OUTSIDE
	bool checkIfPointisInside(std::vector<float> point);
	bool checkIfPointisInsidev2(Vector3f point);
	//SEGMENTATION
	void registerPart(dotObj &input, bool doNormals = false);
	void getPositionsFromChild(dotObj &input);
	//CELLS & OCTREE
	void initCells(float dx);

	//SHAPE ANALYSIS
	void segmentByGeneration(void);
	void segmentByGeneration(dotObj & skel);
	void graphBasedAnalysis(dotObj & skel);
	void graphBased1DModelAnalysis(bool LR=true , bool generations=true, int inlet=-1);
	void graphBasedPruning(void);
	void segmentBySkeleton(void);
	void segmentationApply(void);
	void updateSegmentProperties(void);

	//
	std::string segmentationMode = "SDF";
	int normal_est_neighb = 30;
	int smoothness_iter = 5;
	double smoothness_sharpness_angle = 70;
	int smoothness_neighbouthood_size = 50;

	//OTHER
	std::vector<float> locate_flat_tips_based_on_calculated_normal(bool doFaces);

	//CGAL
	void simplificationEdgeCollapse(int numberOfEdges = 1000000);
	void isotropicRemeshing(void);
	void refine_fair(void);

	void splitSkeletonEdges(float minimumDistance);

	//DECIMATION
	void searchAndRemoveMinimum(void);
	void removeVertex(int i);
	void removeNormal(int i);
	void removeFace(int i);
	void mergeCloseVertices(void);
	void mergeVertices(int i, int j);

	//INDEXING
	int getPointIndex(float x, float y, float z);
	//ADJACENCY MATRICES
	void getAdjacencyMatrices(bool dense, bool vector, bool sparse);
	MatrixXd getAdjacency(void);
	std::vector<std::vector<int>> getAdjacencySimplified(void);
	std::vector<std::vector<int>> getAdjacencyTripletList(void);
	SparseMatrix<double> getAdjacencySparse(void);
	SparseMatrix<double> getAdjacencySparseWeighted(void);

	SparseMatrix<double> getAdjacencySparse(std::vector<int> &selectionList);

	//GEODESIC
	std::vector<int>  geodesicSelectorBetweenPoints(int a, int b);
	std::vector<int>  geodesicSelectorBetweenPoints(std::vector<int> pointlist);
	std::vector<double> getMinDists(int point1);
	std::vector<double> getGeodesicDistanceFromSeedPoint(int index);
	std::vector<double> geodesicBasedWeighting(double Freq);
	std::vector<double> geodesicBasedWeighting(void);
	double getNorm(int p1, int p2);

	////ONE-DIMENTIONAL-GRAPH======================================================================
	void removeVertexFromSkeleton(int vertexIndex);
	void removeVertexFromSkeleton(std::vector<int> toDelete);
	void removeLineFromSkeleton(int vertexIndex);
	void getNeighbouringVerticesOnSkeleton(void); //Initial point not included
	std::vector<int> getOneRingNeighboursOnSkeleton(int skeletonVertexIndex);
	std::vector<int> getSecondRingNeighboursOnSkeleton(int skeletonVertexIndex);

	//PATHFINDER
	std::vector<std::vector<int>> getNeighboursPerVertex(void);
	std::vector<std::vector<int>> pathFinder(std::vector<std::vector<int>> &neighboursPerVertex);
	std::vector<std::vector<int>> pathFinderTwoStep(void);
	void pathFinderTwoStep(std::vector<std::vector<int>> &neighboursPerVertex, std::vector<std::vector<int>> &paths);
	std::vector<std::vector<int>> analyzePathsPerVertex();
	std::vector<int> dotObj::getBrachesToDeleteBasedOnPathSize(std::vector<std::vector<int>> &paths, std::vector<std::vector<int>> &neighboursPerVertex, int lim);
	void mcfskelRefine(void);
	void mcfskelRefineStepOne(void);
	void mcfskelRefineStepOne(int lim);
	void mcfskelRefineStepTwo(void);
	void buildInitialTrees(dotObj &out1, dotObj &out2, int gens = 3);

	//POINT CLOUD======================================================================
	//DENSITY
	double getDensity(int vertexIndex, float Radius);
	double getDensity(int vertexIndex, float Radius, std::vector<int> Indices);

	////NOT REQUIRED======================================================================

	int assign_triangle_to_group(void);
	int models_diff(void);
	int assign_close_triangle_to_group(void);

	//Other
	void graphBasedFlattening(ggraph &sk);
};

class simulation {
public:

	dotObj *tempModel;
	std::vector<dotObj> lungmodel;
	int	 seedPoint;
	double narrowingRatio;
	double localDIameterBeforeNarrowing;
	double localDIameterAfterNarrowing;
	status * stat;

public:
	helper functions;

	//int softBodyDeform(glottis glot, dotObj source, dotObj parent);

	void generateTree(
		int numOfGenerations,
		int density,
		float radiusMultiplier,
		volume &thevol,
		dotObj boundary,
		dotObj hostMesh,
		dotObj &oneDim,
		std::string id,
		bool build1D,
		bool build3D,
		bool useHostMeshForPointCloud,
		int mode);

	void generateTree(
		int numOfGenerations,
		int density,
		float radiusMultiplier,
		dotObj oneDRep,
		dotObj hostMesh,
		dotObj &treePCexp,
		std::string id,
		bool build3D,
		bool useHostMeshForPointCloud,
		int mode);

	void generate1DTree(int numOfGenerations, volume &thevol, dotObj &boundary, dotObj &hostMesh, dotObj &oneDim);

	void generate3DPointCloudFromTree(
		int numOfGenerations,
		int density,
		float radiusMultiplier,
		dotObj oneDRep,
		dotObj hostMesh,
		dotObj &treePCexp,
		std::string id,
		bool useHostMeshForPointCloud);

	void generate3DPointCloudFromTree(
		int density,
		float radiusMultiplier,
		dotObj oneDRep,
		dotObj &treePCexp,
		std::vector<float> diameters,
		std::string id);

	void generate3DPointCloudFromGraph(
		int numOfGenerations,
		int density,
		float radiusMultiplier,
		ggraph * theGraph,
		dotObj * treePCexp
	);

	void generateVolume(
		int numOfGenerations,
		int density,
		volume &thevol,
		dotObj boundary,
		std::string id);
	void segmentTree(dotObj &oneD, dotObj &recModel, std::string result, std::string keyword);
	void segmentTree2(dotObj &oneDL, dotObj &oneDR, dotObj &recModel, std::string result, std::string keyword);
	void flatten(dotObj &oneDL, dotObj &oneDR, dotObj &recModel, std::string result, std::string keyword);
	void narrow(
		double contractionPercentageBox,
		int contractionStrengthBox,
		bool useCustomFunction,
		double frequencyBox,
		bool useInterpolation,
		double interpolationPercentage,
		bool useExtendedSmoothing,
		int  numberOfRaysSDF,
		int numberOfClustersSDF,
		double lamdaSDF);
	void loadFile(std::string modelName);
	void writeFile(void);
	void exportSegmented(void);
	void exportFile(void);
	void exportVerticesToText(void);
	void sdf(
		std::string segmentOrCluster,
		double coneAngleSDF,
		int numberOfRaysSDF,
		bool postProcessingSDF,
		int numberOfClustersSDF,
		double lamdaSDF);
	void loadSelectionFromTXT(std::string filepath);
	void extendBronchialTree(void);
	void extendBronchialTree(
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
		bool verbose);

	void extendBronchialTreeV2(
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
		bool verbose);

	void extendBronchialTreeV2(
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
		volume * volLeft,
		volume * volRight
	);
};

class mIndex {
public:
	int x;
	int y;
	int z;

	mIndex(int x, int y, int z) {
		this->x = x;
		this->y = y;
		this->z = z;
	}

	mIndex(void) {
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}

	template<typename T>
	bool isInsideMatrix(std::vector<std::vector<std::vector<T>>> &A, int padding) {
		if ((this->x - padding >= 0) && (this->x + padding < A.size()))
		{
			if ((this->y - padding >= 0) && (this->y + padding < A[0].size()))
			{
				if ((this->z - padding >= 0) && (this->z + padding < A[0][0].size()))
				{
					return true;
				}
			}
		}

		return false;
	}
};

class mVertex {
public:
	Vector3f position;
	Vector3f normal;
	std::vector<mVertex*> neighbouringVertices;
};

class mEdge {
public:
	std::vector<mVertex*> vertices;
};

class mFace {
public:
	std::vector<mEdge*> edges;
	std::vector<mVertex*> vertices;
	std::vector<mFace*> neighbouringFaces;
	Vector3f centroid;
	Vector3f normal;
};

class mCellSimplified {
public:
	mIndex index;
	float value;

	mCellSimplified(int x, int y, int z, float value) {
		this->index = mIndex(x, y, z);
		this->value = value;
	}

	mCellSimplified(void) {
		this->index = mIndex(0, 0, 0);
		this->value = 0.0;
	}
};

template <class T>
class mProperty {
public:
	int id;
	T value;
	mProperty<T>(int id, T val) {
		this->id = id;
		this->value= val;
	}
	mProperty<T>(void) {
		this->id = 0;
		this->value = 0;
	}

	bool propertyHasID(int N) {
		if (this->id == N) {
			return true;
		}
		else {
			return false;
		}
	}
};

class pCell {
public:
	std::vector<mProperty<float>> properties;

	pCell(void) {
		this->properties.resize(14);
	}

	mProperty<float> maxValIndex(void) {
		int ind = 0;
		float mVal = this->properties[ind].value;
		for (int i = 0; i < this->properties.size(); i++) {
			if (this->properties[i].value > mVal) {
				mVal = this->properties[i].value;
				ind = i;
			}
		}
		return mProperty<float>(ind, mVal);
	}


	mProperty<float> minValIndex(void) {
		int ind = 0;
		//float mVal = this->properties[ind].value;
		float mVal = 20.0;
		for (int i = 0; i < this->properties.size(); i++) {
			if ((this->properties[i].value < mVal)&&(this->properties[i].value >0)) {
				mVal = this->properties[i].value;
				ind = i;
			}
		}
		return mProperty<float>(ind, mVal);
	}

};

class mCell {
public:
	mIndex index;
	Vector3f position;
	bbox boundingBox;
	std::vector<gnode*> nodes;
	int label;
	bool isLung;

	std::vector<mProperty<float>> properties;

	mCell(void) {
		this->properties.resize(10);
	}
};

class mDiscretization {
public:
	float dx;
	float dy;
	float dz;
	int nCX;
	int nCY;
	int nCZ;

	mDiscretization(int X, int Y, int Z) {
		nCX = X;
		nCY = Y;
		nCZ = Z;
	}

	mDiscretization(void) {
		nCX = 0;
		nCY = 0;
		nCZ = 0;
	}
};

class mGrid {
public:
	mDiscretization discretization;
	std::vector<std::vector<std::vector<mCell*>>> cells;

	bbox boundingBox;
};



class voxelSpace {
public:
	bool regionHasNotLabel(mIndex * ind, std::vector<std::vector<std::vector<int>>> &A, int padding, int label) {
		bool found = true;
		for (int i = -padding; i <= padding; i++) {
			for (int j = -padding; j <= padding; j++) {
				for (int k = -padding; k <= padding; k++) {
					if (A[ind->x + i][ind->y + j][ind->z + k] == 0) {
						found = false;
					}
				}
			}
		}
		return found;
	}


	bool regionHasOnlyLabel(mIndex * ind, std::vector<std::vector<std::vector<int>>> &A, int padding, int label) {
		bool found = false;
		for (int i = -padding; i <= padding; i++) {
			for (int j = -padding; j <= padding; j++) {
				for (int k = -padding; k <= padding; k++) {
					if (A[ind->x + i][ind->y + j][ind->z + k] == label) {
						found = true;
					}
				}
			}
		}

		return found;
	}


	void Matrix2SparseRepresentationByValue(std::vector<std::vector<std::vector<int>>> &A, std::vector<mCellSimplified> &Asparse, float value) {
		Asparse.clear();
		for (int i = 0; i < A.size(); i++) {
			for (int j = 0; j < A[i].size(); j++) {
				for (int k = 0; k < A[i][j].size(); k++) {
					if (A[i][j][k] == value) {
						Asparse.push_back(mCellSimplified(i, j, k, value));
					}
				}
			}
		}

		return;
	}
	

	template<typename T>
	void resize3DMat(std::vector<std::vector<std::vector<T>>> &A, mDiscretization d) {
		A.resize(d.nCX);
		for (int i = 0; i < A.size(); i++) {
			A[i].resize(d.nCY);
		}

		for (int i = 0; i < A.size(); i++) {
			for (int j = 0; j < A[i].size(); j++) {
				A[i][j].resize(d.nCZ);
			}
		}

		return;
	}





	int findSeed(std::vector<std::vector<std::vector<int>>> &A, std::vector<mCellSimplified>Asparse, int label) {

		int seed;
		bool found = false;
		int a, x, y, z;
		a = -2;
		int b;
		b = a + 5;

		for (seed = 0; ((seed < Asparse.size()) && !found); seed++) {
			if (!found) {
				x = Asparse[seed].index.x;
				y = Asparse[seed].index.y;
				z = Asparse[seed].index.z;

				if ((x + a >= 0) && (x + b < A.size()))
				{
					if ((y + a >= 0) && (y + b < A[0].size()))
					{
						if ((z + a >= 0) && (z + b < A[0][0].size()))
						{
							found = true;

							for (int i = a; i < b; i++) {
								for (int j = a; j < b; j++) {
									for (int k = a; k < b; k++) {
										if (A[x + i][y + j][z + k] != label) {
											found = false;
										}
									}
								}
							}
						}
					}
				}
			}
		}

		return seed;
	}



};







#endif _DOTOBJ_