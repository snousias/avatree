#include "lungModelling.h"

#ifndef _OBJ_
#define _OBJ_

class status{
public:
	status();
	std::string action;
	float progress;
	bool isComplete;
	bool error;
};

class bbox{
public:
	helper functions;
	float maxY;
	float maxX;
	float maxZ;
	float minY;
	float minX;
	float minZ;
	bbox();
};

class plane{
public:
	helper functions;
	Vector3f pointOnPlane;
	Vector3f normal;
	bool checkWhichSideOfPlaneIsPoint(Vector3f pointOffPlane);
};

class PCA{
public:
	helper functions;
	PCA();  // This is the constructor
	std::vector<std::vector<float>> setOfPoints;
	std::vector<Vector3f> eigenVectors;
	std::vector<float> eigenValues;
	std::vector<Vector3f> eigenVectorsSorted;
	std::vector<float> eigenValuesSorted;
	Vector3f centerOfMass;

	Vector3f getEigenVector(int order);
	void getPlane(int planeOne, int planeTwo);
	Vector3f getAverage();
	Vector3f getAverage(std::vector<std::vector<float>> &points);
	void getPCA(void);
	void getPCA(std::vector<std::vector<float>> &points);
};

class multiBodyObj{
public:
	helper functions;
	std::vector<dotObj> scene;
	void addObject(std::string partName, dotObj &input);
	void removeObject(std::string partName);
	void removeObject(int id);
	std::string exportToFile(std::string name);
};

class normal{
public:
	helper functions;
	Vector3d vect;
};

class ray{
public:
	helper functions;
	Vector3d pointOfOrigin;
	Vector3d direction;
};

class volume{
public:
	helper functions;
	Vector3f startingPoint;
	int startingPointIndex; //Just in case
	Vector3f bifurcationPoint;
	int bifurcationPointIndex;
	Vector3f centerOfMass;
	int parentIndex;
	std::vector<int> childIndex;

	std::vector<std::vector<float>> vertices;
	std::vector<volume> subsets;
	plane splittingPlane;
	//Methods
	void defineSplittingPlane(std::string method="PCA");
	void defineSplittingPlane(Vector3f & p1, Vector3f & p2, Vector3f & p3);
	Vector3f setCenterOfMass(void);
	void splitVolume(plane splittingPlane, volume &subVolumeA, volume &subVolumeB);
	void createVolumeFromObj(dotObj boundaries, int perEdge,int maxPoints=0);
	void setStartingPoint(Vector3f point);
	void splitVolume(Vector3f startingPoint);
	void exportToFile(std::string k);
	void initializeFromFile(std::string thefile);
};

class world{
	helper functions;
};

class sorter{
public:
	helper functions;
	std::vector<float> toBeSorted;
	std::vector<float> Sorted;
	void sort(void);
};




class distance{
public:
	float value;
	int index;
	int nestedIndex;
	int index3D[3];
};

class gnode{
public:
	gnode();
	bool isBifurcation;
	bool isTerminal;
	bool isInlet;
	bool isMedian;
	bool isLeft;
	bool isRight;
	bool doStop;
	float diameter;
	float diameterBranch;
	int generation;
	int segmentIndex;
	int index;
	int line;
	int HorsfieldOrder;
	int StrahlerOrderStage_1;
	int StrahlerOrderStage_2;
	float length;
	float lengthOfBranch;
	float thetaToParent;
	float thetaToSimbling;
	float RB;
	float RD;
	float RL;

	float RBS;
	float RDS;
	float RLS;

	float RBH;
	float RDH;
	float RLH;

	std::vector<float> diametersToProcess;

	

	std::vector<int> meshVerticesIndices;
	std::vector<Vector3f> meshVerticesPositions;
	std::vector<Vector3f> Vertices;


	Vector3f position;
	Vector3f positionOfPreviousNode;


	int previous;
	int previousBifurcation;
	std::vector<int> next;
	std::vector<int> nextBifurcation;
	

	gnode * previousNodePtr;
	std::vector<gnode *> nextNodesPtr;
	gnode * nextMedianPtr;
	gnode * previousBifurcationPtr;
	std::vector<gnode *> nextBifurcationPtr;

	gnode * correspondingNodePtr;

	void copyNodeProperties(gnode * source, gnode * target);

};

class vertexProperty{
public:
	int generation;
	int segmentIndex;
	bool isOutlet = false;
	bool isInlet = false;
	bool isLeft = false;
	bool isRight = false;
};

class vertexPropertyMap{
public:

	std::vector<vertexProperty> map;
};

class faceProperty{
public:

	int generation;
	int segmentIndex;
	bool isOutlet = false;
	bool isInlet  = false;
	bool isLeft = false;
	bool isRight = false;



};

class facePropertyMap{
public:
	std::vector<faceProperty> map;

};



class gline{
public:
	int nodeA;
	int nodeB;
};

class gpath{
public:
	int startingNode;
	int endingNode;
	std::vector<int> nodes;
};

class ggraph{
public:
	std::vector<gnode> nodes;
	gnode * init;

	gnode * RU;
	gnode * RM;
	gnode * RL;
	gnode * LU;
	gnode * LL;

	void ggraph2Model(dotObj * O, bool jump2Bifurcations = false,std::string colorize="None");
	void ggraph2ModelColorize(dotObj * O, bool jump2Bifurcations = false, std::string colorize = "None");
	void ggraph2ModelV(dotObj * result);
	void enrichGraph(float minimumDistance);
	void exportGraphFeatures(std::string outfile);
	void exportGraphFeatures(std::string outfile, gnode* inode);
	ggraph generateGraphFromGraph(void);
	void getLobes(void);

	void initializeDiameters(void);
	void generateGraphFeaturesGeneration(void);
	void generateGraphFeaturesDiameter(void);
	void generateGraphFeaturesRatios(void);
	void generateGraphFeaturesLRDiscrimination(void);
	void generateGraphFeaturesLength(std::string type="branch");
	void propagateGraphFeatures(void);
};

class analysis{
public:
	bool parseTerminals = true;

	ggraph graph;
	gnode * inletPtr;

	int inlet;

	dotObj * model;
	dotObj * skel;

	std::vector<std::vector<int>> neighboursPerVertex;
	std::vector<gpath> paths;
	std::vector<std::vector<int>> branches;
	std::vector<std::vector<int>> fullpaths;
	std::vector<std::vector<int>> fullpathsSFInlet;

	std::vector<int> MeshVertex2SkeletonVertex;
	std::vector<int> MeshVertex2SkeletonEdge;
	std::vector<int> MeshFace2SkeletonVertex;
	std::vector<float> LocalDiameterPerMeshVertex;
	std::vector<std::vector<int>> SkeletonEdge2MeshVertices;
	std::vector<float> LocalDiameterPerSkeletonEdge;
	std::vector<int> SkeletonVertex2Generation;

	analysis(void);
	distance getMaxPath(void);


	void getNeighboursPerVertex(dotObj &skel);
	void pathFinder(std::vector<std::vector<int>> &neighboursPerVertex);
	void pathAnalyzer(std::vector<std::vector<int>> &neighboursPerVertex);



	void pathAnalyzer(dotObj &skel);
	void generateGraph(dotObj &skel, int inlet);
	void analyzeGraph(void);
	//void analyzeGraph2(void);
	void generateGraphFeatures(void);
	
	void generateGraphFeaturesLRDiscrimination(void);
	void generateGraphFeaturesGeneration(void);
	
	void initialize(dotObj * geom, dotObj * skel);
	void getSkeletonVertex2Generation(dotObj &skel);

	void exportGraphFeatures(std::string outfile);
	
};







#endif
