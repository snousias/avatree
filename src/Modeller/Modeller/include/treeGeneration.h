#include "lungModelling.h"

#ifndef _TREEGENERATIONCLASS_
#define _TREEGENERATIONCLASS_

class treeGeneration {
public:
	helper functions;
	multiBodyObj  tree;
	dotObj   oneDimension;
	dotObj   treePC;
	dotObj   reconstructedModel;
	std::vector<std::vector<std::vector<float>>> pointCloudPerLine;
	std::vector<float> diameterPerLineSegment;
	float limCutOffFactor;
	float overshootPercentageFactor;
	int outletCounter;
	std::vector<double> anatomyBasedRadius;

public:
	treeGeneration();
	int getCurrentGeneration(int input);
	int isLeftOrRight(int input);
	void volumeFillingSupplyDemand(volume v, dotObj & p, int numOfBranches);
	void assignDiameters(void);
	void assignDiameters(std::vector<float> diameters);
	void generateCloudForLineSegment(
		std::vector<std::vector<float>> &geom,
		std::vector<std::vector<std::vector<float>>> &geomPerLine,
		int & i,
		Vector3f * Points,
		int pointA,
		int pointB,
		int * v,
		int * vindex,
		int density,
		float radiusMultiplier,
		float dispersion,
		float overshootPercentage,
		bool useUniformDistance
		);
		

	void generateCloudForLineSegment(
		dotObj * geom,
		Vector3f start,
		Vector3f end,
		int density,
		float radiusMultiplier,
		float dispersion,
		float overshootPercentage,
		bool useUniformDistance,
		float R1,
		float R2
	);

	void generateCloudForLineSegment(
		std::vector<std::vector<float>> &geom,
		std::vector<std::vector<std::vector<float>>> &geomPerLine,
		Vector3f start,
		Vector3f end,
		//int pointIndex,
		int lineIndex,
		int density,
		float radiusMultiplier,
		float dispersion,
		float overshootPercentage,
		bool useUniformDistance
	);

	void generateCloudForLineSegment(
		std::vector<std::vector<float>> &geom,
		std::vector<std::vector<std::vector<float>>> &geomPerLine,
		Vector3f start,
		Vector3f end,
		int * v,
		int density,
		float radiusMultiplier,
		float dispersion,
		float overshootPercentage,
		bool useUniformDistance
	);







	void generateCloudForLineTip(
		std::vector<std::vector<float>> &geom,
		std::vector<std::vector<std::vector<float>>> &geomPerLine,
		int & i,
		Vector3f * Points,
		int pointA,
		int pointB,
		int * v,
		int * vindex,
		int density,
		float radiusMultiplier,
		float dispersion,
		float overshootPercentage,
		bool useUniformDistance
	);



	void generateCloudForLineTip(
		std::vector<std::vector<float>> &geom,
		std::vector<std::vector<std::vector<float>>> &geomPerLine,
		Vector3f start,
		Vector3f end,
		//int pointIndex,
		int lineIndex,
		int density,
		float radiusMultiplier,
		float dispersion,
		float overshootPercentage,
		bool useUniformDistance
	);

	void treeGeneration::generateCloudForLineTip(
		dotObj * geom,
		Vector3f start,
		Vector3f end,
		int density,
		float radiusMultiplier,
		float dispersion,
		float overshootPercentage,
		bool useUniformDistance,
		float R2
	);

	void get3DPointCloudv3(int density, float radius);
	void get3DPointCloudv3(bool useHost, dotObj & h, int numgen, int density, float radius);
	void cleanUpv3();
	void cleanUpv3(float radiusMultiplier);
	void naming(std::string result, std::string keyword);
	void naming2(std::string result, std::string keyword);
	void flatten(std::string result, std::string keyword);
};

#endif _TREEGENERATIONCLASS_