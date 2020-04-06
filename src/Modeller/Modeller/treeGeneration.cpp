#include "lungModelling.h"

treeGeneration::treeGeneration() {
	//this->anatomyBasedRadius = { 14.72, 11.68, 9.27, 8.16, 5.84, 2.64, 2.18, 1.82, 1.62, 1.34, 1.16, 0.8, 0.6, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2 };

	this->anatomyBasedRadius = { 12.0, 12.0, 8.63, 5.6, 4.5, 3.5, 2.8, 2.3, 1.86, 1.54, 1.3, 1.09, 0.95, 0.83, 0.74, 0.5, 0.49, 0.4, 0.38, 0.36, 0.34, 0.31, 0.29, 0.25 };

	//this->anatomyBasedRadius = { 12.0, 11.0, 8.0, 7.3, 5.9, 5.9, 2.6, 3.0, 2.2, 1.9, 1.6, 1.25, 1.1, 0.9, 0.8, 0.6, 0.3, 0.2, 0.1 };
}

int treeGeneration::getCurrentGeneration(int input) {
	int k = 0;
	double a = input;
	while (a >= 2) {
		a = a / (2.0);
		k++;
	}

	return k;
}

int treeGeneration::isLeftOrRight(int input) {
	int side = 0;

	int gen = getCurrentGeneration(input);
	int numberOfGenerationAirways = (int)(powf(2.0, (float)gen));
	int residue = input - numberOfGenerationAirways + 1;
	double sideMeasurement = (double)residue / (double)numberOfGenerationAirways;

	if (sideMeasurement > 0.5) { side = 1; }
	else { side = -1; }

	//std::cout << "  airway:" << input << "  side: " << side << "  gen:" << gen << "  numberOfGenerationAirways:" << numberOfGenerationAirways << std::endl;

	return side;
}

void treeGeneration::volumeFillingSupplyDemand(volume _v, dotObj & startingMesh, int numOfBranches = 100) {
	float branchingFraction = 0.40;
	float maxAngleDegrees = 60;
	float minLen = 1.8;
	float maxAngleRad = maxAngleDegrees * 2 * Pi / 360;
	float sinmaxangle = sin(maxAngleRad);

	std::vector <volume>vols;
	int lastPtr = 0;

	//Reset
	dotObj nullObj;
	this->oneDimension = nullObj;
	oneDimension.lines.clear();
	//Initialize starting point
	Vector3f _startingPoint = { startingMesh.vertices.at(0).at(0), startingMesh.vertices.at(0).at(1), startingMesh.vertices.at(0).at(2) };

	//Initialize volume
	_v.startingPoint = { startingMesh.vertices.at(0).at(0), startingMesh.vertices.at(0).at(1), startingMesh.vertices.at(0).at(2) };
	_v.startingPointIndex = 0;
	vols.push_back(_v);

	this->oneDimension.vertices.push_back({ _startingPoint[0], _startingPoint[1], _startingPoint[2] });

	//Init diameters=========

	std::vector<std::vector<int>> examinationTicket;
	examinationTicket.push_back({ startingMesh.lines.at(0).at(0), startingMesh.lines.at(0).at(1) });

	for (int i = 0; i < startingMesh.lines.size(); i++) {
		int lp0 = examinationTicket.at(i).at(0);
		int lp1 = examinationTicket.at(i).at(1);
		int p0 = examinationTicket.at(i).at(0) - 1;
		int p1 = examinationTicket.at(i).at(1) - 1;

		std::vector<int> neighboors;
		neighboors = startingMesh.searchImmediateLineNeighboors(lp1);
		std::sort(neighboors.begin(), neighboors.end());
		Vector3f p0vec = { startingMesh.vertices.at(p0).at(0), startingMesh.vertices.at(p0).at(1), startingMesh.vertices.at(p0).at(2) };
		Vector3f p1vec = { startingMesh.vertices.at(p1).at(0), startingMesh.vertices.at(p1).at(1), startingMesh.vertices.at(p1).at(2) };

		if (neighboors.size() == 3) {
			int p2 = neighboors.at(neighboors.size() - 2) - 1;
			int p3 = neighboors.at(neighboors.size() - 1) - 1;

			examinationTicket.push_back({ p1 + 1, p2 + 1 });
			examinationTicket.push_back({ p1 + 1, p3 + 1 });

			Vector3f p2vec = { startingMesh.vertices.at(p2).at(0), startingMesh.vertices.at(p2).at(1), startingMesh.vertices.at(p2).at(2) };
			Vector3f p3vec = { startingMesh.vertices.at(p3).at(0), startingMesh.vertices.at(p3).at(1), startingMesh.vertices.at(p3).at(2) };

			//Bifurcation Point--------
			Vector3f bifucrationPoint = p1vec;

			//Append vertex
			this->oneDimension.vertices.push_back({ bifucrationPoint[0], bifucrationPoint[1], bifucrationPoint[2] });

			//Append Line------
			oneDimension.lines.push_back({ vols.at(i).startingPointIndex + 1, (int)this->oneDimension.vertices.size() - 1 + 1 });

			//Center of Mass---------
			vols.at(i).setCenterOfMass();

			//Calculate bifurcation
			vols.at(i).bifurcationPointIndex = p1;
			vols.at(i).bifurcationPointIndex = (int)this->oneDimension.vertices.size() - 1;
			Vector3f an = vols.at(i).centerOfMass - p0vec;
			Vector3f bn = p1vec - p0vec;
			float projab = an.dot(bn) / (bn.norm());
			Vector3f newCenterOfMass0 = p0vec + projab * (p1vec - p0vec) / ((p1vec - p0vec).norm());
			//Vector3f newCenterOfMass = p0vec + 2.5*(p1vec - p0vec);
			vols.at(i).centerOfMass = newCenterOfMass0;

			//Splitting plane
			vols.at(i).defineSplittingPlane(p1vec, p2vec, p3vec);

			Vector3f a = vols.at(i).splittingPlane.normal;
			Vector3f b = vols.at(i).centerOfMass - vols.at(i).startingPoint;
			Vector3f a1 = (a.dot((b / b.norm())))*b / b.norm();
			Vector3f a2 = a - a1;
			a2.normalize();
			vols.at(i).splittingPlane.normal = a2;

			//Set subvols
			volume subVolA, subVolB;
			vols.at(i).splitVolume(vols.at(i).splittingPlane, subVolA, subVolB);

			subVolA.startingPoint = bifucrationPoint;
			subVolA.startingPointIndex = vols.at(i).bifurcationPointIndex;
			subVolA.parentIndex = i;

			subVolB.startingPoint = bifucrationPoint;
			subVolB.startingPointIndex = vols.at(i).bifurcationPointIndex;
			subVolB.parentIndex = i;

			vols.push_back(subVolB);
			vols.at(vols.size() - 1).childIndex.push_back(vols.size() - 1);

			vols.push_back(subVolA);
			vols.at(vols.size() - 1).childIndex.push_back(vols.size() - 1);
		}

		if (neighboors.size() == 1) {
			//Bifurcation Point
			Vector3f bifucrationPoint = p1vec;

			this->oneDimension.vertices.push_back({ bifucrationPoint[0], bifucrationPoint[1], bifucrationPoint[2] });
			//oneDimension.lines.push_back({ p0 + 1, (int)this->oneDimension.vertices.size() - 1 + 1 });
			oneDimension.lines.push_back({ vols.at(i).startingPointIndex + 1, (int)this->oneDimension.vertices.size() - 1 + 1 });

			//Center of Mass
			vols.at(i).setCenterOfMass();
			vols.at(i).bifurcationPointIndex = p1;
			vols.at(i).bifurcationPointIndex = (int)this->oneDimension.vertices.size() - 1;

			Vector3f an = vols.at(i).centerOfMass - p0vec;
			Vector3f bn = p1vec - p0vec;
			float projab = an.dot(bn) / (bn.norm());
			Vector3f newCenterOfMass0 = p0vec + projab * (p1vec - p0vec) / ((p1vec - p0vec).norm());

			vols.at(i).centerOfMass = newCenterOfMass0;

			//Splitting plane
			vols.at(i).defineSplittingPlane();
			Vector3f a = vols.at(i).splittingPlane.normal;
			Vector3f b = vols.at(i).centerOfMass - vols.at(i).startingPoint;
			Vector3f a1 = (a.dot((b / b.norm())))*b / b.norm();
			Vector3f a2 = a - a1;
			a2.normalize();
			vols.at(i).splittingPlane.normal = a2;

			//Set subvols
			volume subVolA, subVolB;
			vols.at(i).splitVolume(vols.at(i).splittingPlane, subVolA, subVolB);

			subVolA.startingPoint = bifucrationPoint;
			subVolA.startingPointIndex = vols.at(i).bifurcationPointIndex;
			subVolA.parentIndex = i;

			subVolB.startingPoint = bifucrationPoint;
			subVolB.startingPointIndex = vols.at(i).bifurcationPointIndex;
			subVolB.parentIndex = i;

			vols.push_back(subVolB);
			vols.at(vols.size() - 1).childIndex.push_back(vols.size() - 1);

			vols.push_back(subVolA);
			vols.at(vols.size() - 1).childIndex.push_back(vols.size() - 1);
		}

		lastPtr = i;
	}

	for (int i = lastPtr + 1; i < vols.size(); i++) {
		if ((this->oneDimension.vertices.size() < numOfBranches) && (vols.at(i).vertices.size() > 0)) {
			vols.at(i).setCenterOfMass();

			Vector3f uk = vols.at(i).centerOfMass - vols.at(i).startingPoint;

			Vector3f uknorm = uk;
			uknorm.normalize();
			Vector3f vk = vols.at(i).startingPoint - vols.at(vols.at(i).parentIndex).startingPoint;
			Vector3f vknorm = vk;
			vknorm.normalize();
			Vector3f uknew;
			Vector3f projj;
			float theta = acos(((uknorm).dot(vknorm)) / ((uknorm).norm()*(vknorm).norm()));
			if (theta > maxAngleRad) {
				
				//Vector3f projection = this->functions.projPointOnVector(vols.at(i).centerOfMass, vols.at(vols.at(i).parentIndex).startingPoint, vols.at(i).startingPoint);
				
				
				//uknew = vols.at(vols.at(i).parentIndex).startingPoint +
				//	projection +
				//	tan(maxAngleRad)*(vols.at(i).startingPoint - projection).norm()*(vols.at(i).centerOfMass - projection).normalized();
				projj = (uk.dot(vk) / vk.squaredNorm())*vk;
				uknew = vols.at(i).startingPoint + projj + sin(maxAngleRad) * (uk - (projj));
				vols.at(i).centerOfMass = uknew;
			}


			Vector3f fulldist = (vols.at(i).centerOfMass - vols.at(i).startingPoint);

			Vector3f bifucrationPoint = vols.at(i).startingPoint + branchingFraction * fulldist;



			

			if (((bifucrationPoint - vols.at(i).centerOfMass).norm() > minLen)) {   //<=Caution length
				this->oneDimension.vertices.push_back({ bifucrationPoint[0], bifucrationPoint[1], bifucrationPoint[2] });
				oneDimension.lines.push_back({ vols.at(i).startingPointIndex + 1, (int)this->oneDimension.vertices.size() - 1 + 1 });

				vols.at(i).bifurcationPointIndex = (int)this->oneDimension.vertices.size() - 1;

				//Splitting plane and splitting plane correction
				vols.at(i).defineSplittingPlane();
				Vector3f a = vols.at(i).splittingPlane.normal;
				Vector3f b = vols.at(i).centerOfMass - vols.at(i).startingPoint;
				Vector3f a1 = (a.dot((b / b.norm())))*b / b.norm();
				Vector3f a2 = a - a1;
				a2.normalize();
				vols.at(i).splittingPlane.normal = a2;

				if (vols.at(i).vertices.size() >1) {
					volume subVolA, subVolB;
					vols.at(i).splitVolume(vols.at(i).splittingPlane, subVolA, subVolB);

					subVolA.startingPoint = bifucrationPoint;
					subVolA.startingPointIndex = vols.at(i).bifurcationPointIndex;
					subVolA.parentIndex = i;

					subVolB.startingPoint = bifucrationPoint;
					subVolB.startingPointIndex = vols.at(i).bifurcationPointIndex;
					subVolB.parentIndex = i;

					vols.push_back(subVolB);
					vols.at(vols.size() - 1).childIndex.push_back(vols.size() - 1);

					vols.push_back(subVolA);
					vols.at(vols.size() - 1).childIndex.push_back(vols.size() - 1);
				}
			}
		}
	}

	//Corrections
	for (int i = 0; i < this->oneDimension.lines.size(); i++) {
		int v0 = this->oneDimension.lines.at(i).at(0) - 1;
		int v1 = this->oneDimension.lines.at(i).at(1) - 1;

		//int gen = this->getCurrentGeneration(v1);

		bool isTerminal = true;

		for (int j = 0; j < this->oneDimension.lines.size(); j++) {
			int v2 = this->oneDimension.lines.at(j).at(0) - 1;

			if (v2 == v1) { isTerminal = false; }
		}

		if (isTerminal) {
			Vector3f p0, p1;
			p0[0] = this->oneDimension.vertices.at(v0).at(0);
			p0[1] = this->oneDimension.vertices.at(v0).at(1);
			p0[2] = this->oneDimension.vertices.at(v0).at(2);

			p1[0] = this->oneDimension.vertices.at(v1).at(0);
			p1[1] = this->oneDimension.vertices.at(v1).at(1);
			p1[2] = this->oneDimension.vertices.at(v1).at(2);
			int length = (p1 - p0).norm();
			if (length < 20) {		//<--Caution
				//p1 = p0 + 2 * (p1 - p0);
			}
			else {
				//p1 = p0 + 2 * (p1 - p0);
			}

			this->oneDimension.vertices.at(v1).at(0) = p1[0];
			this->oneDimension.vertices.at(v1).at(1) = p1[1];
			this->oneDimension.vertices.at(v1).at(2) = p1[2];
		}
	}
	return;
}

void treeGeneration::assignDiameters(void) {
	diameterPerLineSegment.clear();
	diameterPerLineSegment.resize(this->oneDimension.lines.size());
	for (int i = 0; i < this->oneDimension.lines.size(); i++) {
		int v[4];
		int vindex[3];
		int vindexptr = 0;
		int vptr = 0;
		vindex[vindexptr++] = i;
		v[vptr++] = oneDimension.lines.at(i).at(0) - 1;
		v[vptr++] = oneDimension.lines.at(i).at(1) - 1;
		for (int g = 0; g < oneDimension.lines.size(); g++) {
			if (oneDimension.lines.at(g).at(0) - 1 == v[1]) {
				v[vptr++] = oneDimension.lines.at(g).at(1) - 1;
				vindex[vindexptr++] = g;
			}
		}
		int genParent = getCurrentGeneration(v[1]);

		this->diameterPerLineSegment[i] = this->anatomyBasedRadius.at(genParent);
	}

	return;
}

void treeGeneration::assignDiameters(std::vector<float> diameters) {
	diameterPerLineSegment.clear();
	diameterPerLineSegment.resize(this->oneDimension.lines.size());
	diameterPerLineSegment = diameters;

	return;
}




void treeGeneration::generateCloudForLineSegment(
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
) {
	Vector3f vecFinal = end - start;
	int numOfPOints = (int)(density * vecFinal.norm());
	float angle;
	float overshoot;
	float dist;
	Vector3f zpositive = { 0.0, 0.0, 1.0 };
	Vector3f ypositive = { 0.0, 1.0, 1.0 };
	Vector3f xpositive = { 1.0, 0.0, 0.0 };
	Vector3f vecOfRef;
	Vector3f axis;
	Quaternion<float> q;
	Vector3f point;
	Vector3f result;

	geom->vertices.resize(numOfPOints);
	for (int j = 0; j < numOfPOints; j++) {
		angle = (2 * Pi)*functions.uniformRandom();  //Angle in rad
		overshoot = (overshootPercentage*R1 + vecFinal.norm()) / vecFinal.norm();
		if (useUniformDistance) {
				dist = overshoot * functions.uniformRandom() - ((overshoot - 1.0) / 2.0);
		}
		else {
			dist = 0.5 + dispersion * functions.normalRandom();
		}
		//SETAXIS
		if (abs(vecFinal.dot(zpositive)) != vecFinal.norm()*zpositive.norm()) {
			vecOfRef = vecFinal.cross(zpositive);
		}
		else if (abs(vecFinal.dot(ypositive)) != vecFinal.norm()*ypositive.norm()) {
			vecOfRef = vecFinal.cross(ypositive);
		}
		else if (abs(vecFinal.dot(xpositive)) != vecFinal.norm()*xpositive.norm()) {
			vecOfRef = vecFinal.cross(xpositive);
		}
		//DEFINEPOINT
		vecOfRef.normalize();
		axis = vecFinal;
		axis.normalize();
		q = AngleAxis<float>(angle, axis);
		point = q * vecOfRef;
		point.normalize();
		result = start + dist * vecFinal + ((R1 + dist * (R2 - R1))*point);
		geom->vertices[j] = { result[0], result[1], result[2] };
	}
	return;
}



void treeGeneration::generateCloudForLineSegment(
	std::vector<std::vector<float>> &geom,
	std::vector<std::vector<std::vector<float>>> &geomPerLine,
	Vector3f start,
	Vector3f end,
	int *v,
	int density,
	float radiusMultiplier,
	float dispersion,
	float overshootPercentage,
	bool useUniformDistance
)
{
	int startingPOint;
	int numOfPOints;
	int startingPOintPerLine;
	float radiusFinal;
	float radiusFinal2;
	float dist;
	float angle;
	float overshoot;
	Vector3f zpositive = { 0.0, 0.0, 1.0 };
	Vector3f ypositive = { 0.0, 1.0, 1.0 };
	Vector3f xpositive = { 1.0, 0.0, 0.0 };
	Vector3f vecOfRef;
	Vector3f axis;
	Quaternion<float> q;
	Vector3f point;
	Vector3f result;
	Vector3f vecFinal = end - start;
	int lineIndex = v[0];

	//std::cout << "Point Index " << pointIndex << "Line Index " << lineIndex << std::endl;

	startingPOint = geom.size();
	numOfPOints = (int)(density * vecFinal.norm());
	geom.resize(geom.size() + numOfPOints);
	startingPOintPerLine = geomPerLine[lineIndex].size(); //pointIndex
	geomPerLine[lineIndex].resize(startingPOintPerLine + numOfPOints);
#pragma omp parallel shared(vecFinal) private(q, point,result,vecOfRef,angle,overshoot,dist,axis,radiusFinal,radiusFinal2)
#pragma omp for
	for (int j = startingPOint; j < startingPOint + numOfPOints; j++) {
		//RADIUS========
		radiusFinal = radiusMultiplier * this->diameterPerLineSegment[v[0]];
		radiusFinal2 = radiusMultiplier * this->diameterPerLineSegment[v[1]];
		//ANGLE=========
		angle = (2 * Pi)*functions.uniformRandom();  //Angle in rad
		//DISTANCE=======
		overshoot = (overshootPercentage*radiusFinal + vecFinal.norm()) / vecFinal.norm();
		if (useUniformDistance) {
			if (lineIndex == 0) {
				//dist = 0.55*overshoot*functions.uniformRandom() - ((overshoot - 1.0) / 2.0) + 0.45;  //0.5 is important
				dist = overshoot * functions.uniformRandom() - ((overshoot - 1.0) / 2.0);
			}
			if (lineIndex > 0) {
				dist = overshoot * functions.uniformRandom() - ((overshoot - 1.0) / 2.0);
			}
		}
		else {
			dist = 0.5 + dispersion * functions.normalRandom();
		}
		//SETAXIS
		if (abs(vecFinal.dot(zpositive)) != vecFinal.norm()*zpositive.norm()) {
			vecOfRef = vecFinal.cross(zpositive);
		}
		else if (abs(vecFinal.dot(ypositive)) != vecFinal.norm()*ypositive.norm()) {
			vecOfRef = vecFinal.cross(ypositive);
		}
		else if (abs(vecFinal.dot(xpositive)) != vecFinal.norm()*xpositive.norm()) {
			vecOfRef = vecFinal.cross(xpositive);
		}
		//DEFINEPOINT
		vecOfRef.normalize();
		axis = vecFinal;
		axis.normalize();
		q = AngleAxis<float>(angle, axis);
		point = q * vecOfRef;
		point.normalize();
		result = start + dist * vecFinal + ((radiusFinal + dist * (radiusFinal2 - radiusFinal))*point);
		geom[j] = { result[0], result[1], result[2] };
		geomPerLine[lineIndex][j - startingPOint + startingPOintPerLine] = { result[0], result[1], result[2] };
	}

	return;
}



void treeGeneration::generateCloudForLineSegment(
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
)
{
	int startingPOint;
	int numOfPOints;
	int startingPOintPerLine;
	float radiusFinal;
	float radiusFinal2;
	float dist;
	float angle;
	float overshoot;
	Vector3f zpositive = { 0.0, 0.0, 1.0 };
	Vector3f ypositive = { 0.0, 1.0, 1.0 };
	Vector3f xpositive = { 1.0, 0.0, 0.0 };
	Vector3f vecOfRef;
	Vector3f axis;
	Quaternion<float> q;
	Vector3f point;
	Vector3f result;
	Vector3f vecFinal = end - start;
	//std::cout << "Point Index " << pointIndex << "Line Index " << lineIndex << std::endl;
	startingPOint = geom.size();
	numOfPOints = (int)(density * vecFinal.norm());
	geom.resize(geom.size() + numOfPOints);
	startingPOintPerLine = geomPerLine[lineIndex].size(); //pointIndex
	geomPerLine[lineIndex].resize(startingPOintPerLine + numOfPOints);
#pragma omp parallel shared(vecFinal) private(q, point,result,vecOfRef,angle,overshoot,dist,axis,radiusFinal,radiusFinal2)
#pragma omp for
	for (int j = startingPOint; j < startingPOint + numOfPOints; j++) {
		//RADIUS========
		radiusFinal = radiusMultiplier * this->diameterPerLineSegment[lineIndex];
		radiusFinal2 = radiusMultiplier * this->diameterPerLineSegment[lineIndex];


		//ANGLE=========
		angle = (2 * Pi)*functions.uniformRandom();  //Angle in rad
		//DISTANCE=======
		overshoot = (overshootPercentage*radiusFinal + vecFinal.norm()) / vecFinal.norm();
		if (useUniformDistance) {
			if (lineIndex == 0) {
				//dist = 0.55*overshoot*functions.uniformRandom() - ((overshoot - 1.0) / 2.0) + 0.45;  //0.5 is important
				dist = overshoot * functions.uniformRandom() - ((overshoot - 1.0) / 2.0);
			}
			if (lineIndex > 0) {
				dist = overshoot * functions.uniformRandom() - ((overshoot - 1.0) / 2.0);
			}
		}
		else {
			dist = 0.5 + dispersion * functions.normalRandom();
		}
		//SETAXIS
		if (abs(vecFinal.dot(zpositive)) != vecFinal.norm()*zpositive.norm()) {
			vecOfRef = vecFinal.cross(zpositive);
		}
		else if (abs(vecFinal.dot(ypositive)) != vecFinal.norm()*ypositive.norm()) {
			vecOfRef = vecFinal.cross(ypositive);
		}
		else if (abs(vecFinal.dot(xpositive)) != vecFinal.norm()*xpositive.norm()) {
			vecOfRef = vecFinal.cross(xpositive);
		}
		//DEFINEPOINT
		vecOfRef.normalize();
		axis = vecFinal;
		axis.normalize();
		q = AngleAxis<float>(angle, axis);
		point = q * vecOfRef;
		point.normalize();
		result = start + dist * vecFinal + ((radiusFinal + dist * (radiusFinal2 - radiusFinal))*point);
		geom[j] = { result[0], result[1], result[2] };
		geomPerLine[lineIndex][j - startingPOint + startingPOintPerLine] = { result[0], result[1], result[2] };
	}

	return;
}


void treeGeneration::generateCloudForLineTip(
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
)
{
	int startingPOint;
	int numOfPOints;
	int startingPOintPerLine;
	float radiusFinal;
	float radiusFinal2;
	float dist;
	float angle;
	float overshoot;
	Vector3f zpositive = { 0.0, 0.0, 1.0 };
	Vector3f ypositive = { 0.0, 1.0, 1.0 };
	Vector3f xpositive = { 1.0, 0.0, 0.0 };
	Vector3f vecOfRef;
	Vector3f axis;
	Quaternion<float> q;
	Vector3f point;
	Vector3f result;
	Vector3f vecFinal = end - start;

	startingPOint = geom.size();
	numOfPOints = (int)(density * vecFinal.norm());
	//numOfPOints = (int)(density);

	geom.resize(geom.size() + numOfPOints);
	//GeomPerLineSection
	startingPOintPerLine = geomPerLine[lineIndex].size(); //pointIndex
	//geomPerLine[v[1] - 1].resize(startingPOintPerLine + numOfPOints);
	geomPerLine[lineIndex].resize(startingPOintPerLine + numOfPOints);
#pragma omp parallel shared(vecFinal) private(q, point,result,vecOfRef,angle,overshoot,dist,axis,radiusFinal,radiusFinal2)
#pragma omp for
	for (int j = startingPOint; j < startingPOint + numOfPOints; j++) {
		//RADIUS========
		radiusFinal = (2 * functions.uniformRandom() - 1) *radiusMultiplier*this->diameterPerLineSegment[lineIndex];
		radiusFinal2 = (2 * functions.uniformRandom() - 1) *radiusMultiplier*this->diameterPerLineSegment[lineIndex];
		//ANGLE=========
		angle = (2 * Pi)*functions.uniformRandom();  //Angle in rad
		//DISTANCE=======
		overshoot = (overshootPercentage*radiusFinal + vecFinal.norm()) / vecFinal.norm();
		if (useUniformDistance) {
			dist = 1;
		}
		else {
			dist = 0.5 + dispersion * functions.normalRandom();
		}
		//SETAXIS
		if (abs(vecFinal.dot(zpositive)) != vecFinal.norm()*zpositive.norm()) {
			vecOfRef = vecFinal.cross(zpositive);
		}
		else if (abs(vecFinal.dot(ypositive)) != vecFinal.norm()*ypositive.norm()) {
			vecOfRef = vecFinal.cross(ypositive);
		}
		else if (abs(vecFinal.dot(xpositive)) != vecFinal.norm()*xpositive.norm()) {
			vecOfRef = vecFinal.cross(xpositive);
		}
		//DEFINEPOINT
		vecOfRef.normalize();
		axis = vecFinal;
		axis.normalize();
		q = AngleAxis<float>(angle, axis);
		point = q * vecOfRef;
		point.normalize();
		result = start + dist * vecFinal + (radiusFinal*point);
		geom[j] = { result[0], result[1], result[2] };
		geomPerLine[lineIndex][j - startingPOint + startingPOintPerLine] = { result[0], result[1], result[2] };
	}

	return;
}



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
)
{
	Vector3f vecFinal = end - start;
	int numOfPOints = (int)(density * vecFinal.norm());
	float angle;
	float overshoot;
	float dist;
	Vector3f zpositive = { 0.0, 0.0, 1.0 };
	Vector3f ypositive = { 0.0, 1.0, 1.0 };
	Vector3f xpositive = { 1.0, 0.0, 0.0 };
	Vector3f vecOfRef;
	Vector3f axis;
	Quaternion<float> q;
	Vector3f point;
	Vector3f result;

	geom->vertices.resize(numOfPOints);
	for (int j = 0; j < numOfPOints; j++) {
		
		float radiusFinal = (2 * functions.uniformRandom() - 1) *radiusMultiplier*R2;
		
		//ANGLE=========
		angle = (2 * Pi)*functions.uniformRandom();  //Angle in rad
		//DISTANCE=======
		overshoot = (overshootPercentage*radiusFinal + vecFinal.norm()) / vecFinal.norm();
		if (useUniformDistance) {
			dist = 1;
		}
		else {
			dist = 0.5 + dispersion * functions.normalRandom();
		}
		//SETAXIS
		if (abs(vecFinal.dot(zpositive)) != vecFinal.norm()*zpositive.norm()) {
			vecOfRef = vecFinal.cross(zpositive);
		}
		else if (abs(vecFinal.dot(ypositive)) != vecFinal.norm()*ypositive.norm()) {
			vecOfRef = vecFinal.cross(ypositive);
		}
		else if (abs(vecFinal.dot(xpositive)) != vecFinal.norm()*xpositive.norm()) {
			vecOfRef = vecFinal.cross(xpositive);
		}
		//DEFINEPOINT
		vecOfRef.normalize();
		axis = vecFinal;
		axis.normalize();
		q = AngleAxis<float>(angle, axis);
		point = q * vecOfRef;
		point.normalize();
		result = start + dist * vecFinal + (radiusFinal*point);
		geom->vertices[j] = { result[0], result[1], result[2] };
	}
	return;
}











void treeGeneration::get3DPointCloudv3(int density, float radiusMultiplier) {
	std::vector<std::vector<float>> geom;
	std::vector<std::vector<std::vector<float>>> geomPerLine;
	dotObj solid;
	float dispersion = 0.27;
	float overshootPercentage = 0; //= 0.92;
	bool useUniformDistance = true;

	geomPerLine.resize(oneDimension.lines.size());

	if (oneDimension.lines.size() > 0 && oneDimension.vertices.size() > 0) {
		for (int i = 0; i < oneDimension.lines.size(); i++) {
			int bindex = oneDimension.lines[i][1] - 1;
			int aindex = oneDimension.lines[i][0] - 1;
			Vector3f b = Vector3f(oneDimension.vertices[bindex][0], oneDimension.vertices[bindex][1], oneDimension.vertices[bindex][2]);
			Vector3f a = Vector3f(oneDimension.vertices[aindex][0], oneDimension.vertices[aindex][1], oneDimension.vertices[aindex][2]);
			generateCloudForLineSegment(geom, geomPerLine, a, b, i, density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);
			solid.vertices.insert(std::end(solid.vertices), std::begin(geom), std::end(geom));
		}
	}

	this->treePC = solid;
	this->pointCloudPerLine = geomPerLine;

	return;
}

void treeGeneration::get3DPointCloudv3(bool useHostMesh, dotObj & hostMesh, int numgen, int density, float radiusMultiplier) {
	Vector3f Points[4];
	std::vector<std::vector<float>> geom;
	std::vector<std::vector<std::vector<float>>> geomPerLine;
	dotObj solid;
	Vector3f  vecParent;
	Vector3f  vecLeftChild;
	Vector3f  vecRightChild;
	int numOfPOints = 700;
	int hostMeshSize;
	int genParent = 0;
	int genreport = 0;
	int genLeftChild;
	int genRightChild;
	Vector3f point;
	Vector3f vecOfRef;
	int curgen = 0;
	float dispersion = 0.27;
	float overshootPercentage = 0; //= 0.92;
	bool useUniformDistance = true;
	int newsize = (int)(powf(2.0, ((float)numgen + 1)) - 1);

	//GeomPerLineSection============================

	geomPerLine.resize(oneDimension.lines.size());

	//Hostmesh==================================
	if (useHostMesh) {
		//hostMeshSize = hostMesh.lines.size() - 1; //Check only size
		curgen = getCurrentGeneration(hostMesh.lines.size());
		hostMeshSize = (int)powf(2.0, curgen + 1) - 1;
	}
	else
	{
		hostMeshSize = 0;
	}

	//INITALIZE
	if (oneDimension.lines.size() > 0 && oneDimension.vertices.size() > 0) {
		//for (int i = hostMeshSize; i < oneDimension.lines.size(); i++){
		for (int i = hostMeshSize; i < newsize; i++) {
			int v[4];
			int vindex[3];
			int vindexptr = 0;
			int vptr = 0;
			vindex[vindexptr++] = i;
			v[vptr++] = oneDimension.lines.at(i).at(0) - 1;
			v[vptr++] = oneDimension.lines.at(i).at(1) - 1;
			for (int g = 0; g < oneDimension.lines.size(); g++) {
				if (oneDimension.lines.at(g).at(0) - 1 == v[1]) {
					v[vptr++] = oneDimension.lines.at(g).at(1) - 1;
					vindex[vindexptr++] = g;
				}
			}
			genParent = getCurrentGeneration(v[1]); // Match branch id to generation starting from 0  //v[0]+1
			if (genreport != genParent) {
				genreport = genParent;
				std::cout << "Current generation under process : " << genreport << std::endl;
			}
			genLeftChild = genParent + 1;
			genRightChild = genParent + 1;
			for (int d = 0; d < vptr; d++) {
				Points[d][0] = oneDimension.vertices.at(v[d]).at(0);
				Points[d][1] = oneDimension.vertices.at(v[d]).at(1);
				Points[d][2] = oneDimension.vertices.at(v[d]).at(2);
			}
			vecParent = Points[1] - Points[0];
			if (vptr == 4) {
				vecLeftChild = Points[2] - Points[1];
				vecRightChild = Points[3] - Points[1];
			}

			/*Vector3f uk;
			Vector3f vk;
			uk = Points[1] - Points[0];
			vk = Points[2] - Points[1];
			float theta1 = ((uk).dot(vk)) / ((uk).norm()*(vk).norm());
			uk = Points[1] - Points[0];
			vk = Points[3] - Points[1];
			float theta2 = ((uk).dot(vk)) / ((uk).norm()*(vk).norm());*/

			std::vector<std::vector<float>> geom;

			if (vptr == 4) {
				//generateCloudForLineSegment(geom, geomPerLine, Points[0], Points[1], vindex[0], density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);
				//generateCloudForLineTip(geom, geomPerLine, Points[0], Points[1], vindex[0], density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);
				if (i == 0) {
					generateCloudForLineSegment(geom, geomPerLine, Points[0], Points[1], vindex[0], density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);
				}
				else {
					generateCloudForLineSegment(geom, geomPerLine, Points[0], Points[1], vindex, density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);
				}
				//generateCloudForLineSegment(geom, geomPerLine, i, Points, 0, 1,  v, vindex, density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);
				//generateCloudForLineSegment(geom, geomPerLine, i, Points, 1, 2,  v, vindex, density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);
				//generateCloudForLineSegment(geom, geomPerLine, i, Points, 1, 3,  v, vindex, density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);
				//generateCloudForLineTip(geom, geomPerLine, i, Points, 0, 1, v, vindex, density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);
			}

			if (vptr != 4) {
				generateCloudForLineSegment(geom, geomPerLine, Points[0], Points[1], vindex[0], density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);
				generateCloudForLineTip(geom, geomPerLine, Points[0], Points[1], vindex[0], density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);

				//generateCloudForLineSegment(geom, geomPerLine, i, Points, 0, 1, v, vindex, density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);
				//generateCloudForLineTip(geom, geomPerLine, i, Points, 0, 1,  v, vindex, density, radiusMultiplier, dispersion, overshootPercentage, useUniformDistance);
			}

			solid.vertices.insert(std::end(solid.vertices), std::begin(geom), std::end(geom));
		}
	}
	this->treePC = solid;
	this->pointCloudPerLine = geomPerLine;
	return;
}

void treeGeneration::cleanUpv3(void) {
	std::vector<std::vector<float>> newvectemp;
	Vector3f p0, p1, point, proj;
	bool found;
	float radiusMultiplier = 0.50;
	float limCutOff = 1.05;
	int countup = 0;
	float previousValue = -1.0;

	std::cout << "Tree generation clean up process commencing" << std::endl;
	std::cout << "Before clean up: " << this->treePC.vertices.size() << std::endl;

	for (int ic = 0; ic < this->pointCloudPerLine.size(); ic++) {
		std::vector<int> v;
		v.clear();
		int start = this->oneDimension.lines.at(ic).at(0) - 1;
		int end = this->oneDimension.lines.at(ic).at(1) - 1;
		for (int j = 0; j < this->oneDimension.lines.size(); j++) {
			if (this->oneDimension.lines.at(j).at(0) - 1 == start) {
				v.push_back(j);
			}
			if (this->oneDimension.lines.at(j).at(0) - 1 == end) {
				v.push_back(j);
			}
			if (this->oneDimension.lines.at(j).at(1) - 1 == start) {
				v.push_back(j);
			}
		}

		for (int i = 0; i < this->pointCloudPerLine[ic].size(); i++) {
			point[0] = this->pointCloudPerLine[ic][i][0];
			point[1] = this->pointCloudPerLine[ic][i][1];
			point[2] = this->pointCloudPerLine[ic][i][2];
			bool appendPoint = true;
			float percent = (((float)countup / (float)this->treePC.vertices.size()));
			countup++;
			float percentb = static_cast<int>(percent * 1000) / 1000.0;
			percentb = percentb * 100;
			if (percentb != previousValue) {
				previousValue = percentb;
				std::cout << "Completed : " << percentb << "%" << std::endl;
			}

#pragma omp parallel shared(appendPoint,radiusMultiplier,v) private(p0,p1,proj)
#pragma omp for
			for (int j = 0; j < this->oneDimension.lines.size(); j++) { //for (int j = 0; j < v.size(); j++){
				int k = j; //int k = v[j];
				std::vector<std::vector<float>> * vm = &this->oneDimension.vertices;
				std::vector<std::vector<int>> * lm = &this->oneDimension.lines;
				p0 = Vector3f((*vm)[(*lm)[k][0] - 1][0], (*vm)[(*lm)[k][0] - 1][1], (*vm)[(*lm)[k][0] - 1][2]);
				p1 = Vector3f((*vm)[(*lm)[k][1] - 1][0], (*vm)[(*lm)[k][1] - 1][1], (*vm)[(*lm)[k][1] - 1][2]);
				proj = this->functions.projPointOnVector(point, p0, p1);

				int vk[4];
				int vindex[3];
				int vindexptr = 0;
				int vptr = 0;
				vindex[vindexptr++] = j;
				vk[vptr++] = oneDimension.lines.at(j).at(0) - 1;
				vk[vptr++] = oneDimension.lines.at(j).at(1) - 1;
				for (int g = 0; g < oneDimension.lines.size(); g++) {
					if (oneDimension.lines.at(g).at(0) - 1 == vk[1]) {
						vk[vptr++] = oneDimension.lines.at(g).at(1) - 1;
						vindex[vindexptr++] = g;
					}
				}

				float lim0, lim1;
				if ((vindexptr > 1) && (j > 0)) {
					//int gen = getCurrentGeneration(this->oneDimension.lines.at(k).at(1) - 1);  //j+1  <!---CAUTION
					//float lim0 = radiusMultiplier*this->anatomyBasedRadius.at(gen);
					//float lim1 = radiusMultiplier*this->anatomyBasedRadius.at(gen + 1);
					lim0 = radiusMultiplier * this->diameterPerLineSegment[vindex[0]];
					lim1 = radiusMultiplier * this->diameterPerLineSegment[vindex[1]];
				}
				else {
					//int gen = getCurrentGeneration(this->oneDimension.lines.at(k).at(1) - 1);  //j+1  <!---CAUTION
					//float lim0 = radiusMultiplier*this->anatomyBasedRadius.at(gen);
					//float lim1 = radiusMultiplier*this->anatomyBasedRadius.at(gen + 1);
					lim0 = radiusMultiplier * this->diameterPerLineSegment[j];
					lim1 = radiusMultiplier * this->diameterPerLineSegment[j];
				}

				if (appendPoint) {
					if ((p0 - proj).norm() + (proj - p1).norm() <= limCutOff * (p1 - p0).norm()) {
						if (ic == 0) {
							lim1 = lim0;
						}
						float dist = ((p0 - proj).norm() / (p0 - p1).norm());
						//float lim = radiusMultiplier*this->diameters.at(j);
						float pointToProj = (proj - point).norm();
						if (pointToProj < lim0 + dist * (lim1 - lim0)) {
							//toBeDeleted.push_back(i);
							//toBeDeleted[i] = 1;
							appendPoint = false;
							//std::cout << "Deleted point" << i << std::endl;
						}
						//if ((proj - p1).norm() < 0.05){ appendPoint = true; }
						//if ((proj - p0).norm() < 0.05){ appendPoint = false; }
					}
					if ((v.size() == 3) && ((proj - p1).norm() < 0.1*(p0 - p1).norm())) { appendPoint = true; }
				}
			}
			if (appendPoint) {
#pragma omp critical
				newvectemp.push_back(this->pointCloudPerLine[ic][i]);
			}
		}
	}
	if (newvectemp.size() > 0) {
		this->treePC.vertices = newvectemp;
	}
	std::cout << "After clean up: " << this->treePC.vertices.size() << std::endl;
	return;
}

void treeGeneration::cleanUpv3(float radiusMultiplier) {
	std::vector<std::vector<float>> newvectemp;
	Vector3f p0, p1, point, proj;
	bool found;
	float limCutOff = 1.05;
	int countup = 0;
	float previousValue = -1.0;

	std::cout << "Tree generation clean up process commencing" << std::endl;
	std::cout << "Before clean up: " << this->treePC.vertices.size() << std::endl;

	for (int ic = 0; ic < this->pointCloudPerLine.size(); ic++) {
		std::vector<int> v;
		v.clear();
		int start = this->oneDimension.lines.at(ic).at(0) - 1;
		int end = this->oneDimension.lines.at(ic).at(1) - 1;
		for (int j = 0; j < this->oneDimension.lines.size(); j++) {
			if (this->oneDimension.lines.at(j).at(0) - 1 == start) {
				v.push_back(j);
			}
			if (this->oneDimension.lines.at(j).at(0) - 1 == end) {
				v.push_back(j);
			}
			if (this->oneDimension.lines.at(j).at(1) - 1 == start) {
				v.push_back(j);
			}
		}

		for (int i = 0; i < this->pointCloudPerLine[ic].size(); i++) {
			point[0] = this->pointCloudPerLine[ic][i][0];
			point[1] = this->pointCloudPerLine[ic][i][1];
			point[2] = this->pointCloudPerLine[ic][i][2];
			bool appendPoint = true;
			float percent = (((float)countup / (float)this->treePC.vertices.size()));
			countup++;
			float percentb = static_cast<int>(percent * 1000) / 1000.0;
			percentb = percentb * 100;
			if (percentb != previousValue) {
				previousValue = percentb;
				std::cout << "Completed : " << percentb << "%" << std::endl;
			}

#pragma omp parallel shared(appendPoint,radiusMultiplier,v) private(p0,p1,proj)
#pragma omp for
			for (int j = 0; j < this->oneDimension.lines.size(); j++) { //for (int j = 0; j < v.size(); j++){
				int k = j; //int k = v[j];
				std::vector<std::vector<float>> * vm = &this->oneDimension.vertices;
				std::vector<std::vector<int>> * lm = &this->oneDimension.lines;
				p0 = Vector3f((*vm)[(*lm)[k][0] - 1][0], (*vm)[(*lm)[k][0] - 1][1], (*vm)[(*lm)[k][0] - 1][2]);
				p1 = Vector3f((*vm)[(*lm)[k][1] - 1][0], (*vm)[(*lm)[k][1] - 1][1], (*vm)[(*lm)[k][1] - 1][2]);
				proj = this->functions.projPointOnVector(point, p0, p1);

				int vk[4];
				int vindex[3];
				int vindexptr = 0;
				int vptr = 0;
				vindex[vindexptr++] = j;
				vk[vptr++] = oneDimension.lines.at(j).at(0) - 1;
				vk[vptr++] = oneDimension.lines.at(j).at(1) - 1;
				for (int g = 0; g < oneDimension.lines.size(); g++) {
					if (oneDimension.lines.at(g).at(0) - 1 == vk[1]) {
						vk[vptr++] = oneDimension.lines.at(g).at(1) - 1;
						vindex[vindexptr++] = g;
					}
				}

				float lim0, lim1;
				if ((vindexptr > 1) && (j > 0)) {
					//int gen = getCurrentGeneration(this->oneDimension.lines.at(k).at(1) - 1);  //j+1  <!---CAUTION
					//float lim0 = radiusMultiplier*this->anatomyBasedRadius.at(gen);
					//float lim1 = radiusMultiplier*this->anatomyBasedRadius.at(gen + 1);
					lim0 = radiusMultiplier * this->diameterPerLineSegment[vindex[0]];
					lim1 = radiusMultiplier * this->diameterPerLineSegment[vindex[1]];
				}
				else {
					//int gen = getCurrentGeneration(this->oneDimension.lines.at(k).at(1) - 1);  //j+1  <!---CAUTION
					//float lim0 = radiusMultiplier*this->anatomyBasedRadius.at(gen);
					//float lim1 = radiusMultiplier*this->anatomyBasedRadius.at(gen + 1);
					lim0 = radiusMultiplier * this->diameterPerLineSegment[j];
					lim1 = radiusMultiplier * this->diameterPerLineSegment[j];
				}

				if (appendPoint) {
					if ((p0 - proj).norm() + (proj - p1).norm() <= limCutOff * (p1 - p0).norm()) {
						if (ic == 0) {
							lim1 = lim0;
						}
						float dist = ((p0 - proj).norm() / (p0 - p1).norm());
						//float lim = radiusMultiplier*this->diameters.at(j);
						float pointToProj = (proj - point).norm();
						if (pointToProj < lim0 + dist * (lim1 - lim0)) {
							//toBeDeleted.push_back(i);
							//toBeDeleted[i] = 1;
							appendPoint = false;
							//std::cout << "Deleted point" << i << std::endl;
						}
						//if ((proj - p1).norm() < 0.05){ appendPoint = true; }
						//if ((proj - p0).norm() < 0.05){ appendPoint = false; }
					}
					if ((v.size() == 3) && ((proj - p1).norm() < 0.1*(p0 - p1).norm())) { appendPoint = true; }
				}
			}
			if (appendPoint) {
#pragma omp critical
				newvectemp.push_back(this->pointCloudPerLine[ic][i]);
			}
		}
	}
	if (newvectemp.size() > 0) {
		this->treePC.vertices = newvectemp;
	}
	std::cout << "After clean up: " << this->treePC.vertices.size() << std::endl;
	return;
}

void treeGeneration::flatten(std::string result, std::string  keyword) {
	Vector3f p, projection;
	int gen = 0;
	this->anatomyBasedRadius;
	this->oneDimension;
	this->treePC;
	this->reconstructedModel;
	int v0, v1, v2, v3;
	bool isTerminal = false;
	int pp = 0;
	Vector3f Points[4];
	int sidenum = 0;
	if (keyword == "L") { sidenum = 1; }
	if (keyword == "R") { sidenum = 2; }
	float completePercentage;

	for (int i = 0; i < this->oneDimension.lines.size(); i++) {
		v0 = this->oneDimension.lines.at(i).at(0) - 1;
		v1 = this->oneDimension.lines.at(i).at(1) - 1;

		/*for (int ik = 0; ik < this->oneDimension.lines.size(); ik++){
		}*/

		completePercentage = (100.0*(float)i) / this->oneDimension.lines.size();
		completePercentage = (int)(completePercentage * 10);

		if ((int)completePercentage % 10 == 0) {
			std::cout << "Complete: " << completePercentage / 10 << "%" << std::endl;
		}
		gen = this->getCurrentGeneration(v1);
		pp = 0;
		Points[0][0] = this->oneDimension.vertices.at(v0).at(0);
		Points[0][1] = this->oneDimension.vertices.at(v0).at(1);
		Points[0][2] = this->oneDimension.vertices.at(v0).at(2);

		pp++;
		Points[1][0] = this->oneDimension.vertices.at(v1).at(0);
		Points[1][1] = this->oneDimension.vertices.at(v1).at(1);
		Points[1][2] = this->oneDimension.vertices.at(v1).at(2);

		isTerminal = true;
		for (int j = 0; j < this->oneDimension.lines.size(); j++) {
			v2 = this->oneDimension.lines.at(j).at(0) - 1;

			if (v2 == v1) {
				isTerminal = false;
				pp++;
				Points[pp][0] = this->oneDimension.vertices.at(v1).at(0);
				Points[pp][1] = this->oneDimension.vertices.at(v1).at(0);
				Points[pp][2] = this->oneDimension.vertices.at(v1).at(0);
			}
		}

		if (isTerminal) {
			this->outletCounter++;
			std::vector<int> verticesOnOutlet;
			for (int j = 0; j < this->reconstructedModel.vertices.size(); j++) {
				p[0] = this->reconstructedModel.vertices.at(j).at(0);
				p[1] = this->reconstructedModel.vertices.at(j).at(1);
				p[2] = this->reconstructedModel.vertices.at(j).at(2);
				//projection = this->functions.projPointOnVector(p, p0, p1);

				if ((p - Points[1]).norm() <= this->anatomyBasedRadius.at(gen)) { //<------Caution Flattening boundaries
					//Calculate face normal
					Vector3f pnormal;
					Vector3f direction = Points[1] - Points[0];
					bool found = false;
					bool doFlat = false;
					for (int jf = 0; jf < this->reconstructedModel.faces.size(); jf++) {
						if (!found) {
							if ((this->reconstructedModel.faces.at(jf).at(0) - 1 == j) ||
								(this->reconstructedModel.faces.at(jf).at(3) - 1 == j) ||
								(this->reconstructedModel.faces.at(jf).at(6) - 1 == j)) {
								int vertex0 = this->reconstructedModel.faces.at(jf).at(0) - 1;
								int vertex1 = this->reconstructedModel.faces.at(jf).at(3) - 1;
								int vertex2 = this->reconstructedModel.faces.at(jf).at(6) - 1;
								Vector3f v0, v1, v2;
								v0[0] = this->reconstructedModel.vertices.at(vertex0).at(0);
								v0[1] = this->reconstructedModel.vertices.at(vertex0).at(1);
								v0[2] = this->reconstructedModel.vertices.at(vertex0).at(2);
								v1[0] = this->reconstructedModel.vertices.at(vertex1).at(0);
								v1[1] = this->reconstructedModel.vertices.at(vertex1).at(1);
								v1[2] = this->reconstructedModel.vertices.at(vertex1).at(2);
								v2[0] = this->reconstructedModel.vertices.at(vertex2).at(0);
								v2[1] = this->reconstructedModel.vertices.at(vertex2).at(1);
								v2[2] = this->reconstructedModel.vertices.at(vertex2).at(2);
								Vector3f edge1 = v1 - v0;
								edge1.normalize();
								Vector3f edge2 = v1 - v2;
								edge2.normalize();
								pnormal = edge1.cross(edge2);
								pnormal.normalize();
								float angleBetweenCos = (direction.dot(pnormal)) / (direction.norm()*pnormal.norm());
								found = true;
								if (abs(angleBetweenCos) > 0.1) { doFlat = true; } //Caution flattening limit
							}
						}
					}
					if (doFlat) {
						verticesOnOutlet.push_back(j);
						Vector3f n = (Points[1] - Points[0]);
						n.normalize();
						Vector3f pf = Points[0] + 1.0*(Points[1] - Points[0]);				 //<------Caution Flattening boundaries
						Vector3f r = p - ((p - pf).dot(n))*n;				//pf<--p1 //<------Caution Flattening boundaries
						//p = r;
						this->reconstructedModel.vertices.at(j).at(0) = r[0];
						this->reconstructedModel.vertices.at(j).at(1) = r[1];
						this->reconstructedModel.vertices.at(j).at(2) = r[2];
						//this->reconstructedModel.segment_property_map_per_vertex[j] = 99 * EXPMAXNUMOFBRANCH + outletCounter; ///<----Caution
					}
					/*if (1){
					Vector3f n = (Points[1] - Points[0]);
					n.normalize();
					Vector3f pf = Points[1];
					Vector3f r = p - ((p - pf).dot(n))*n;				//pf<--p1 //<------Caution Flattening boundaries
					//p = r;
					this->reconstructedModel.vertices.at(j).at(0) = r[0];
					this->reconstructedModel.vertices.at(j).at(1) = r[1];
					this->reconstructedModel.vertices.at(j).at(2) = r[2];
					}*/
					/*
					//this->reconstructedModel.segment_property_map_per_vertex[j] = 99 * EXPMAXNUMOFBRANCH + outletCounter;// v1 - (powf(2.0, (float)gen) - 1); ///<----Caution
					*/
				}
			}

			//Do something with verticesOnOutlet
		}

		if (false) {
			for (int j = 0; j < this->reconstructedModel.vertices.size(); j++) {
				p[0] = this->reconstructedModel.vertices.at(j).at(0);
				p[1] = this->reconstructedModel.vertices.at(j).at(1);
				p[2] = this->reconstructedModel.vertices.at(j).at(2);

				projection = this->functions.projPointOnVector(p, Points[0], Points[1]);

				if ((projection - Points[0]).norm() + (projection - Points[1]).norm() <= 1.05*(Points[0] - Points[1]).norm()) {
					if ((projection - p).norm() <= 0.55*this->anatomyBasedRadius.at(gen)) {
						//Assign to generation
						if (this->reconstructedModel.segment_property_map_per_vertex[j] < 0) {
							this->reconstructedModel.segment_property_map_per_vertex[j] = gen * EXPMAXNUMOFBRANCH + (sidenum * EXPMAXNUMOFBRANCH / 10) + v0 - (powf(2.0, (float)gen) - 1); ///<----Caution
						}
					}
				}
			}
		}
	}

	return;
}

void treeGeneration::naming(std::string result, std::string  keyword) {
	Vector3f p, p0, p1, projection;
	int gen = 0;
	this->anatomyBasedRadius;
	this->oneDimension;
	this->treePC;
	this->reconstructedModel;
	int v0, v1, v2, v3;
	bool isTerminal = false;

	this->reconstructedModel.segment_property_map_per_vertex.resize(this->reconstructedModel.vertices.size());

	for (int i = 0; i < this->reconstructedModel.segment_property_map_per_vertex.size(); i++) {
		this->reconstructedModel.segment_property_map_per_vertex[i] = -5;
	}

	for (int i = 0; i < this->oneDimension.lines.size(); i++) {
		v0 = this->oneDimension.lines.at(i).at(0) - 1;
		v1 = this->oneDimension.lines.at(i).at(1) - 1;

		gen = this->getCurrentGeneration(v1);

		p0[0] = this->oneDimension.vertices.at(v0).at(0);
		p0[1] = this->oneDimension.vertices.at(v0).at(1);
		p0[2] = this->oneDimension.vertices.at(v0).at(2);

		p1[0] = this->oneDimension.vertices.at(v1).at(0);
		p1[1] = this->oneDimension.vertices.at(v1).at(1);
		p1[2] = this->oneDimension.vertices.at(v1).at(2);

		//p1 = p0 + 0.9*(p1 - p0);

		isTerminal = true;

		for (int j = 0; j < this->oneDimension.lines.size(); j++) {
			v2 = this->oneDimension.lines.at(j).at(0) - 1;

			if (v2 == v1) { isTerminal = false; }
		}

		if (isTerminal) {
			for (int j = 0; j < this->reconstructedModel.vertices.size(); j++) {
				p[0] = this->reconstructedModel.vertices.at(j).at(0);
				p[1] = this->reconstructedModel.vertices.at(j).at(1);
				p[2] = this->reconstructedModel.vertices.at(j).at(2);

				if ((p - p1).norm() < 0.52*this->anatomyBasedRadius.at(gen)) {  //<------Caution
					projection = this->functions.projPointOnVector(p, p0, p1);
					if ((p - projection).norm() < 0.53*this->anatomyBasedRadius.at(gen)) { //<------Caution Flattening boundaries
						if ((projection - p0).norm() > (p1 - p0).norm()) {		//<------Caution Flattening boundaries
							if ((projection - p1).norm() < (p1 - p0).norm()) { //<------Caution Flattening boundaries
								//project on plane vertical to p0-p1

								Vector3f n = (p1 - p0);
								n.normalize();
								Vector3f pf = p0 + 1.0*(p1 - p0);				 //<------Caution Flattening boundaries
								Vector3f r = p - ((p - pf).dot(n))*n;				//pf<--p1 //<------Caution Flattening boundaries
								//p = r;
								this->reconstructedModel.vertices.at(j).at(0) = r[0];
								this->reconstructedModel.vertices.at(j).at(1) = r[1];
								this->reconstructedModel.vertices.at(j).at(2) = r[2];
								this->reconstructedModel.segment_property_map_per_vertex[j] = 99 * EXPMAXNUMOFBRANCH + v1 - (powf(2.0, (float)gen) - 1); ///<----Caution
							}
						}
					}
				}
			}
		}

		//#pragma omp parallel  private(p,projection,p1,p0,gen)
		//#pragma omp for
		for (int j = 0; j < this->reconstructedModel.vertices.size(); j++) {
			p[0] = this->reconstructedModel.vertices.at(j).at(0);
			p[1] = this->reconstructedModel.vertices.at(j).at(1);
			p[2] = this->reconstructedModel.vertices.at(j).at(2);

			projection = this->functions.projPointOnVector(p, p0, p1);

			if ((projection - p0).norm() + (projection - p1).norm() <= 1.05*(p0 - p1).norm()) {
				if ((projection - p).norm() <= 1.05*this->anatomyBasedRadius.at(gen)) {
					//Assign to generation
					if (this->reconstructedModel.segment_property_map_per_vertex[j] < 0) {
						this->reconstructedModel.segment_property_map_per_vertex[j] = gen * EXPMAXNUMOFBRANCH + v1 - (powf(2.0, (float)gen) - 1); ///<----Caution
					}
				}
			}
		}
	}

	this->reconstructedModel.segment_property_map.resize(this->reconstructedModel.faces.size());
	for (int i = 0; i < this->reconstructedModel.faces.size(); i++) {
		this->reconstructedModel.segment_property_map[i] = 1000;
	}
	for (int i = 0; i < this->reconstructedModel.faces.size(); i++) {
		int vm1 = this->reconstructedModel.faces.at(i).at(0) - 1;
		int vm2 = this->reconstructedModel.faces.at(i).at(3) - 1;
		int vm3 = this->reconstructedModel.faces.at(i).at(6) - 1;

		int group1 = this->reconstructedModel.segment_property_map_per_vertex.at(vm1);
		int group2 = this->reconstructedModel.segment_property_map_per_vertex.at(vm2);
		int group3 = this->reconstructedModel.segment_property_map_per_vertex.at(vm3);

		if ((group1 == group2) && (group3 == group2)) { this->reconstructedModel.segment_property_map.at(i) = group1; }
		else {
			if ((group1 != 99) && (group3 != 99) && (group2 != 99)) {
				if (group1 == group2) { this->reconstructedModel.segment_property_map.at(i) = group1; }

				if (group2 == group3) { this->reconstructedModel.segment_property_map.at(i) = group2; }

				if (group3 == group1) { this->reconstructedModel.segment_property_map.at(i) = group3; }

				if ((group1 != group2) && (group2 != group3)) { this->reconstructedModel.segment_property_map.at(i) = group1; }
			}
			/*
			else{
			if ((group1 == 99) && (group3 == 99)){ this->reconstructedModel.segment_property_map.at(i) = group2; }
			if ((group2 == 99) && (group3 == 99)){ this->reconstructedModel.segment_property_map.at(i) = group1; }
			if ((group1 == 99) && (group2 == 99)){ this->reconstructedModel.segment_property_map.at(i) = group3; }
			}	*/
		}
	}

	this->reconstructedModel.exportToFile(result + keyword + "no");
	//this->reconstructedModel.exportToFileSegmented(result, keyword);
	return;
}

void treeGeneration::naming2(std::string result, std::string  keyword) {
	Vector3f p, projection;
	int gen = 0;
	this->anatomyBasedRadius;
	this->oneDimension;
	this->treePC;
	this->reconstructedModel;
	int v0, v1, v2, v3;
	bool isTerminal = false;
	int pp = 0;
	Vector3f Points[4];
	int sidenum = 0;
	if (keyword == "L") { sidenum = 1; }
	if (keyword == "R") { sidenum = 2; }
	float completePercentage;

	for (int i = 0; i < this->oneDimension.lines.size(); i++) {
		v0 = this->oneDimension.lines.at(i).at(0) - 1;
		v1 = this->oneDimension.lines.at(i).at(1) - 1;

		/*for (int ik = 0; ik < this->oneDimension.lines.size(); ik++){
		}*/

		completePercentage = (100.0*(float)i) / this->oneDimension.lines.size();
		completePercentage = (int)(completePercentage * 10);

		if ((int)completePercentage % 10 == 0) {
			std::cout << "Complete: " << completePercentage / 10 << "%" << std::endl;
		}
		gen = this->getCurrentGeneration(v1);
		pp = 0;
		Points[0][0] = this->oneDimension.vertices.at(v0).at(0);
		Points[0][1] = this->oneDimension.vertices.at(v0).at(1);
		Points[0][2] = this->oneDimension.vertices.at(v0).at(2);

		pp++;
		Points[1][0] = this->oneDimension.vertices.at(v1).at(0);
		Points[1][1] = this->oneDimension.vertices.at(v1).at(1);
		Points[1][2] = this->oneDimension.vertices.at(v1).at(2);

		isTerminal = true;
		for (int j = 0; j < this->oneDimension.lines.size(); j++) {
			v2 = this->oneDimension.lines.at(j).at(0) - 1;

			if (v2 == v1) {
				isTerminal = false;
				pp++;
				Points[pp][0] = this->oneDimension.vertices.at(v1).at(0);
				Points[pp][1] = this->oneDimension.vertices.at(v1).at(0);
				Points[pp][2] = this->oneDimension.vertices.at(v1).at(0);
			}
		}

		if (isTerminal) {
			this->outletCounter++;
			for (int j = 0; j < this->reconstructedModel.vertices.size(); j++) {
				p[0] = this->reconstructedModel.vertices.at(j).at(0);
				p[1] = this->reconstructedModel.vertices.at(j).at(1);
				p[2] = this->reconstructedModel.vertices.at(j).at(2);
				//projection = this->functions.projPointOnVector(p, p0, p1);

				if ((p - Points[1]).norm() <= this->anatomyBasedRadius.at(gen)) { //<------Caution Flattening boundaries
					//Calculate face normal
					Vector3f pnormal;
					Vector3f direction = Points[1] - Points[0];
					bool found = false;
					bool doFlat = false;
					for (int jf = 0; jf < this->reconstructedModel.faces.size(); jf++) {
						if (!found) {
							if ((this->reconstructedModel.faces.at(jf).at(0) - 1 == j) ||
								(this->reconstructedModel.faces.at(jf).at(3) - 1 == j) ||
								(this->reconstructedModel.faces.at(jf).at(6) - 1 == j)) {
								int vertex0 = this->reconstructedModel.faces.at(jf).at(0) - 1;
								int vertex1 = this->reconstructedModel.faces.at(jf).at(3) - 1;
								int vertex2 = this->reconstructedModel.faces.at(jf).at(6) - 1;
								Vector3f v0p, v1p, v2p;
								v0p[0] = this->reconstructedModel.vertices.at(vertex0).at(0);
								v0p[1] = this->reconstructedModel.vertices.at(vertex0).at(1);
								v0p[2] = this->reconstructedModel.vertices.at(vertex0).at(2);
								v1p[0] = this->reconstructedModel.vertices.at(vertex1).at(0);
								v1p[1] = this->reconstructedModel.vertices.at(vertex1).at(1);
								v1p[2] = this->reconstructedModel.vertices.at(vertex1).at(2);
								v2p[0] = this->reconstructedModel.vertices.at(vertex2).at(0);
								v2p[1] = this->reconstructedModel.vertices.at(vertex2).at(1);
								v2p[2] = this->reconstructedModel.vertices.at(vertex2).at(2);
								Vector3f edge1 = v1p - v0p;
								edge1.normalize();
								Vector3f edge2 = v1p - v2p;
								edge2.normalize();
								pnormal = edge1.cross(edge2);
								pnormal.normalize();
								float angleBetweenCos = (direction.dot(pnormal)) / (direction.norm()*pnormal.norm());
								found = true;
								if (abs(angleBetweenCos) > 0.8) { doFlat = true; } //Caution flattening limit
							}
						}
					}
					if (doFlat) {
						/*Vector3f n = (Points[1] - Points[0]);
						n.normalize();
						Vector3f pf = Points[0] + 1.0*(Points[1] - Points[0]);				 //<------Caution Flattening boundaries
						Vector3f r = p - ((p - pf).dot(n))*n;				//pf<--p1 //<------Caution Flattening boundaries
						//p = r;
						this->reconstructedModel.vertices.at(j).at(0) = r[0];
						this->reconstructedModel.vertices.at(j).at(1) = r[1];
						this->reconstructedModel.vertices.at(j).at(2) = r[2];*/
						this->reconstructedModel.segment_property_map_per_vertex[j] = 99 * EXPMAXNUMOFBRANCH + outletCounter; ///<----Caution
					}
				}
			}
		}

		if (true) {
			for (int j = 0; j < this->reconstructedModel.vertices.size(); j++) {
				p[0] = this->reconstructedModel.vertices.at(j).at(0);
				p[1] = this->reconstructedModel.vertices.at(j).at(1);
				p[2] = this->reconstructedModel.vertices.at(j).at(2);

				projection = this->functions.projPointOnVector(p, Points[0], Points[1]);

				if ((projection - Points[0]).norm() + (projection - Points[1]).norm() <= 1.05*(Points[0] - Points[1]).norm()) {
					if ((projection - p).norm() <= 0.59*this->anatomyBasedRadius.at(gen)) {
						//Assign to generation
						if (this->reconstructedModel.segment_property_map_per_vertex[j] < 0) {
							this->reconstructedModel.segment_property_map_per_vertex[j] = gen * EXPMAXNUMOFBRANCH + (sidenum * EXPMAXNUMOFBRANCH / 10) + (EXPMAXNUMOFBRANCH - (v0 - (powf(2.0, (float)gen) - 1 + 1))); ///<----Caution //(sidenum * EXPMAXNUMOFBRANCH / 10) //-( v1 - (powf(2.0, (float)gen-1) - 1));
						}
					}
				}
			}
		}
	}

	return;
}

/*
void treeGeneration::generateCloudForLineTip(
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
)
{
int startingPOint;
int numOfPOints;
int startingPOintPerLine;
float radiusFinal;
float radiusFinal2;
float dist;
float angle;
float overshoot;
Vector3f zpositive = { 0.0, 0.0, 1.0 };
Vector3f ypositive = { 0.0, 1.0, 1.0 };
Vector3f xpositive = { 1.0, 0.0, 0.0 };
Vector3f vecOfRef;
Vector3f axis;
Quaternion<float> q;
Vector3f point;
Vector3f result;
Vector3f vecFinal = Points[pointB] - Points[pointA];

startingPOint = geom.size();
numOfPOints = (int)(density * (Points[pointA] - Points[pointB]).norm());
//numOfPOints = (int)(density);

geom.resize(geom.size() + numOfPOints);
//GeomPerLineSection
startingPOintPerLine = geomPerLine[v[pointB] - 1].size();
//geomPerLine[v[1] - 1].resize(startingPOintPerLine + numOfPOints);
geomPerLine[vindex[pointB - 1]].resize(startingPOintPerLine + numOfPOints);
#pragma omp parallel shared(vecFinal) private(q, point,result,vecOfRef,angle,overshoot,dist,axis,radiusFinal,radiusFinal2)
#pragma omp for
for (int j = startingPOint; j < startingPOint + numOfPOints; j++){
//RADIUS========
//radiusFinal = radiusMultiplier*this->anatomyBasedRadius.at(genParent);
//radiusFinal2 = radiusMultiplier*this->anatomyBasedRadius.at(genParent + 1);
radiusFinal = (2 * functions.uniformRandom() - 1) *radiusMultiplier*this->diameterPerLineSegment[vindex[pointB - 1]];
radiusFinal2 = (2 * functions.uniformRandom() - 1) *radiusMultiplier*this->diameterPerLineSegment[vindex[pointB - 1]];
//ANGLE=========
angle = (2 * Pi)*functions.uniformRandom();  //Angle in rad
//DISTANCE=======
overshoot = (overshootPercentage*radiusFinal + (Points[pointA] - Points[pointB]).norm()) / (Points[pointA] - Points[pointB]).norm();
if (useUniformDistance){
dist = 1;
}
else{
dist = 0.5 + dispersion*functions.normalRandom();
}
//SETAXIS
if (abs(vecFinal.dot(zpositive)) != vecFinal.norm()*zpositive.norm()){
vecOfRef = vecFinal.cross(zpositive);
}
else if (abs(vecFinal.dot(ypositive)) != vecFinal.norm()*ypositive.norm()){
vecOfRef = vecFinal.cross(ypositive);
}
else if (abs(vecFinal.dot(xpositive)) != vecFinal.norm()*xpositive.norm()){
vecOfRef = vecFinal.cross(xpositive);
}
//DEFINEPOINT
vecOfRef.normalize();
axis = vecFinal;
axis.normalize();
q = AngleAxis<float>(angle, axis);
point = q*vecOfRef;
point.normalize();
result = Points[pointA] + dist*vecFinal + (radiusFinal*point);
//result = Points[pointA] + dist*vecFinal + ((radiusFinal + dist*(radiusFinal2 - radiusFinal))*point);
geom[j] = { result[0], result[1], result[2] };
geomPerLine[vindex[pointB - 1]][j - startingPOint + startingPOintPerLine] = { result[0], result[1], result[2] };
}

return;
}

*/

/*
void treeGeneration::generateCloudForLineSegment(
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
)
{
int startingPOint;
int numOfPOints;
int startingPOintPerLine;
float radiusFinal;
float radiusFinal2;
float dist;
float angle;
float overshoot;
Vector3f zpositive = { 0.0, 0.0, 1.0 };
Vector3f ypositive = { 0.0, 1.0, 1.0 };
Vector3f xpositive = { 1.0, 0.0, 0.0 };
Vector3f vecOfRef;
Vector3f axis;
Quaternion<float> q;
Vector3f point;
Vector3f result;
Vector3f vecFinal = Points[pointB] - Points[pointA];

startingPOint = geom.size();
numOfPOints = (int)(density * (Points[pointA] - Points[pointB]).norm());
//numOfPOints = (int)(density);

geom.resize(geom.size() + numOfPOints);
//GeomPerLineSection
startingPOintPerLine = geomPerLine[v[pointB] - 1].size();
//geomPerLine[v[1] - 1].resize(startingPOintPerLine + numOfPOints);
geomPerLine[vindex[pointB - 1]].resize(startingPOintPerLine + numOfPOints);
#pragma omp parallel shared(vecFinal) private(q, point,result,vecOfRef,angle,overshoot,dist,axis,radiusFinal,radiusFinal2)
#pragma omp for
for (int j = startingPOint; j < startingPOint + numOfPOints; j++){
//RADIUS========
//radiusFinal = radiusMultiplier*this->anatomyBasedRadius.at(genParent);
//radiusFinal2 = radiusMultiplier*this->anatomyBasedRadius.at(genParent + 1);
radiusFinal = radiusMultiplier*this->diameterPerLineSegment[vindex[pointB - 1]];
radiusFinal2 = radiusMultiplier*this->diameterPerLineSegment[vindex[pointB - 1]];
if (i == 0){
radiusFinal2 = radiusFinal;
}
//ANGLE=========
angle = (2 * Pi)*functions.uniformRandom();  //Angle in rad
//DISTANCE=======
overshoot = (overshootPercentage*radiusFinal + (Points[pointA] - Points[pointB]).norm()) / (Points[pointA] - Points[pointB]).norm();
if (useUniformDistance){
if (i == 0){
dist = 0.55*overshoot*functions.uniformRandom() - ((overshoot - 1.0) / 2.0) + 0.45;  //0.5 is important
}
if (i > 0){
dist = overshoot*functions.uniformRandom() - ((overshoot - 1.0) / 2.0);
}
}
else{
dist = 0.5 + dispersion*functions.normalRandom();
}
//SETAXIS
if (abs(vecFinal.dot(zpositive)) != vecFinal.norm()*zpositive.norm()){
vecOfRef = vecFinal.cross(zpositive);
}
else if (abs(vecFinal.dot(ypositive)) != vecFinal.norm()*ypositive.norm()){
vecOfRef = vecFinal.cross(ypositive);
}
else if (abs(vecFinal.dot(xpositive)) != vecFinal.norm()*xpositive.norm()){
vecOfRef = vecFinal.cross(xpositive);
}
//DEFINEPOINT
vecOfRef.normalize();
axis = vecFinal;
axis.normalize();
q = AngleAxis<float>(angle, axis);
point = q*vecOfRef;
point.normalize();
result = Points[pointA] + dist*vecFinal + ((radiusFinal + dist*(radiusFinal2 - radiusFinal))*point);
geom[j] = { result[0], result[1], result[2] };
geomPerLine[vindex[pointB - 1]][j - startingPOint + startingPOintPerLine] = { result[0], result[1], result[2] };
}

return;
}

*/