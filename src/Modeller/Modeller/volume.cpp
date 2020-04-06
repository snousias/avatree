#include "lungModelling.h"

void volume::initializeFromFile(std::string thefile) {
	dotObj mdfile;
	mdfile.initializeFromFile(thefile);

	this->vertices = mdfile.vertices;

	return;
}

void volume::createVolumeFromObj(dotObj boundaries, int uniformPerEdge, int maxPoints) {
	int perEdge = uniformPerEdge;
	float x;
	float y;
	float z;
	bbox limits;
	bool pointInPolygonCleanup = false;
	bool laplacianCleanUp = false;
	helper *fns;
	fns = new helper();
	float distanceInitial = 0, distanceFinal = 0;
	std::vector<int> outs;
	Vector3f p1, p2, p3;
	boundaries.getBoundingBox();
	limits = boundaries.boundingBox;
	if (maxPoints == 0) {
		maxPoints = perEdge * perEdge*perEdge;
	}

#pragma omp parallel shared(perEdge) private(x,y,z)
#pragma omp for
	for (int i = 0; i < perEdge; i++)
	{
		for (int j = 0; j < perEdge; j++)
		{
			for (int k = 0; k < perEdge; k++)
			{
				x = 2 * (1.0 / 6.0) * (boundaries.boundingBox.maxX - boundaries.boundingBox.minX)*fns->normalRandom() + ((boundaries.boundingBox.maxX + boundaries.boundingBox.minX) / 2);
				y = 2 * (1.0 / 6.0) * (boundaries.boundingBox.maxY - boundaries.boundingBox.minY)*fns->normalRandom() + ((boundaries.boundingBox.maxY + boundaries.boundingBox.minY) / 2);
				z = 2 * (1.0 / 6.0) * (boundaries.boundingBox.maxZ - boundaries.boundingBox.minZ)*fns->normalRandom() + ((boundaries.boundingBox.maxZ + boundaries.boundingBox.minZ) / 2);

				//x = boundaries.boundingBox.minX + i*(boundaries.boundingBox.maxX - boundaries.boundingBox.minX) / perEdge;
				//y = boundaries.boundingBox.minY + j*(boundaries.boundingBox.maxY - boundaries.boundingBox.minY) / perEdge;
				//z = boundaries.boundingBox.minZ + k*(boundaries.boundingBox.maxZ - boundaries.boundingBox.minZ) / perEdge;
				//boundaries.checkIfPointisInside({ x, y, z })
				if (this->vertices.size() < maxPoints) {
					if ((boundaries.checkIfPointisInsidev2(Vector3f(x, y, z)))) {
						#pragma omp critical
						this->vertices.push_back({ x, y, z });
					}
				}
			}
		}
	}

	return;
}

bool dotObj::checkIfPointisInside(std::vector<float> point) {
	bool isInside = true;
	int vertex0, vertex1, vertex2;
	int normal0, normal1, normal2;
	Vector3f p0, p1, v0, v1, v2, n0, n1, n2;
	Vector3f n, w, u, v, r, pi;
	double pri, s, t;

	std::vector<Vector3f> rays;

	//for (int k = 0; k < this->normals.size(); k++){
	//	rays.push_back(Vector3f(this->normals.at(k).at(0), this->normals.at(k).at(1), this->normals.at(k).at(2)));
	//}

	for (int i = 0; i < 10; i++) {
		rays.push_back(Vector3f(2 * this->functions.uniformRandom() - 1, 2 * this->functions.uniformRandom() - 1, 2 * this->functions.uniformRandom() - 1));
	}

	//rays.push_back(Vector3f(0.5, 0.5, 0.5));
	//rays.push_back(Vector3f(0.0, 1.0, 0.0));

	//rays.push_back(Vector3f(-0.5, -0.5, 0.5));
	//rays.push_back(Vector3f(0.5, -0.5, -0.5));
	//rays.push_back(Vector3f(-0.5, 0.5, -0.5));
	//rays.push_back(Vector3f(-0.5, 0.5, 0.5));
	//rays.push_back(Vector3f(0.5, -0.5, 0.5));
	//rays.push_back(Vector3f(0.5, 0.5, -0.5));

	for (int h = 0; h < rays.size(); h++) {
		if (isInside) {
			r = rays.at(h);
			p0 = Vector3f(point.at(0), point.at(1), point.at(2));
			p1 = p0 + r;
			unsigned int intersectionsCount = 0;

			for (unsigned int i = 0; i < this->faces.size(); i++) {
				vertex0 = this->faces.at(i).at(0) - 1;
				vertex1 = this->faces.at(i).at(3) - 1;
				vertex2 = this->faces.at(i).at(6) - 1;

				normal0 = this->faces.at(i).at(2) - 1;
				normal1 = this->faces.at(i).at(5) - 1;
				normal2 = this->faces.at(i).at(8) - 1;

				n0[0] = this->normals.at(normal0).at(0);
				n0[1] = this->normals.at(normal0).at(1);
				n0[2] = this->normals.at(normal0).at(2);

				n1[0] = this->normals.at(normal1).at(0);
				n1[1] = this->normals.at(normal1).at(1);
				n1[2] = this->normals.at(normal1).at(2);

				n2[0] = this->normals.at(normal2).at(0);
				n2[1] = this->normals.at(normal2).at(1);
				n2[2] = this->normals.at(normal2).at(2);

				n = (n0 + n1 + n2) / 3;

				v0[0] = this->vertices.at(vertex0).at(0);
				v0[1] = this->vertices.at(vertex0).at(1);
				v0[2] = this->vertices.at(vertex0).at(2);

				v1[0] = this->vertices.at(vertex1).at(0);
				v1[1] = this->vertices.at(vertex1).at(1);
				v1[2] = this->vertices.at(vertex1).at(2);

				v2[0] = this->vertices.at(vertex2).at(0);
				v2[1] = this->vertices.at(vertex2).at(1);
				v2[2] = this->vertices.at(vertex2).at(2);

				if (n.dot(p1 - p0) > 0) {
					pri = (n.dot(v0 - p0)) / (n.dot(p1 - p0));
					if (pri > 0) {
						pi = p0 + pri * (p1 - p0);

						w = pi - v0;
						v = v2 - v0;
						u = v1 - v0;

						double ulen = u.norm();
						double vlen = v.norm();

						double stsum = ulen + vlen / 2;
						stsum = 1;

						s = ((u.dot(v))*(w.dot(v)) - (v.dot(v))*(w.dot(u))) / ((u.dot(v))*(u.dot(v)) - (u.dot(u))*(v.dot(v)));
						t = ((u.dot(v))*(w.dot(u)) - (u.dot(u))*(w.dot(v))) / ((u.dot(v))*(u.dot(v)) - (u.dot(u))*(v.dot(v)));

						if ((pi - p1).dot(p1 - p0) > 0) { //Check if ray is ahead or behind triangle
							if ((s > 0) && (t > 0)) {
								if (s + t <= stsum) {
									if ((s < ulen) && (t < vlen)) {
										intersectionsCount++;
									}
								}
							}
						}
					}
				}
			}
			if (intersectionsCount > 0) {
				if (intersectionsCount % 2 == 0)
				{
					isInside = false;
				}
				else
				{
					isInside = true;
				}
			}
			else {
				isInside = false;
			}
		}
	}

	return isInside;
}

bool dotObj::checkIfPointisInsidev2(Vector3f point) {  //Test this
	this->getBoundingBox();
	float maxside = max(max(this->boundingBox.maxX - this->boundingBox.minX, this->boundingBox.maxY - this->boundingBox.minY), this->boundingBox.maxZ - this->boundingBox.minZ);
	Vector3f M = this->findCenterOfMass();
	Vector3f A = point;
	//Vector3f B = M + 40 * maxside*(M - point).normalized();
	Vector3f B = 150 * Vector3f(1, 1, 1);
	//Vector3f Rn = Vector3f(2 * this->functions.uniformRandom() - 1, 2 * this->functions.uniformRandom() - 1, 2 * this->functions.uniformRandom() - 1);
	//B = 40 * maxside*Rn;

	bool isInside = false;
	int vertex0, vertex1, vertex2;
	int normal0, normal1, normal2;
	Vector3f v0, v1, v2, n, w, u, v, r, pointOnPlane;
	double pri, s, t;

	int intersectionsCount = 0;

	for (unsigned int i = 0; i < this->faces.size(); i++) {
		vertex0 = this->faces.at(i).at(0) - 1;
		vertex1 = this->faces.at(i).at(3) - 1;
		vertex2 = this->faces.at(i).at(6) - 1;
		v0 = Vector3f(this->vertices[vertex0][0], this->vertices[vertex0][1], this->vertices[vertex0][2]);
		v1 = Vector3f(this->vertices[vertex1][0], this->vertices[vertex1][1], this->vertices[vertex1][2]);
		v2 = Vector3f(this->vertices[vertex2][0], this->vertices[vertex2][1], this->vertices[vertex2][2]);
		n = (v2 - v0).cross(v1 - v0);
		n.normalize();
		n = -n;

		if (n.dot(B - A) != 0) {
			pri = (n.dot(v0 - A)) / (n.dot(B - A));
			if (pri > 0) {
				pointOnPlane = A + pri * (B - A);

				w = pointOnPlane - v0;
				v = v2 - v0;
				u = v1 - v0;
				s = ((u.dot(v))*(w.dot(v)) - (v.dot(v))*(w.dot(u))) / ((u.dot(v))*(u.dot(v)) - (u.dot(u))*(v.dot(v)));
				t = ((u.dot(v))*(w.dot(u)) - (u.dot(u))*(w.dot(v))) / ((u.dot(v))*(u.dot(v)) - (u.dot(u))*(v.dot(v)));

				if ((s > 0) && (t > 0)) {
					if (s + t < 1) {
						if ((s < u.norm()) && (t < v.norm())) {
							float BA = (B - A).norm();
							float pA = (pointOnPlane - A).norm();
							float Bp = (B - pointOnPlane).norm();
							float diff = BA - pA - Bp;
							if (diff < 0.001) {
								intersectionsCount++;
							}
						}
					}
				}
			}
		}
	}
	if (intersectionsCount > 0) {
		if (intersectionsCount % 2 == 1)
		{
			isInside = true;
		}
		else
		{
			isInside = false;
		}
	}
	else {
		isInside = false;
	}

	//if (intersectionsCount > 2){ std::cout << intersectionsCount << std::endl; }
	return isInside;
	//return true;
}

void volume::exportToFile(std::string output) {
	std::ofstream outfile;
	outfile.open(output);
	for (unsigned int i = 0; i < this->vertices.size(); i++)
	{
		outfile << "v" << " " << this->vertices.at(i).at(0) << " " << this->vertices.at(i).at(1) << " " << this->vertices.at(i).at(2) << std::endl;
	}

	outfile.close();

	std::cout << "File Export Complete" << std::endl;

	return;
}

Vector3f volume::setCenterOfMass(void) {
	Vector3f center = { 0.0, 0.0, 0.0 };

	int numOfPoints = this->vertices.size();

	if (numOfPoints > 0) {
		for (int i = 0; i < numOfPoints; i++) {
			center = center + Vector3f({ this->vertices.at(i).at(0), this->vertices.at(i).at(1), this->vertices.at(i).at(2) });
		}
		center = center / numOfPoints;
	}

	this->centerOfMass = center;

	return center;
}

void volume::splitVolume(plane splittingPlane, volume &subVolumeA, volume &subVolumeB) {
	Vector3f somePoint;
	volume subsetA;
	volume subsetB;
#pragma omp parallel private(somePoint)
#pragma omp for
	for (int i = 0; i < this->vertices.size(); i++) {
		somePoint[0] = this->vertices.at(i).at(0);
		somePoint[1] = this->vertices.at(i).at(1);
		somePoint[2] = this->vertices.at(i).at(2);
		if (splittingPlane.checkWhichSideOfPlaneIsPoint(somePoint))
		{
#pragma omp critical
			subVolumeA.vertices.push_back(this->vertices.at(i));
		}
		else
		{
#pragma omp critical
			subVolumeB.vertices.push_back(this->vertices.at(i));
		}
	}

	return;
}

void volume::defineSplittingPlane(std::string method) {
	if (method == "PCA") {
		PCA * thePCA;
		thePCA = new PCA();
		thePCA->setOfPoints = this->vertices;
		//std::cout << "Number of vertices per PCA operation: "<<this->vertices.size()<< std::endl;
		thePCA->getPCA();
		this->splittingPlane.pointOnPlane = this->centerOfMass;
		this->splittingPlane.normal = thePCA->eigenVectorsSorted.at(thePCA->eigenVectorsSorted.size() - 1);
	}



	if (method == "Random") {
		this->splittingPlane.pointOnPlane = this->centerOfMass;
		this->splittingPlane.normal.setRandom();
	}



	return;
}

void volume::defineSplittingPlane(Vector3f & p1, Vector3f & p2, Vector3f & p3) {
	//Vector3f center = (p3 + p2) / 2;
	//Vector3f up = (p3 - p2).normalized();

	this->splittingPlane.pointOnPlane = (p3 + p2) / 2;
	this->splittingPlane.normal = (p3 - p2).normalized();

	return;
}