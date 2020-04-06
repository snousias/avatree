#include "lungModelling.h"

void ggraph::ggraph2ModelColorize(dotObj * result, bool jump2Bifurcations, std::string colorize) {
	result->segment_property_map_per_line.clear();
	result->vertices.push_back({
		this->nodes[this->init->index].position.x(),
		this->nodes[this->init->index].position.y(),
		this->nodes[this->init->index].position.z()
		});
	std::vector<gnode*> A;
	//A.push_back(&this->nodes[this->init->index]);
	//this->nodes[this->init->index].index = result->vertices.size() - 1;
	A.push_back(this->init);
	this->init->index = result->vertices.size() - 1;
	this->nodes[this->init->index].index = result->vertices.size() - 1;
	if (jump2Bifurcations) { this->init->nextBifurcationPtr[0]->previousBifurcationPtr = this->init; }
	while (!A.empty()) {
		std::vector<gnode*> B;
		if (!jump2Bifurcations) { if (!A.back()->doStop)B = A.back()->nextNodesPtr; }
		if (jump2Bifurcations) { if (!A.back()->doStop)B = A.back()->nextBifurcationPtr; }
		for (int i = 0; i < B.size(); i++) {
			//gnode * n = B[i];
			//Process

			result->vertices.push_back({ B[i]->position.x(), B[i]->position.y(), B[i]->position.z() });

			B[i]->index = result->vertices.size() - 1;
			if (!jump2Bifurcations) {
				result->lines.push_back({ A.back()->index + 1, B[i]->index + 1 });
				result->segment_property_map_per_line.push_back(B[i]->HorsfieldOrder);
				//generation
				//HorsfieldOrder
				//StrahlerOrder
			}
			else {
				result->lines.push_back({ B[i]->previousBifurcationPtr->index + 1, B[i]->index + 1 });
				result->segment_property_map_per_line.push_back(B[i]->HorsfieldOrder);
				//generation
				//HorsfieldOrder
				//StrahlerOrder
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}
	return;
}

void ggraph::ggraph2Model(dotObj * result, bool jump2Bifurcations, std::string colorize) {
	result->vertices.push_back({
		this->nodes[this->init->index].position.x(),
		this->nodes[this->init->index].position.y(),
		this->nodes[this->init->index].position.z()
		});
	std::vector<gnode*> A;
	//A.push_back(&this->nodes[this->init->index]);
	//this->nodes[this->init->index].index = result->vertices.size() - 1;
	A.push_back(this->init);
	this->init->index = result->vertices.size() - 1;
	this->nodes[this->init->index].index = result->vertices.size() - 1;
	if (jump2Bifurcations) { this->init->nextBifurcationPtr[0]->previousBifurcationPtr = this->init; }
	while (!A.empty()) {
		std::vector<gnode*> B;
		if (!jump2Bifurcations) { if (!A.back()->doStop)B = A.back()->nextNodesPtr; }
		if (jump2Bifurcations) { if (!A.back()->doStop)B = A.back()->nextBifurcationPtr; }
		for (int i = 0; i < B.size(); i++) {
			//gnode * n = B[i];
			//Process

			result->vertices.push_back({ B[i]->position.x(), B[i]->position.y(), B[i]->position.z() });

			B[i]->index = result->vertices.size() - 1;
			if (!jump2Bifurcations) {
				//result->lines.push_back({ B[i]->previousNodePtr->index + 1, B[i]->index + 1 });
				result->lines.push_back({ A.back()->index + 1, B[i]->index + 1 });
			}
			else {
				result->lines.push_back({ B[i]->previousBifurcationPtr->index + 1, B[i]->index + 1 });
				//result->lines.push_back({ A.back()->index + 1, B[i]->index + 1 });
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}
	return;
}

void ggraph::ggraph2ModelV(dotObj * result) {
	std::cout << "ggraph::ggraph2ModelV" << std::endl;
	result->vertices.push_back({
		this->nodes[this->init->index].position.x(),
		this->nodes[this->init->index].position.y(),
		this->nodes[this->init->index].position.z()
		});
	std::vector<gnode*> A;
	A.push_back(this->init);
	this->init->index = result->vertices.size() - 1;
	while (!A.empty()) {
		std::vector<gnode*> B;
		B = A.back()->nextNodesPtr;
		for (int i = 0; i < B.size(); i++) {
			if (B[i]->doStop) {
				B.erase(B.begin() + i);
			}
		}
		A.pop_back();
		for (int i = 0; i < B.size(); i++) {
			gnode * n = B[i];
			//Process
			for (int j = 0; j < B[i]->meshVerticesPositions.size(); j++) {
				result->vertices.push_back({ B[i]->meshVerticesPositions[j].x(), B[i]->meshVerticesPositions[j].y(), B[i]->meshVerticesPositions[j].z() });
			}
		}
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}
	return;
}

void dotObj::splitSkeletonEdges(float minimumDistance) {
	std::cout << "dotObj::splitSkeletonEdges" << std::endl;
	for (int i = 0; i < this->lines.size(); i++) {
		int vindex0 = this->lines[i][0] - 1;
		int vindex1 = this->lines[i][1] - 1;

		Vector3f v0 = Vector3f(this->vertices[vindex0][0], this->vertices[vindex0][1], this->vertices[vindex0][2]);
		Vector3f v1 = Vector3f(this->vertices[vindex1][0], this->vertices[vindex1][1], this->vertices[vindex1][2]);

		if ((v1 - v0).norm() > minimumDistance) {
			Vector3f newPos = v0 + minimumDistance * (v1 - v0);
			this->vertices.push_back({ newPos[0], newPos[1], newPos[2] });
			int ni = this->vertices.size() - 1;
			this->lines[i][1] = ni + 1;
			this->lines.push_back({ ni + 1, vindex1 + 1 });
			i = 0;
		}
	}

	std::cout << "Split skeleton complete" << std::endl;

	return;
}

void ggraph::exportGraphFeatures(std::string outfile) {
	std::cout << " ggraph::exportGraphFeatures" << std::endl;

	std::vector<gnode*> A;
	A.push_back(this->init);
	std::ofstream mFile;
	mFile.open(outfile);
	mFile << "[{}";
	while (!A.empty()) {
		gnode* q = A.back();
		//std::vector<gnode*> B = q->nextBifurcationPtr;
		std::vector<gnode*> M = q->nextBifurcationPtr;

		Vector3f b1, b2, p;

		std::vector<float> theta = { 0.0 , 0.0 };

		float thetatot = 0.0;

		/*if (B.size() == 2) {
			p = -(q->position - q->previousBifurcationPtr->position);
			b1 = B[0]->position - q->position;
			b2 = B[1]->position - q->position;
			thetatot = acos(b1.dot(b2) / (b1.norm()*b2.norm())) * 360 / (2 * Pi);
			theta[0] = 180 - (acos(b1.dot(p) / (b1.norm()*p.norm())) * 360 / (2 * Pi));
			theta[1] = 180 - (acos(b2.dot(p) / (b2.norm()*p.norm())) * 360 / (2 * Pi));
		}*/

		for (int i = 0; i < M.size(); i++) {
			if (!M.empty()) {
				//Process
				//len = (M[i]->position - M[i]->previousNodePtr->position).norm();
				std::string LR = M[i]->isLeft ? "R" : "L";

				mFile << ",{"
					<< "\"horsfieldOrder\":" << M[i]->HorsfieldOrder << ","
					<< "\"strahlerOrder\":" << M[i]->StrahlerOrderStage_1 << ","
					<< "\"generation\":" << M[i]->generation - 1 << ","
					<< "\"diameter\":" << M[i]->diameterBranch << ","
					<< "\"lengthOfNode\":" << M[i]->length << ","
					<< "\"lengthOfBranch\":" << M[i]->lengthOfBranch << ","
					<< "\"thetaToParent\":" << M[i]->thetaToParent << ","
					<< "\"thetaToSibling\":" << M[i]->thetaToSimbling << ","
					<< "\"RB\":" << M[i]->RB << ","
					<< "\"RD\":" << M[i]->RD << ","
					<< "\"RL\":" << M[i]->RL << ","
					<< "\"RBH\":" << M[i]->RBH << ","
					<< "\"RDH\":" << M[i]->RDH << ","
					<< "\"RLH\":" << M[i]->RLH << ","
					<< "\"RBS\":" << M[i]->RBS << ","
					<< "\"RDS\":" << M[i]->RDS << ","
					<< "\"RLS\":" << M[i]->RLS << ","
					<< "\"leftOrRight\":" << "\"" << LR << "\"" << "}" << std::endl;
			}
		}
		A.pop_back();
		A.insert(A.end(), M.begin(), M.end());
	}
	mFile << "]";
	mFile.close();
	return;
}

void ggraph::exportGraphFeatures(std::string outfile, gnode* inode) {
	std::cout << " ggraph::exportGraphFeatures" << std::endl;

	std::vector<gnode*> A;
	A.push_back(inode);
	std::ofstream mFile;
	mFile.open(outfile);
	mFile << "[{}";
	while (!A.empty()) {
		gnode* q = A.back();
		//std::vector<gnode*> B = q->nextBifurcationPtr;
		std::vector<gnode*> M = q->nextBifurcationPtr;

		Vector3f b1, b2, p;

		std::vector<float> theta = { 0.0 , 0.0 };

		float thetatot = 0.0;

		/*if (B.size() == 2) {
			p = -(q->position - q->previousBifurcationPtr->position);
			b1 = B[0]->position - q->position;
			b2 = B[1]->position - q->position;
			thetatot = acos(b1.dot(b2) / (b1.norm()*b2.norm())) * 360 / (2 * Pi);
			theta[0] = 180 - (acos(b1.dot(p) / (b1.norm()*p.norm())) * 360 / (2 * Pi));
			theta[1] = 180 - (acos(b2.dot(p) / (b2.norm()*p.norm())) * 360 / (2 * Pi));
		}*/

		for (int i = 0; i < M.size(); i++) {
			if (!M.empty()) {
				//Process
				//len = (M[i]->position - M[i]->previousNodePtr->position).norm();
				std::string LR = M[i]->isLeft ? "R" : "L";

				mFile << ",{"
					<< "\"horsfieldOrder\":" << M[i]->HorsfieldOrder << ","
					<< "\"strahlerOrder\":" << M[i]->StrahlerOrderStage_1 << ","
					<< "\"generation\":" << M[i]->generation - 1 << ","
					<< "\"diameter\":" << M[i]->diameterBranch << ","
					<< "\"lengthOfNode\":" << M[i]->length << ","
					<< "\"lengthOfBranch\":" << M[i]->lengthOfBranch << ","
					<< "\"thetaToParent\":" << M[i]->thetaToParent << ","
					<< "\"thetaToSibling\":" << M[i]->thetaToSimbling << ","
					<< "\"RB\":" << M[i]->RB << ","
					<< "\"RD\":" << M[i]->RD << ","
					<< "\"RL\":" << M[i]->RL << ","
					<< "\"leftOrRight\":" << "\"" << LR << "\"" << "}" << std::endl;
			}
		}
		A.pop_back();
		A.insert(A.end(), M.begin(), M.end());
	}
	mFile << "]";
	mFile.close();
	return;
}

void ggraph::initializeDiameters(void) {
	std::cout << " ggraph::exportGraphFeatures" << std::endl;

	/*std::vector<gnode*> A;
	A.push_back(this->init);
	while (!A.empty()) {
		gnode* q = A.back();
		std::vector<gnode*> M = q->nextNodesPtr;
		A.pop_back();
		for (int i = 0; i < M.size(); i++) {
			if ((isnan(M[i]->diameter) || M[i]->diameter < 0) && q->diameter >= 0) {
				M[i]->diameter = q->diameter;
			}
		}
		A.insert(A.end(), M.begin(), M.end());
	}
	A.clear();*/

	std::vector<gnode*> A;
	A.clear();
	A.push_back(this->init);
	while (!A.empty()) {
		gnode* q = A.back();
		std::vector<gnode*> M = q->nextNodesPtr;
		//std::vector<gnode*> B = q->nextNodesPtr;

		for (int i = 0; i < M.size(); i++) {
			if (!M[i]->isBifurcation && !isnan(M[i]->diameter)) {
				if (M[i]->diameter > 0) {
					if (M[i]->nextBifurcationPtr.size() > 0) {
						M[i]->nextBifurcationPtr[0]->diametersToProcess.push_back(M[i]->diameter);
					}
				}
			}
		}
		A.pop_back();
		A.insert(A.end(), M.begin(), M.end());
	}

	A.clear();
	A.push_back(this->init);
	while (!A.empty()) {
		gnode* q = A.back();
		std::vector<gnode*> B = q->nextBifurcationPtr;
		A.pop_back();
		for (int i = 0; i < B.size(); i++) {
			if (B[i]->diametersToProcess.size() > 0) {
				B[i]->diameterBranch = std::accumulate(B[i]->diametersToProcess.begin(), B[i]->diametersToProcess.end(), 0.0) / B[i]->diametersToProcess.size();
				std::cout << B[i]->diameterBranch << std::endl;
				B[i]->diametersToProcess.clear();
			}
		}
		A.insert(A.end(), B.begin(), B.end());
	}

	A.clear();
	A.push_back(this->init);
	while (!A.empty()) {
		gnode* q = A.back();
		std::vector<gnode*> B = q->nextBifurcationPtr;

		for (int i = 0; i < B.size(); i++) {
			if (B.size() == 2) {
				Vector3f p = (q->position - q->previousBifurcationPtr->position);
				Vector3f b1 = B[0]->position - q->position;
				Vector3f b2 = B[1]->position - q->position;
				Vector3f b = B[i]->position - q->position;
				b.normalize();
				p.normalize();
				b1.normalize();
				b2.normalize();

				B[i]->thetaToParent = (acos(b.dot(p) / (b.norm()*p.norm())) * 360 / (2 * Pi));
				B[i]->thetaToSimbling = acos(b1.dot(b2) / (b1.norm()*b2.norm())) * 360 / (2 * Pi);;
			}
		}
		A.pop_back();
		A.insert(A.end(), B.begin(), B.end());
	}
	A.clear();

	/*A.push_back(this->init);
	while (!A.empty()) {
		gnode* q = A.back();
		std::vector<gnode*> B = q->nextBifurcationPtr;
		A.pop_back();
		Vector3f b1, b2, p;
		std::vector<float> theta = { 0.0 , 0.0 };
		float thetatot = 0.0;
		if (B.size() == 2) {
			p = -(q->position - q->previousBifurcationPtr->position);
			b1 = B[0]->position - q->position;
			b2 = B[1]->position - q->position;
			thetatot = acos(b1.dot(b2) / (b1.norm()*b2.norm())) * 360 / (2 * Pi);
			theta[0] = 180 - (acos(b1.dot(p) / (b1.norm()*p.norm())) * 360 / (2 * Pi));
			theta[1] = 180 - (acos(b2.dot(p) / (b2.norm()*p.norm())) * 360 / (2 * Pi));

			q->nextBifurcationPtr[0]->thetaToSimbling = thetatot;
			q->nextBifurcationPtr[0]->thetaToParent = theta[0];

			q->nextBifurcationPtr[1]->thetaToSimbling = thetatot;
			q->nextBifurcationPtr[1]->thetaToParent = theta[1];
		}
	}*/

	return;
}

////==================================================================///

void gnode::copyNodeProperties(gnode * source, gnode * target) {
	target->isBifurcation = source->isBifurcation;
	target->isTerminal = source->isTerminal;
	target->isInlet = source->isInlet;
	target->isMedian = source->isMedian;
	target->isLeft = source->isLeft;
	target->isRight = source->isRight;
	target->diameter = source->diameter;
	target->generation = source->generation;
	target->segmentIndex = source->segmentIndex;
	target->HorsfieldOrder = source->HorsfieldOrder;
	target->StrahlerOrderStage_1 = source->StrahlerOrderStage_1;
	target->StrahlerOrderStage_2 = source->StrahlerOrderStage_2;
	target->doStop = false;
	target->position = source->position;
	//target->index = source->index;
	//target->line;
	return;
}

ggraph ggraph::generateGraphFromGraph(void) {
	ggraph output;

	gnode * n = new gnode;
	n->copyNodeProperties(this->init, n);
	output.nodes.push_back(*n);
	output.init = &output.nodes[0];
	delete n;

	std::cout << "Init start" << std::endl;

	std::vector<gnode*> A;
	A.push_back(this->init);
	A.back()->correspondingNodePtr = &output.nodes[0];
	while (!A.empty()) {
		std::vector<gnode*> B;
		if (!A.back()->doStop)B = A.back()->nextNodesPtr;
		for (int i = 0; i < B.size(); i++) {
			gnode * n = new gnode;
			n->copyNodeProperties(B[i], n);
			output.nodes.push_back(*n);
			delete n;
			B[i]->correspondingNodePtr = &output.nodes[output.nodes.size() - 1];
			A.back()->correspondingNodePtr->nextNodesPtr.push_back(&output.nodes[output.nodes.size() - 1]);
			//B[i]->correspondingNodePtr->previousNodePtr = A.back()->correspondingNodePtr;
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}

	/*

	result->vertices.push_back({
		this->nodes[this->init->index].position.x(),
		this->nodes[this->init->index].position.y(),
		this->nodes[this->init->index].position.z()
		});

	//A.push_back(&this->nodes[this->init->index]);
	//this->nodes[this->init->index].index = result->vertices.size() - 1;

	this->init->index = result->vertices.size() - 1;
	this->nodes[this->init->index].index = result->vertices.size() - 1;
	if (jump2Bifurcations) { this->init->nextBifurcationPtr[0]->previousBifurcationPtr = this->init; }
	while (!A.empty()) {
		std::vector<gnode*> B;
		if (!jump2Bifurcations) { if (!A.back()->doStop)B = A.back()->nextNodesPtr; }
		if (jump2Bifurcations) { if (!A.back()->doStop)B = A.back()->nextBifurcationPtr; }
		for (int i = 0; i < B.size(); i++) {
			//gnode * n = B[i];
			//Process

			result->vertices.push_back({ B[i]->position.x(), B[i]->position.y(), B[i]->position.z() });

			B[i]->index = result->vertices.size() - 1;
			if (!jump2Bifurcations) {
				//result->lines.push_back({ B[i]->previousNodePtr->index + 1, B[i]->index + 1 });
				result->lines.push_back({ A.back()->index + 1, B[i]->index + 1 });
			}
			else {
				result->lines.push_back({ B[i]->previousBifurcationPtr->index + 1, B[i]->index + 1 });
				//result->lines.push_back({ A.back()->index + 1, B[i]->index + 1 });
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}

	*/

	return output;
}

void ggraph::generateGraphFeaturesGeneration(void) {
	std::vector<gnode*> A;

	//this->graph.nodes[this->inlet].nextBifurcationPtr[0]->generation = 0;

	this->init->generation = 0;
	this->init->segmentIndex = 0;
	this->init->HorsfieldOrder = -1;
	this->init->StrahlerOrderStage_1 = -1;
	A.clear();
	A.push_back(this->init);
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
	A.push_back(this->init);
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

	while (this->init->HorsfieldOrder < 1) {
		A.clear();
		A.push_back(this->init);
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
	A.push_back(this->init);
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

void ggraph::generateGraphFeaturesLRDiscrimination(void) {
	gnode * g1 = new gnode();

	g1 = this->init->nextBifurcationPtr[0];

	std::vector < gnode *> g2 = g1->nextBifurcationPtr;

	std::cout << "Crack here 1" << std::endl;
	if ((g2[0]->position - g1->position).norm() > (g2[1]->position - g1->position).norm()) {
		g1->nextNodesPtr[0]->isRight = true;
		g1->nextNodesPtr[0]->isLeft = false;

		g1->nextNodesPtr[1]->isLeft = true;
		g1->nextNodesPtr[1]->isRight = false;
	}
	else {
		g1->nextNodesPtr[0]->isLeft = true;
		g1->nextNodesPtr[0]->isRight = false;

		g1->nextNodesPtr[1]->isRight = true;
		g1->nextNodesPtr[1]->isLeft = false;
	}

	std::cout << "Crack here 2" << std::endl;
	std::vector<gnode*> A;
	A.clear();
	A.push_back(this->init);
	std::cout << "Crack here 3" << std::endl;
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextNodesPtr;

		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				//Process
				if (A.back()->isLeft != A.back()->isRight) {
					B[i]->isLeft = A.back()->isLeft;
					B[i]->isRight = A.back()->isRight;
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

void ggraph::generateGraphFeaturesDiameter(void) {
	std::vector<gnode*> A;
	A.clear();
	A.push_back(this->init);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextBifurcationPtr;
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				//Process

				if (!B[i]->diameterBranch > 0) {
					//B[i]->diameterBranch = pow(2.0, -1.0 / 3.0)*A.back()->diameterBranch;
					B[i]->diameterBranch = sqrt(sin(B[i]->thetaToParent*Pi / 180) / sin(B[i]->thetaToSimbling*Pi / 180))*A.back()->diameterBranch;
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

void ggraph::generateGraphFeaturesRatios(void) {
	std::vector<gnode*> A;
	A.clear();
	A.push_back(this->init);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextBifurcationPtr;
		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				//Process
				if ((A.back()->diameterBranch > 0) && (B[i]->diameterBranch > 0)) {
					A.back()->RD = B[i]->diameterBranch / A.back()->diameterBranch;
				}
				if ((A.back()->lengthOfBranch > 0) && (B[i]->lengthOfBranch > 0)) {
					A.back()->RL = B[i]->lengthOfBranch / A.back()->lengthOfBranch;
				}
				if ((A.back()->thetaToParent > 0) && (B[i]->thetaToParent > 0)) {
					A.back()->RB = B[i]->thetaToParent / A.back()->thetaToParent;
				}
				if (A.back()->StrahlerOrderStage_1 != B[i]->StrahlerOrderStage_1) {
					if ((A.back()->diameterBranch > 0) && (B[i]->diameterBranch > 0)) {
						A.back()->RDS = 1 / (B[i]->diameterBranch / A.back()->diameterBranch);
					}
					if ((A.back()->lengthOfBranch > 0) && (B[i]->lengthOfBranch > 0)) {
						A.back()->RLS = 1 / (B[i]->lengthOfBranch / A.back()->lengthOfBranch);
					}
					if ((A.back()->thetaToParent > 0) && (B[i]->thetaToParent > 0)) {
						A.back()->RBS = 1 / (B[i]->thetaToParent / A.back()->thetaToParent);
					}
				}
				if (A.back()->HorsfieldOrder != B[i]->HorsfieldOrder) {
					if ((A.back()->diameterBranch > 0) && (B[i]->diameterBranch > 0)) {
						A.back()->RDH = 1 / (B[i]->diameterBranch / A.back()->diameterBranch);
					}
					if ((A.back()->lengthOfBranch > 0) && (B[i]->lengthOfBranch > 0)) {
						A.back()->RLH = 1 / (B[i]->lengthOfBranch / A.back()->lengthOfBranch);
					}
					if ((A.back()->thetaToParent > 0) && (B[i]->thetaToParent > 0)) {
						A.back()->RBH = 1 / (B[i]->thetaToParent / A.back()->thetaToParent);
					}
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

void ggraph::generateGraphFeaturesLength(std::string type) {
	std::vector<gnode*> A;

	if (type == "branch") {
		A.clear();
		A.push_back(this->init);
		while (!A.empty()) {
			std::vector<gnode*> B = A.back()->nextBifurcationPtr;
			for (int i = 0; i < B.size(); i++) {
				if (!B.empty()) {
					B[i]->lengthOfBranch = (B[i]->position - A.back()->position).norm();
				}
			}
			A.pop_back();
			if (!B.empty()) {
				A.insert(A.end(), B.begin(), B.end());
			}
		}
	}
	else {
		if (type == "node") {
			A.clear();
			A.push_back(this->init);
			while (!A.empty()) {
				std::vector<gnode*> B = A.back()->nextNodesPtr;
				for (int i = 0; i < B.size(); i++) {
					if (!B.empty()) {
						B[i]->length = (B[i]->position - A.back()->position).norm();
					}
				}
				A.pop_back();
				if (!B.empty()) {
					A.insert(A.end(), B.begin(), B.end());
				}
			}
		}
	}
	return;
}

void ggraph::propagateGraphFeatures(void) {
	this->init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->diameter =
		this->init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->diameter;

	this->init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->diameter =
		this->init->nextBifurcationPtr[0]->nextNodesPtr[0]->nextMedianPtr->diameter;

	this->init->nextBifurcationPtr[0]->nextBifurcationPtr[1]->diameterBranch =
		this->init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->diameter;

	this->init->nextBifurcationPtr[0]->nextBifurcationPtr[1]->diameterBranch =
		this->init->nextBifurcationPtr[0]->nextNodesPtr[1]->nextMedianPtr->diameter;

	std::vector<gnode*> A;
	A.clear();
	A.push_back(this->init);
	while (!A.empty()) {
		gnode* q = A.back();
		std::vector<gnode*> M = q->nextNodesPtr;
		A.pop_back();
		for (int i = 0; i < M.size(); i++) {
			if (M[i]->isTerminal) {
				std::cout << "Terminal:" << M[i]->previousBifurcationPtr->nextBifurcationPtr[i]->diameterBranch << std::endl;
				//M[i]->diameterBranch = M[i]->previousBifurcationPtr->diameterBranch;
			}
		}
		A.insert(A.end(), M.begin(), M.end());
	}

	A.clear();
	A.push_back(this->init);
	while (!A.empty()) {
		gnode* q = A.back();
		std::vector<gnode*> M = q->nextNodesPtr;
		for (int i = 0; i < M.size(); i++) {
			if (M[i]->nextBifurcationPtr.size() > 0 && !M[i]->isBifurcation) {
				if (M[i]->nextBifurcationPtr[0]->diameterBranch > 0) {
					M[i]->diameter = M[i]->nextBifurcationPtr[0]->diameterBranch;
				}
			}
			if (M[i]->isBifurcation) {
				M[i]->diameter = M[i]->diameterBranch;
			}
			if (M[i]->isTerminal) {
				M[i]->diameter = M[i]->diameterBranch;
			}
		}
		A.pop_back();
		A.insert(A.end(), M.begin(), M.end());
	}

	std::cout << "Propagate LR" << std::endl;
	A.clear();
	A.push_back(this->init);
	while (!A.empty()) {
		std::vector<gnode*> B = A.back()->nextNodesPtr;

		for (int i = 0; i < B.size(); i++) {
			if (!B.empty()) {
				//Process
				if (A.back()->isLeft != A.back()->isRight) {
					B[i]->isLeft = A.back()->isLeft;
					B[i]->isRight = A.back()->isRight;
				}
			}
		}
		A.pop_back();
		if (!B.empty()) {
			A.insert(A.end(), B.begin(), B.end());
		}
	}

	std::cout << "Propagate Ratios" << std::endl;
	A.clear();
	A.push_back(this->init);
	while (!A.empty()) {
		gnode* q = A.back();
		std::vector<gnode*> M = q->nextNodesPtr;

		for (int i = 0; i < M.size(); i++) {
			if (M[i]->nextBifurcationPtr.size() > 0 && !M[i]->isBifurcation) {
				if (M[i]->nextBifurcationPtr[0]->RB > 0) {
					M[i]->RB = M[i]->nextBifurcationPtr[0]->RB;
				}
				if (M[i]->nextBifurcationPtr[0]->RD > 0) {
					M[i]->RD = M[i]->nextBifurcationPtr[0]->RD;
				}
				if (M[i]->nextBifurcationPtr[0]->RL > 0) {
					M[i]->RL = M[i]->nextBifurcationPtr[0]->RL;
				}
			}
		}
		A.pop_back();
		A.insert(A.end(), M.begin(), M.end());
	}

	return;
}

void ggraph::enrichGraph(float minimumDistance) { //<==Pretty fucked up.
	std::cout << "ggraph::enrichGraph" << std::endl;
	std::vector<gnode*> A;

	int iterationNumber = 0;

	bool  processComplete = false;
	while (!processComplete) {
		std::cout << "Graph enrichment iteration : " << iterationNumber++ << std::endl;
		processComplete = true;
		A.clear();
		A.push_back(this->init);
		while (!A.empty()) {
			std::vector<gnode*> B;
			B = A.back()->nextNodesPtr;
			for (int i = 0; i < B.size(); i++) {
				//Process
				if ((B[i]->position - A.back()->position).norm() > 2 * minimumDistance) {
					processComplete = false;
					//gnode n;
					//this->nodes.push_back(n);
					//gnode * m = &this->nodes[this->nodes.size() - 1];

					gnode * m = new gnode();
					m->position = A.back()->position + minimumDistance * (B[i]->position - A.back()->position) / (B[i]->position - A.back()->position).norm();
					m->previousNodePtr = A.back();
					if (A.back()->isBifurcation) {
						m->previousBifurcationPtr = A.back();
					}
					else {
						m->previousBifurcationPtr = A.back()->previousBifurcationPtr;
					}
					m->nextNodesPtr.push_back(B[i]);
					if (B[i]->isBifurcation) {
						m->nextBifurcationPtr.push_back(B[i]);
					}
					else {
						m->nextBifurcationPtr = B[i]->nextBifurcationPtr;
					}
					A.back()->nextNodesPtr[i] = m;
					B[i]->previousNodePtr = m;
				}
			}
			A.pop_back();
			if (!B.empty()) {
				A.insert(A.end(), B.begin(), B.end());
			}
		}
	}
	return;
}

void ggraph::getLobes(void) {
	helper * functions = new helper();
	std::vector<gnode*> A;
	Vector3f x;
	if (this->init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->isRight) {
		A = this->init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->nextBifurcationPtr;
		x = functions->projPointOnVector(A[0]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
		float d1 = (x - this->init->position).norm();
		x = functions->projPointOnVector(A[1]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
		float d2 = (x - this->init->position).norm();
		if (d1 <= d2) {
			this->LU = A[0];
			this->LL = A[1];
		}
		else {
			this->LU = A[0];
			this->LL = A[1];
		}
	}

	if (this->init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->isLeft) {
		A = this->init->nextBifurcationPtr[0]->nextBifurcationPtr[1]->nextBifurcationPtr;
		x = functions->projPointOnVector(A[0]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
		float d1 = (x - this->init->position).norm();
		x = functions->projPointOnVector(A[1]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
		float d2 = (x - this->init->position).norm();
		if (d1 <= d2) {
			this->LU = A[0];
			this->LL = A[1];
		}
		else {
			this->LU = A[1];
			this->LL = A[0];
		}
	}
	A.clear();
	if (this->init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->isLeft) {
		A = this->init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->nextBifurcationPtr;
		x = functions->projPointOnVector(A[0]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
		float d1 = (x - this->init->position).norm();
		x = functions->projPointOnVector(A[1]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
		float d2 = (x - this->init->position).norm();
		if (d1 <= d2) {
			this->RU = A[0];
			A.clear();
			A = this->init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->nextBifurcationPtr[1]->nextBifurcationPtr;
			x = functions->projPointOnVector(A[0]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
			float f1 = (x - this->init->position).norm();
			x = functions->projPointOnVector(A[1]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
			float f2 = (x - this->init->position).norm();
			if (f1 <= f2) {
				this->RM = A[0];
				this->RL = A[1];
			}
			else {
				this->RM = A[1];
				this->RL = A[0];
			}
		}
		else {
			this->RU = A[1];
			A.clear();
			A = this->init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->nextBifurcationPtr[0]->nextBifurcationPtr;
			x = functions->projPointOnVector(A[0]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
			float f1 = (x - this->init->position).norm();
			x = functions->projPointOnVector(A[1]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
			float f2 = (x - this->init->position).norm();
			if (f1 <= f2) {
				this->RM = A[0];
				this->RL = A[1];
			}
			else {
				this->RM = A[1];
				this->RL = A[0];
			}
		}
	}
	if (this->init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->isRight) {
		A = this->init->nextBifurcationPtr[0]->nextBifurcationPtr[1]->nextBifurcationPtr;
		x = functions->projPointOnVector(A[0]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
		float d1 = (x - this->init->position).norm();
		x = functions->projPointOnVector(A[1]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
		float d2 = (x - this->init->position).norm();
		if (d1 <= d2) {
			this->RU = A[0];
			A.clear();
			A = this->init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->nextBifurcationPtr[1]->nextBifurcationPtr;
			x = functions->projPointOnVector(A[0]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
			float f1 = (x - this->init->position).norm();
			x = functions->projPointOnVector(A[1]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
			float f2 = (x - this->init->position).norm();
			if (f1 <= f2) {
				this->RM = A[0];
				this->RL = A[1];
			}
			else {
				this->RM = A[1];
				this->RL = A[0];
			}
		}
		else {
			this->RU = A[1];
			A.clear();
			A = this->init->nextBifurcationPtr[0]->nextBifurcationPtr[0]->nextBifurcationPtr[0]->nextBifurcationPtr;
			x = functions->projPointOnVector(A[0]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
			float f1 = (x - this->init->position).norm();
			x = functions->projPointOnVector(A[1]->position, this->init->position, this->init->nextBifurcationPtr[0]->position);
			float f2 = (x - this->init->position).norm();
			if (f1 <= f2) {
				this->RM = A[0];
				this->RL = A[1];
			}
			else {
				this->RM = A[1];
				this->RL = A[0];
			}
		}
	}
	return;
}