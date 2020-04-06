#include "lungModelling.h"


//gnode constructor
gnode::gnode(void){
	this->isBifurcation = false;
	this->isTerminal = false;
	this->isInlet = false;
	this->isMedian = false;
	this->isLeft = false;
	this->isRight = false;
	this->doStop = false;
	this->index = -1;
	int line = -1;
	this->generation = -1;

	this->previousBifurcation = -1;
	this->previous = -1;
	this->HorsfieldOrder = -1;
	this->diameter = -1;
	this->diameterBranch = 0;
	this->StrahlerOrderStage_1 = -1;
	this->StrahlerOrderStage_2 = -1;
	this->thetaToParent = -20.0;
	this->thetaToSimbling = -10.0;
	this->length = 0.0;
	this->lengthOfBranch = 0.0;
	this->RB = -1;
	this->RD = -1;
	this->RL = -1;
	this->RBS = -1;
	this->RDS = -1;
	this->RLS = -1;
	this->RBH = -1;
	this->RDH = -1;
	this->RLH = -1;
}
