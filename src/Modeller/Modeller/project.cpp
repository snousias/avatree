#include "lungModelling.h"


void project::update(void){
	*this->tempModel = this->lungmodel.back();

	return;
}