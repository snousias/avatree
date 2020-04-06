#include "lungModelling.h"


void dotObj::tic(){
	ticTocTimerBegin = new clock_t();
	*ticTocTimerBegin = clock();

	return;
}

void dotObj::toc(){
	ticTocTimerCheckPoint = new clock_t();

	*ticTocTimerCheckPoint = clock();
	this->ticTocTimerElapsedSecs = double(*ticTocTimerCheckPoint - *ticTocTimerBegin) / CLOCKS_PER_SEC;
	std::cout << this->ticTocTimerElapsedSecs << std::endl;

	delete ticTocTimerCheckPoint;
	return;
}