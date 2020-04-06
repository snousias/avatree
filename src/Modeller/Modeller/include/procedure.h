#include "lungModelling.h"

#ifndef _PROCEDURECLASS_
#define _PROCEDURECLASS_

class procedure{
public:
	int seedPoint;
	int _brushSizeBox;
	bool	brushSelectionFunctionality;
	bool	partSelectionFunctionality;
	std::vector<int> currentSelection;
	dotObj *tempModel;
	std::vector<dotObj> lungmodel;
	std::string id;
	
};

#endif _PROCEDURECLASS_