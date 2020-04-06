#include "lungModelling.h"


#ifndef _PROJECT_
#define _PROJECT_

class project{

public:
std::string token;
std::string buffer;

int seedPoint;
int startingVertex; //same as seedpoint
int _brushdistanceModeIndex;
int _brushSizeBox;
float _brushDistance;

bool	brushSelectionFunctionality;
bool	partSelectionFunctionality;

std::vector<int> currentSelection;
double narrowingRatio;
double localDIameterBeforeNarrowing;
double localDIameterAfterNarrowing;


dotObj *tempModel;
std::vector<dotObj> lungmodel;
dotObj leftLungBoundary;
dotObj rightLungBoundary;

status stat;

void update(void);

};

#endif