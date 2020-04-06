#ifndef _OBJEDITOR_
#define _OBJEDITOR_
#include "definitions.h"
#include "cellArray.h"
#include "helper.h"
#include "cgal.h"
#include "obj.h"
#include "dotObj.h"
#include "treeGeneration.h"
#include "simulation.h"
#include "procedure.h"
#include "project.h"
#include "vox.h"


#ifdef _POISSON_MODE_1 
int PoissonRecon(std::string in, std::string out, int depth);
#endif

#ifdef _POISSON_MODE_2
int SSDRecon(std::string in, std::string out, int depth, dotObj * mModel);
#endif




#endif
