#include <QMenuBar>
#include <QMenu>
#include <QMessageBox>
#include <QFileDialog>
#include <QDialog>
#include <QTimer>
#include <QDebug>
#include "lungModelling.h"


#ifndef _MODELLERTOVIEWER_H_
#define _MODELLERTOVIEWER_H_

extern QString buffer; 

extern int seedPoint;
extern int _brushSizeBox;
extern float _brushDistance;
extern std::string _brushdistanceMode;
extern int _brushdistanceModeIndex;
extern bool	brushSelectionFunctionality;
extern bool	partSelectionFunctionality;


extern dotObj *tempModel;
extern std::vector<dotObj> lungmodel;

#endif