#include <QMenuBar>
#include <QMenu>
#include <QMessageBox>
#include <QFileDialog>
#include <QDialog>
#include <QTimer>
#include <QDebug>
#include "modellerToViewer.h"


QString buffer; 
bool	brushSelectionFunctionality = true;
bool	partSelectionFunctionality = false;
int		_brushSizeBox = 200;
float _brushDistance=1.0;
std::string _brushdistanceMode="Nearest Neighbours";
int _brushdistanceModeIndex = 0;
dotObj *tempModel;
std::vector<dotObj> lungmodel; 
int	 seedPoint;


