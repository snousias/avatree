#ifndef _MAINWINDOW_H_
#define _MAINWINDOW_H_

#include "ui_mainwindow.h"
#include "modellerToViewer.h"
#include <QMainWindow>
#include <QBasicTimer>
#include <QMenuBar>
#include <QMenu>
#include <QMessageBox>
#include <QFileDialog>
#include <QDialog>
#include <QTimer>
#include <QDebug>

/*******************************************************************************
  MainWindow
  ******************************************************************************/

class MainWindow : public QMainWindow, public Ui::MainWindow
{
Q_OBJECT

public:
    // ctor
    MainWindow();
    ~MainWindow();

public:
    //QTimer m_timer;

    bool m_animOn;
    int m_timeWindowframe;
    int m_lastFrame;
	QMessageBox msgBox;

private slots:
	void narrow(void);
	void threadtester(void);
	void writeObj(void);
	void exportVerticesToText(void);
	void reload(void);
	void sdf(void);
	void loadSegmentPropertyList(void);
	void exportSegmented(void);
	void exportObj(void);
	void updateBrushType(void);
	void loadFile(void);
	void extendSelection(void);

	//void workspaceFolderPathToLineEdit(void);
	//void existingGeomSetPath(void);
	//void leftLungGeomSetPath(void);
	//void rightLungGeomSetPath(void);
	//void initializeWorkspace(void); 
	void doSimulate(void);
	void updateBrushSize(void);
	void updateDistanceMode(void);
	void updateBrushDistance(void);
	void segmentBySkeleton(void);
	void segmentByGeneration(void);

public:
	std::thread member1Thread() {
		return std::thread([=] { 	threadtester(); });
	}

};

#endif
