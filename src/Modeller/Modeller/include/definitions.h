#ifndef _DEFS_
#define _DEFS_








#if _MSC_VER
#include <Windows.h>
#include <conio.h>
#endif
#include <ctime>
#include <math.h>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <thread> 
#include <filesystem>
#include <vector>
#include <algorithm>
#include <stdio.h> 
#include <Eigen/Eigenvalues>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/SparseQR>
#include <Eigen/SparseLU>
#include <Eigen/OrderingMethods>
#include <Eigen/Geometry>
#include <boost/math/distributions/normal.hpp>
#include <regex> 
#include <boost/thread.hpp>
#include <unordered_map> 

#define EXPMAXNUMOFBRANCH 10000.0
#define PRECISION 5

#ifdef _OPENMP
#include <omp.h>
#define MAX_THREADS 4
#else
//int omp_get_thread_num(void) { return 0; }
//int omp_get_num_threads(void) { return 1; }
#endif


#ifdef _BULLET_
#include "softDemo.h"
#include "GlutStuff.h"
#include "GLDebugDrawer.h"
#include "btBulletDynamicsCommon.h"
#include "hacdCircularList.h"
#include "hacdVector.h"
#include "hacdICHull.h"
#include "hacdGraph.h"
#include "hacdHACD.h"
#include "cd_wavefront.h"
#include "ConvexBuilder.h"
#include "btBulletDynamicsCommon.h"
#include "LinearMath/btQuickprof.h"
#include "LinearMath/btIDebugDraw.h"
#include "LinearMath/btGeometryUtil.h"
#include "BulletCollision/CollisionShapes/btShapeHull.h"
#include "GLDebugDrawer.h"
#include "btBulletDynamicsCommon.h"
#include "BulletSoftBody/btSoftRigidDynamicsWorld.h"
#include "BulletCollision/CollisionDispatch/btSphereSphereCollisionAlgorithm.h"
#include "BulletCollision/NarrowPhaseCollision/btGjkEpa2.h"
#include "LinearMath/btQuickprof.h"
#include "LinearMath/btIDebugDraw.h"
#include "BunnyMesh.h"
#include "TorusMesh.h"
#include "LinearMath/btConvexHull.h"
#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"
#include "BulletSoftBody/btSoftBodyHelpers.h"
#include "SoftDemo.h"
#include "GL_ShapeDrawer.h"
#include "GLDebugFont.h"
#include "GlutStuff.h"
#endif
//#include "armadillo"


/*=====================================*/
//SOFTBODY
#define APPLIEDFORCEONGLOTTIS 250
#define timeToModelExtraction 200
#define KLST 0.01f
#define KAST 0.0001f
#define KVST 0.001f
#define KDP 0.8
#define KMT 0.05
#define VERTICALRESTRICTION 0.1		//
#define HORIZONTALRESTRICTION 0.4	//
#define SETOFPOINTS 3000

//LAPLACIAN MESH CONTRACTION
#define SKELETONIZATIONPOINTS 1000 //1000 initially
#define WEIGHTER 0.65
#define LINEAR true
#define SKELETONIZATIONITERATIONS 8
#define Pi 3.1415926
#define FREQ Pi

#define MODELSDIRECTORY "../../../../models/"
#define OUTPUTDITECTORY "C:\\Users\\Stavros\\Desktop\\"




using namespace Eigen;
//using namespace arma;
class status;
class PCA;
class dotObj;
class multiBodyObj;
class ray;
class volume;
class glottis;
class world;
class sorter;
class plane;
class bbox;
class treeGeneration;
class helper;
class simulation;
class cellArray;
class procedure;
class analysis;
class ggraph;
class gpath;
class gbranch;
class gline;
class gnode;
class distance;


class mDiscretization;
class mIndex;
class mEdge;
class mPosition;
class mVertex;
class mFace;
class mCell;
class mGrid;



#define _POISSON_MODE_2


#endif
