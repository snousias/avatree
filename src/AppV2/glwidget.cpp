#include  "glwidget.h"
#include  "ray.h"
#include  "modellerToViewer.h"
#include  <QtWidgets>
#include  <QWidget>
#include  <QMouseEvent>
#include  <QOpenGLShaderProgram>
#include  <QCoreApplication>
#include  <math.h>


GLWidget::GLWidget(QWidget *parent)
	: QOpenGLWidget(parent),
	m_model(NULL),
	m_wireframe(0),
	m_bbox(false),
	m_transparent(true)

{
	setFocusPolicy(Qt::WheelFocus);

	if (m_transparent) setAttribute(Qt::WA_TranslucentBackground);

	//connect(&m_camera, SIGNAL(updated()), this, SLOT(update()));
	connect(&m_camera, &Camera::updated, this, static_cast<void (GLWidget::*)()>(&GLWidget::update));
}

GLWidget::~GLWidget()
{
	cleanup();
}

void GLWidget::cleanup()
{
	// Make sure the context is current when deleting the buffers.
	makeCurrent();
	delete m_model;
	m_model = NULL;
	doneCurrent();
}

QSize GLWidget::minimumSizeHint() const
{
	return QSize(50, 50);
}

QSize GLWidget::sizeHint() const
{
	return QSize(800, 600);
}

void GLWidget::saveScreenshot()
{
	qDebug() << "Save screenshot...";
	// Capture
	qApp->beep();
	m_originalPixmap = QPixmap(); // clear image for low memory situations on embedded devices.
	QScreen *screen = QGuiApplication::primaryScreen();
	if (screen)
		m_originalPixmap = screen->grabWindow(0);

	// Save
	QString format = "png";
	QString initialPath = QDir::currentPath() + tr("/untitled.") + format;

	QString fileName = QFileDialog::getSaveFileName(this, tr("Save As"), initialPath,
		tr("%1 Files (*.%2);;All Files (*)")
		.arg(format.toUpper())
		.arg(format));
	if (!fileName.isEmpty())
		m_originalPixmap.save(fileName, format.toLatin1().constData());
	qDebug() << "Screenshot saved";
}

void GLWidget::castRay(const QPoint& pos, Ray& r)
{
	int w = width();
	int h = height();
	// Constuct ray in camera space

	// Convert selected pixel from screen space to normalized [-1, 1] cliping space
	QVector3D target_clipingspace(((pos.x() + 0.5) / (0.5 * w)) - 1,
		-(((pos.y() + 0.5) / (0.5 * h)) - 1),
		0);
	// Convert target to camera space
	QVector3D target_cameraspace = m_camera.projectionMatrix().inverted() * target_clipingspace;

	// Convert ray from camera to model space
	QVector3D eye_modelspace = m_mvMatrix.inverted() * QVector3D(0, 0, 0);
	QVector3D target_modelspace = m_mvMatrix.inverted() * target_cameraspace;

	// Get the ray in model space
	r.set(eye_modelspace, target_modelspace);
}

void GLWidget::initializeGL()
{
	// In this example the widget's corresponding top-level window can change
	// several times during the widget's lifetime. Whenever this happens, the
	// QOpenGLWidget's associated context is destroyed and a new one is created.
	// Therefore we have to be prepared to clean up the resources on the
	// aboutToBeDestroyed() signal, instead of the destructor. The emission of
	// the signal will be followed by an invocation of initializeGL() where we
	// can recreate all resources.
	connect(context(), &QOpenGLContext::aboutToBeDestroyed, this, &GLWidget::cleanup);

	initializeOpenGLFunctions();

	// White backround
	glClearColor(0, 0, 0, m_transparent ? 0 : 1);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);

	// Cull triangles which normal is not towards the camera
	glEnable(GL_CULL_FACE);

	glLineWidth(1);

	initShaders();



	QString filename = QString(":/Models/cube.obj");
	std::string modelstring = filename.toStdString();
	dotObj theModel;
	theModel.initializeFromFile(modelstring);
	lungmodel.push_back(theModel);


	initModel(filename, Model::VIEW::FRONT);
	tempModel = new dotObj;

	// Our camera changes in this example.
	m_camera.calibrate(m_model->bbox());

	GLfloat dist = m_camera.distance();

	// Bind the program
	m_program.bind();

	// Light position is fixed.
	m_program.setUniformValue("lightPosition_cameraspace", QVector3D(0, 0, 0));
	m_program.setUniformValue("lightColor", QVector3D(1, 1, 1));
	m_program.setUniformValue("lightPower", dist * dist);
	m_program.setUniformValue("isBlack", 0);

	m_program.release();

	// Use QBasicTimer because its faster than QTimer
	//m_timer.start(12, this);
}

void GLWidget::paintGL()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// Get the current world Matrix transformation from model
	QMatrix4x4 worldMatrix = m_model->worldMatrix();
	// Get the current view Matrix transformation from camera
	QMatrix4x4 viewMatrix = m_camera.viewMatrix();
	// Calculate modelview Matrix
	m_mvMatrix = viewMatrix * worldMatrix;

	m_program.bind();
	m_program.setUniformValue("projMatrix", m_camera.projectionMatrix());
	m_program.setUniformValue("mvMatrix", m_mvMatrix);
	m_program.setUniformValue("viewMatrix", viewMatrix);
	m_program.setUniformValue("modelMatrix", worldMatrix);
	m_program.setUniformValue("normalMatrix", m_mvMatrix.normalMatrix());

	m_model->draw(&m_program, m_wireframe, m_bbox);

	m_program.release();
}

void GLWidget::resizeGL(int w, int h)
{
	// Calculate aspect ratio
	qreal aspect = qreal(w) / qreal(h ? h : 1);
	// Set perspective projection
	m_camera.setAspect(aspect);
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
	std::vector<int> currentSelection;

	if ((event->buttons() & Qt::RightButton) && !(event->modifiers() & Qt::ControlModifier))
	{
		Ray ray_modelspace;
		castRay(event->pos(), ray_modelspace);

		// Get a list with intersecting triangles
		int closestPointID = m_lastClosestPointID;
		QList<int> triainter;
		m_model->intersectWithRay(ray_modelspace, closestPointID, triainter);

		if (triainter.size() && closestPointID != m_lastClosestPointID)
		{
			m_lastClosestPointID = closestPointID;
			seedPoint = closestPointID;
			qDebug() << "Last clicked vertex : " << seedPoint;

			

			if (partSelectionFunctionality)
			{
				std::vector<int> v = lungmodel.at(lungmodel.size() - 1).segment_property_map_per_vertex;
				if (v.size() > 0){
					int group = v.at(seedPoint);
					lungmodel.at(lungmodel.size() - 1).selectedVertices.clear();
					for (int i = 0; i < v.size(); i++){
						if (v.at(i) == group){
							lungmodel.at(lungmodel.size() - 1).selectedVertices.push_back(i);
						}
					}
				}
			}

			if (brushSelectionFunctionality)
			{
				currentSelection = lungmodel.at(lungmodel.size() - 1).selectedVertices;
				lungmodel.at(lungmodel.size() - 1).selectedVertices.clear();

				if (_brushdistanceModeIndex==0)
				{
					lungmodel.at(lungmodel.size() - 1).selector(seedPoint, _brushSizeBox);
				}
				if (_brushdistanceModeIndex==1)
				{
					lungmodel.at(lungmodel.size() - 1).selectorBasedOnCenterline(seedPoint, _brushDistance);
				}
				lungmodel.at(lungmodel.size() - 1).selectedVertices.insert(lungmodel.at(lungmodel.size() - 1).selectedVertices.end(), currentSelection.begin(), currentSelection.end());
				sort(lungmodel.at(lungmodel.size() - 1).selectedVertices.begin(), lungmodel.at(lungmodel.size() - 1).selectedVertices.end());
				lungmodel.at(lungmodel.size() - 1).selectedVertices.erase(unique(lungmodel.at(lungmodel.size() - 1).selectedVertices.begin(), lungmodel.at(lungmodel.size() - 1).selectedVertices.end()), lungmodel.at(lungmodel.size() - 1).selectedVertices.end());
			}
			
			m_model->populateColorsPerVertexUniform(0.8);
			if (lungmodel.at(lungmodel.size() - 1).segment_property_map.size() > 0){
				m_model->populateColorsPerVertexWithSegmentProperty(lungmodel.at(lungmodel.size() - 1).segment_property_map);
			}
			m_model->populateColorsPerVertexList(lungmodel.at(lungmodel.size() - 1).selectedVertices);

			//m_model->populateColorsPerVertex(closestPointID);
			//m_model->populateColorsWithCorrelations(closestPointID);

		}
		else {
			lungmodel.at(lungmodel.size() - 1).selectedVertices.clear();
			m_model->populateColorsPerVertexUniform(0.8);
			if (lungmodel.at(lungmodel.size() - 1).segment_property_map.size() > 0){
				m_model->populateColorsPerVertexWithSegmentProperty(lungmodel.at(lungmodel.size() - 1).segment_property_map);
			}
			qDebug() << "No intersection";
		}
	}

	

	m_lastPos = QVector2D(event->localPos());
}


void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
	// Mouse release position - mouse press position
	QVector2D diff = QVector2D(event->localPos()) - m_lastPos;

	if (event->buttons() & Qt::MiddleButton){
		m_camera.translate(QVector3D(diff.x()/10, -diff.y()/10, 0));

	}

	if (event->buttons() & Qt::LeftButton)
	{
		// Rotation axis is perpendicular to the mouse position difference vector
		QVector3D rotationAxis(diff.y(), diff.x(), 0.0);
		// Accelerate angular speed relative to the length of the mouse sweep
		qreal angle = diff.length() / 10.0;

		// Update rotation
		m_camera.rotate(QQuaternion::fromAxisAndAngle(rotationAxis.normalized(), angle));
	}
	if (event->buttons() & Qt::RightButton)
	{
		Ray ray_modelspace;
		castRay(event->pos(), ray_modelspace);

		// Get a list with intersecting triangles
		int closestPointID;
		QList<int> triainter;
		m_model->intersectWithRay(ray_modelspace, closestPointID, triainter);

		if (triainter.size())
		{
			if (closestPointID != m_lastClosestPointID)
			{
				m_lastClosestPointID = closestPointID;
				seedPoint = closestPointID;
				lungmodel.at(lungmodel.size() - 1).selector(seedPoint, _brushSizeBox);
				m_model->populateColorsPerVertexList(lungmodel.at(lungmodel.size() - 1).selectedVertices);
			}
		}
		else qDebug() << "No intersection";
	}

	m_lastPos = QVector2D(event->localPos());
}


void GLWidget::wheelEvent(QWheelEvent *event)
{
	m_camera.setDistance(event->delta());
}

void GLWidget::keyPressEvent(QKeyEvent *event)
{
	int modif = mkModif(event);
	QString c = event->text();
	unsigned char key = c.toLatin1()[0];
	switch (isprint(key) ? tolower(key) : key)
	{
	case 'p':
		saveScreenshot();
		break;
	case 'w':
	{
		m_wireframe = ++m_wireframe % 3;
		emit wireframeModeChanged((Qt::CheckState) m_wireframe);
		update();
	}
	break;
	case 'b':
	{
		m_bbox = !m_bbox;
		emit bboxModeChanged(m_bbox);
		update();
	}
	break;
	case 'c':
		m_camera.calibrate(m_model->bbox());
		break;
	case 'l':
	{
		QString format = "obj";
		QString initialPath = QDir::currentPath() + tr("/untitled.") + format;

		QString fileName = QFileDialog::getOpenFileName(this, tr("Load Model"), initialPath,
			tr("%1 Files (*.%2);;All Files (*)").arg(format.toUpper()).arg(format));

		if (!fileName.isEmpty())
			loadNewOBJ(fileName);
		// Our camera changes in this example.
		m_camera.calibrate(m_model->bbox());
	}
	break;
	case 's':
	{
		QString format = "obj";
		QString initialPath = QDir::currentPath() + tr("/untitled.") + format;

		QString fileName = QFileDialog::getSaveFileName(this, tr("Save Model"), initialPath,
			tr("%1 Files (*.%2);;All Files (*)").arg(format.toUpper()).arg(format));

		if (!fileName.isEmpty())
			saveOBJ(fileName);
		// Our camera changes in this example.
		m_camera.calibrate(m_model->bbox());
	}
	break;
	}

	switch (event->key())
	{
	case Qt::Key_Up:
		m_camera.translate(QVector3D(0, 1, 0));
		break;
	case Qt::Key_Down:
		m_camera.translate(QVector3D(0, -1, 0));
		break;
	case Qt::Key_Right:
		m_camera.translate(QVector3D(1, 0, 0));
		break;
	case Qt::Key_Left:
		m_camera.translate(QVector3D(-1, 0, 0));
		break;
	case Qt::Key_Escape:
		qApp->quit();
		break;
	}
}

void GLWidget::keyReleaseEvent(QKeyEvent *event)
{}

int GLWidget::mkModif(QInputEvent *event)
{
	int ctrl = event->modifiers() & Qt::ControlModifier ? 1 : 0;
	int shift = event->modifiers() & Qt::ShiftModifier ? 1 : 0;
	int alt = event->modifiers() & Qt::AltModifier ? 1 : 0;
	int modif = (ctrl << 0) | (shift << 1) | (alt << 2);
	return modif;
}

/*
void GLWidget::timerEvent(QTimerEvent *)
{
// Decrease angular speed (friction)
m_angularSpeed *= 0.99;

// Stop rotation when speed goes below threshold
if (m_angularSpeed < 0.01) {
m_angularSpeed = 0.0;
}
else {
// Update rotation
m_rotation = QQuaternion::fromAxisAndAngle(m_rotationAxis, m_angle) * m_rotation;

// Request an update
update();
}
}
//*/

void GLWidget::initShaders()
{
	if (!m_program.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/vertShader.glsl"))
	{
		GLenum err = glGetError();
		qDebug() << "OpenGL ERROR No:" << err;
		close();
	}
	if (!m_program.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/fragShader.glsl"))
	{
		GLenum err = glGetError();
		qDebug() << "OpenGL ERROR No:" << err;
		close();
	}
	if (!m_program.link())
	{
		GLenum err = glGetError();
		qDebug() << "OpenGL ERROR No:" << err;
		close();
	}
	if (!m_program.bind())
	{
		GLenum err = glGetError();
		qDebug() << "OpenGL ERROR No:" << err;
		close();
	}
	else m_program.release();
}

void GLWidget::initModel(QString& filename, Model::VIEW type)
{
	if (m_model)
	{
		qDebug() << "Model already exist!! Wrong use of this API";
		// Disconnect previous signals from slots
		//disconnect(m_model, SIGNAL(verticesChanged()), this, SLOT(update()));
		//disconnect(m_model, SIGNAL(normalsChanged()), this, SLOT(update()));
		//disconnect(m_model, SIGNAL(colorsChanged()), this, SLOT(update()));
		delete m_model;
		m_model = NULL;
	}

	// Load the model
	m_model = new Model(filename, &m_program, type);

	// Connect new signals to slots
	//connect(m_model, SIGNAL(verticesChanged()), this, SLOT(update()));
	//connect(m_model, SIGNAL(normalsChanged()), this, SLOT(update()));
	//connect(m_model, SIGNAL(colorsChanged()), this, SLOT(update()));
	connect(m_model, &Model::verticesChanged, this, static_cast<void (GLWidget::*)()>(&GLWidget::update));
	connect(m_model, &Model::normalsChanged, this, static_cast<void (GLWidget::*)()>(&GLWidget::update));
	connect(m_model, &Model::colorsChanged, this, static_cast<void (GLWidget::*)()>(&GLWidget::update));
}

void GLWidget::loadNewOBJ(QString& filename)
{
	// Load the model
	m_model->loadNewOBJ(filename);
}

void GLWidget::loadNewOBJ(dotObj * input)
{
	// Load the model
	m_model->loadNewOBJ(input);
}

void GLWidget::saveOBJ(QString& filename)
{
	// Load the model
	m_model->saveOBJ(filename);
}
