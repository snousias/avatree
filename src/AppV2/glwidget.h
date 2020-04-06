#ifndef _GLWIDGET_H_
#define _GLWIDGET_H_


#include <QPixmap>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QInputEvent>
#include <QKeyEvent>
#include <QWheelEvent>
#include <QPoint>
#include <QVector2D>
#include <QVector3D>
#include <QQuaternion>
#include <QMatrix4x4>
#include <QString>
#include <QBasicTimer>
#include "camera.h"
#include "model.h"
/*******************************************************************************
  GLWidget
  ******************************************************************************/

class GLWidget : public QOpenGLWidget, protected QOpenGLFunctions
{
Q_OBJECT

public:
    // ctor
    GLWidget(QWidget* parent = 0);
    ~GLWidget();
    
    // Methods
    QSize minimumSizeHint() const Q_DECL_OVERRIDE;
    QSize sizeHint() const Q_DECL_OVERRIDE;
    Model* getModel(){ return m_model; }
    void saveScreenshot();
    void castRay(const QPoint& pos, Ray& r);
    int m_wireframe;
    bool m_bbox;
    // The Timer
    QBasicTimer m_timer;
    void loadNewOBJ(QString& filename);
	void loadNewOBJ(dotObj *input);
    void saveOBJ(QString& filename);

	//protected:

    // OpenGL framework Methods 
    void initializeGL() Q_DECL_OVERRIDE;
    void paintGL() Q_DECL_OVERRIDE;
	
    void resizeGL(int width, int height) Q_DECL_OVERRIDE;
    // IO handle Methods
    void mousePressEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
    void mouseMoveEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
    //void mouseReleaseEvent(QMouseEvent *event) Q_DECL_OVERRIDE;
    void wheelEvent(QWheelEvent *event) Q_DECL_OVERRIDE;
    void keyPressEvent(QKeyEvent *event) Q_DECL_OVERRIDE;
    void keyReleaseEvent(QKeyEvent *event) Q_DECL_OVERRIDE;
    int  mkModif(QInputEvent *event);
    // Timer Event handle
    //void timerEvent(QTimerEvent *e) Q_DECL_OVERRIDE;

private:
    // Methods
    void initShaders();
    void initModel(QString& filename, Model::VIEW type = Model::VIEW::FRONT);
    

public:
    // Screen pixmap
    QPixmap m_originalPixmap;
    // Shader program
    QOpenGLShaderProgram m_program;
    // Camera
    Camera m_camera;
    // Model
    Model* m_model;
    // Transformation Matrices
    QMatrix4x4 m_mvMatrix;
    // Window Properties
    QVector2D m_lastPos;
    int m_lastClosestPointID;
    bool m_transparent;
    
public slots:
    void cleanup();

signals:
    void wireframeModeChanged(Qt::CheckState state);
    void bboxModeChanged(bool drawBbox);
};

#endif
