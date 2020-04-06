#ifndef _Model_H_
#define _Model_H_

#include "ray.h"
#include "bbox.h"
#include <QObject>
#include <QList>
#include <QVector>
#include <QVector3D>
#include <QString>
#include <QOpenGLFunctions>
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QMatrix4x4>
#include "modellerToViewer.h"




/*******************************************************************************
  Model
  ******************************************************************************/

class Model : public QObject, protected QOpenGLFunctions
{
Q_OBJECT

public:
    enum VIEW
    {
        FRONT = 1,
        BACK = 1 << 1,
        RIGHT = 1 << 2,
        LEFT = 1 << 3,
        UP = 1 << 4,
        DOWN = 1 << 5
    };

    // Ctor
    Model(QString objFilename, QOpenGLShaderProgram* program, VIEW type = VIEW::FRONT);
    ~Model();

    // Methods
    void loadOBJ(const QString& objFilename);
    void loadNewOBJ(const QString& objFilename);
	void loadNewOBJ(dotObj * input);
    void saveOBJ(const QString& objFilename);
    void draw(QOpenGLShaderProgram* program, int wireframe = 0, bool bbox = false);

    void populateColorsPerVertex(int targetPointID);
    void populateColorsWithCorrelations(int targetPointID);
    void populateColorsWithBar();
	void populateColorsPerVertexList(std::vector<int> list);
	void populateColorsPerVertexList(std::vector<int> list, double HSVValue);
	void populateColorsPerVertexWithSegmentProperty(std::vector<int> & segmentPropertyMapVertex);
	void populateColorsPerVertexUniform(double HSVValue);

    void setWorldMatrix(VIEW type);
    QMatrix4x4& worldMatrix();
    void intersectWithRay(Ray& r, int& closest, QList<int>& triList);

    BBox& bbox();
    void setRandomCorrelations();

public:
    // Data members
    QList< QVector<float> > m_correlations;
    int m_offset;

private:
    // Methods
    void loadOBJ();
	void loadOBJ(dotObj * input);
    void calculateNormals();
    void calculateBbox();
    bool intersectTriWithRay(Ray& r, int& triNum, float& tau, int& closest);
    void setupVertexAttribs();
    void setupBboxVertexAttribs();
    void updateVertices();
    void updateNormals();
    void updateColors();
    void updateAll();

private:
    // Data members
    QString m_objFilename;
    QVector<QVector3D> m_vertices;
    QVector<QVector3D> m_normals;
    QVector<QVector3D> m_colors;
    QVector<GLuint> m_indices;

    BBox m_bbox;
    QVector<QVector3D> m_bboxdata;
    QVector<GLuint> m_bboxindices;

    // Transformation matrix that converts local coordinates to world coordinates
    QMatrix4x4 m_worldMatrix;

    // OpenGL stuff
    QOpenGLVertexArrayObject m_vao;
    QOpenGLBuffer m_vertexBuf;
    QOpenGLBuffer m_normalBuf;
    QOpenGLBuffer m_colorBuf;
    QOpenGLBuffer m_indexBuf;

    QOpenGLVertexArrayObject m_bboxvao;
    QOpenGLBuffer m_bboxBuf;
    QOpenGLBuffer m_bboxindexBuf;

    //int m_selectedPoint;

public slots:
    void cleanup();

signals:
    void verticesChanged();
    void normalsChanged();
    void colorsChanged();
};






typedef struct {
	double r;       // percent
	double g;       // percent
	double b;       // percent
} rgb;

typedef struct {
	double h;       // angle in degrees
	double s;       // percent
	double v;       // percent
} hsv;

static hsv   rgb2hsv(rgb in);
static rgb   hsv2rgb(hsv in);


#endif
