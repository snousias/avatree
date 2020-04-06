#ifndef _Camera_H_
#define _Camera_H_

#include "bbox.h"
#include <QObject>
#include <QVector3D>
#include <QQuaternion>
#include <QMatrix4x4>

class Camera : public QObject
{
Q_OBJECT

public:
    // Ctor
    Camera();
    ~Camera();

public:
    // Methods
    const QVector3D& eye();
    const QVector3D& focal();
    const QVector3D& upvector();
    const QVector3D& translation();
    const QMatrix4x4& projectionMatrix();
    const QMatrix4x4& viewMatrix();
    qreal distance();
    qreal sensitivity();

    void translate(QVector3D  tran);
    void rotate(QQuaternion  rot);
	void calibrate(BBox box);
    void setAspect(qreal aspect);
    void setDistance(int wheel);

private:
    // Data members
    // Perspective
    qreal m_fov;
    qreal m_aspect;
    qreal m_zNear;
    qreal m_zFar;
    // View
    QVector3D m_eye;
    QVector3D m_focal;
    QVector3D m_upvec;
    // Transforamtions
    QVector3D m_translation;
    QQuaternion m_rotation;
    // Matrices
    QMatrix4x4 m_projectionMatrix;
    QMatrix4x4 m_viewMatrix;
    // Interaction
    qreal m_cDist;
    qreal m_zoomStep;
    qreal m_sensitivity;

private:
    // Update the corresponding Matrices
    void updateView();
    void updateProjection();

signals:
    void updated();
};

#endif
