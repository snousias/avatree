#include "camera.h"
#include <QVector3D>
#include <QQuaternion>
#include <QMatrix4x4>
#include <QtMath>
#include <QDebug>

/*******************************************************************************
Implementation of Camera methods
  ******************************************************************************/

Camera::Camera()
    : QObject(),
    m_fov(45),
    m_aspect(1),
    m_zNear(1),
    m_zFar(500),
    m_zoomStep(5),
    m_sensitivity(1),
    m_eye(0,0,300),
    m_focal(0,0,0),
    m_upvec(0,1,0)
{}

Camera::~Camera()
{}

const QVector3D& Camera::eye()
{
    return m_eye;
}

const QVector3D& Camera::focal()
{
    return m_focal;
}

const QVector3D& Camera::upvector()
{
    return m_upvec;
}

const QVector3D& Camera::translation()
{
    return m_translation;
}

const QMatrix4x4& Camera::projectionMatrix()
{
    return m_projectionMatrix;
}

const QMatrix4x4& Camera::viewMatrix()
{
    return m_viewMatrix;
}

qreal Camera::distance()
{
    return m_cDist;
}

qreal Camera::sensitivity()
{
    return m_sensitivity;
}

void Camera::translate(QVector3D trans)
{
    m_translation += m_sensitivity * trans;
    updateView();
}

void Camera::rotate(QQuaternion rot)
{
    m_rotation = rot * m_rotation;
    updateView();
}

void Camera::calibrate(BBox box)
{
    m_focal = box.center();
    (qFabs(box.max().x()) > qFabs(box.max().y())) ?
        m_cDist = qFabs(box.max().x()) / qTan(qDegreesToRadians(m_fov / 2)) + qFabs(box.max().z()) :
        m_cDist = qFabs(box.max().y()) / qTan(qDegreesToRadians(m_fov / 2)) + qFabs(box.max().z());
    m_eye = QVector3D(0, 0, 2*m_cDist);//2*
    // Set sensitivity
    m_sensitivity = m_cDist / 100;
    // Reset rotation
    m_rotation = QQuaternion();
    // Reset traslation
    m_translation = QVector3D();
    updateView();
    updateProjection();
}

void Camera::setAspect(qreal aspect)
{
    m_aspect = aspect;
    updateProjection();
}

void Camera::setDistance(int wheel)
{
    m_cDist -= m_zoomStep * m_sensitivity * (wheel / 120.);
    if (m_cDist < 10 * m_sensitivity) m_cDist = 10 * m_sensitivity;
    if (m_cDist > 800 * m_sensitivity) m_cDist = 800 * m_sensitivity;
    m_eye.setZ(2*m_cDist);
    updateView();
}

void Camera::updateView()
{
    m_viewMatrix.setToIdentity();
    // Translate the whole system
    m_viewMatrix.translate(m_translation);
    // Look at the (0,0,0) point
    m_viewMatrix.lookAt(m_eye, QVector3D(0,0,0), m_upvec);
    // Rotate based on quaternion
    m_viewMatrix.rotate(m_rotation);
    // Move to the center of the model
    m_viewMatrix.translate(-m_focal);
    emit updated();
}

void Camera::updateProjection()
{
    // Reset projection
    m_projectionMatrix.setToIdentity();
    // Set perspective projection
    m_projectionMatrix.perspective(m_fov, m_aspect, m_zNear * m_sensitivity, m_zFar * m_sensitivity);
    emit updated();
}
