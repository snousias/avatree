#ifndef _Ray_H_
#define _Ray_H_

#include <QVector3D>

/*******************************************************************************
  Ray
  ******************************************************************************/

class Ray
{
public:
    // Ctor
    Ray() : m_start(), m_end(), m_dir(){}
    Ray(QVector3D& start, QVector3D& end) : m_start(start),
                                            m_end(end),
                                            m_dir((end - start).normalized()){}

    ~Ray(){}

public:
    // Methods
    const QVector3D& start() { return m_start; }
    const QVector3D& end() { return m_end; }
    const QVector3D& dir() { return m_dir; }
    void set(const QVector3D& start, const QVector3D& end)
    {
        m_start = start;
        m_end = end;
        m_dir = (end - start).normalized();
    }

private:
    // Data members
    QVector3D m_start;
    QVector3D m_end;
    QVector3D m_dir;
};

#endif