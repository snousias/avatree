#ifndef _BBox_H_
#define _BBox_H_

#include <QVector3D>

/*******************************************************************************
BBox
******************************************************************************/

class BBox
{
public:
	// Ctor
	BBox() : m_center(), m_max(){}
	BBox(QVector3D center, QVector3D max) : m_center(center), m_max(max){}
	~BBox(){}

public:
	// Methods
	const QVector3D& center() { return m_center; }
	const QVector3D& max(){return m_max;}
	void set(const QVector3D& center, const QVector3D& max)
	{
		m_center = center;
		m_max = max;
	}

public:
	// Data members
	QVector3D m_center;
	QVector3D m_max;
};

#endif