#include "model.h"
#include <QOpenGLShaderProgram>
#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QFile>
#include <QTextStream>
#include <QDebug>

/*******************************************************************************
  Implementation of Model methods
  ******************************************************************************/

Model::Model(QString objFilename, QOpenGLShaderProgram* program, VIEW type)
    : QObject(),
    m_vertexBuf(QOpenGLBuffer::VertexBuffer),
    m_normalBuf(QOpenGLBuffer::VertexBuffer),
    m_colorBuf(QOpenGLBuffer::VertexBuffer),
    m_indexBuf(QOpenGLBuffer::IndexBuffer),
    m_bboxBuf(QOpenGLBuffer::VertexBuffer),
    m_bboxindexBuf(QOpenGLBuffer::IndexBuffer),
    m_offset(0)
{
    initializeOpenGLFunctions();

    // Set the world transformation matrix
    setWorldMatrix(type);

    // Create a vertex array object. In OpenGL ES 2.0 and OpenGL 2.x
    // implementations this is optional and support may not be present
    // at all. Nonetheless the below code works in all cases and makes
    // sure there is a VAO when one is needed.
    m_vao.create();
    QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);

    // Generate 4 VBOs for Model
    m_vertexBuf.create();
    m_normalBuf.create();
    m_colorBuf.create();
    m_indexBuf.create();

    loadOBJ(objFilename);

    setupVertexAttribs();

    // Tell OpenGL which VBOs to use
    m_vertexBuf.bind();
    // Tell OpenGL programmable pipeline how to locate vertex position data
    int vertexLocation = program->attributeLocation("vertexPosition_modelspace");
    program->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3, sizeof(QVector3D));
    program->enableAttributeArray(vertexLocation);
    m_vertexBuf.release();

    // Tell OpenGL which VBOs to use
    m_normalBuf.bind();
    // Tell OpenGL programmable pipeline how to locate vertex position data
    int normalLocation = program->attributeLocation("vertexNormal_modelspace");
    program->setAttributeBuffer(normalLocation, GL_FLOAT, 0, 3, sizeof(QVector3D));
    program->enableAttributeArray(normalLocation);
    m_normalBuf.release();

    // Tell OpenGL which VBOs to use
    m_colorBuf.bind();
    // Tell OpenGL programmable pipeline how to locate vertex texture coordinate data
    int colorLocation = program->attributeLocation("vertexColor");
    program->setAttributeBuffer(colorLocation, GL_FLOAT, 0, 3, sizeof(QVector3D));
    program->enableAttributeArray(colorLocation);
    m_colorBuf.release();

    // Generate 2 VBOs for Bbox
    m_bboxvao.create();
    QOpenGLVertexArrayObject::Binder bboxvaoBinder(&m_bboxvao);
    m_bboxBuf.create();
    m_bboxindexBuf.create();

    setupBboxVertexAttribs();
    //*
    // Tell OpenGL which VBOs to use
    m_bboxBuf.bind();
    // Tell OpenGL programmable pipeline how to locate vertex position data
    program->setAttributeBuffer(vertexLocation, GL_FLOAT, 0, 3, 2 * sizeof(QVector3D));
    program->enableAttributeArray(vertexLocation);
    program->disableAttributeArray(normalLocation);
    program->setAttributeBuffer(colorLocation, GL_FLOAT, sizeof(QVector3D), 3, 2 * sizeof(QVector3D));
    program->enableAttributeArray(colorLocation);
    m_bboxBuf.release();
    //*/
}

Model::~Model()
{
    cleanup();
}

void Model::cleanup()
{
    // Destroy 4 VBOs
    m_vertexBuf.destroy();
    m_normalBuf.destroy();
    m_colorBuf.destroy();
    m_indexBuf.destroy();
    // Destroy 2 VBOs
    m_bboxBuf.destroy();
    m_bboxindexBuf.destroy();
}

void Model::loadOBJ(const QString& objFilename)
{
    m_objFilename = objFilename;
    loadOBJ();
    calculateBbox();
    //setRandomCorrelations();
}

void Model::loadNewOBJ(dotObj * input)
{
	loadOBJ(input);
	calculateBbox();
	setupVertexAttribs();
	// No need setupBboxVertexAttribs(); becouse the number of bbox verts are always the same 
	updateAll();
}

void Model::loadNewOBJ(const QString& objFilename)
{
    m_objFilename = objFilename;
    loadOBJ();
    calculateBbox();
    setupVertexAttribs();
    // No need setupBboxVertexAttribs(); becouse the number of bbox verts are always the same 
    updateAll();
}

void Model::saveOBJ(const QString& objFilename)
{
    qDebug() << "Exporting OBJ...";

    QFile file(objFilename);
    if (file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        QTextStream out(&file);
        out << "# " + (objFilename) << endl;
        out.setRealNumberNotation(QTextStream::FixedNotation);
        out.setRealNumberPrecision(6);
		hsv colorHSV;
		rgb colorRGB;


		//hsv colorBuffer;
		//colorBuffer.h = 0;
		//colorBuffer.s = 0;
		//colorBuffer.v = 0;

		int g = 0;

        for (int i = 0; i < m_vertices.size(); i++)
        {
			

			colorHSV.h = 360*m_colors[i].x();
			colorHSV.s = m_colors[i].y();
			colorHSV.v = m_colors[i].z();

			//if ((colorHSV.h != colorBuffer.h) || (colorHSV.s != colorBuffer.s) || (colorHSV.v != colorBuffer.v))
			//{
			//	g++;
			//	colorBuffer = colorHSV;
			//	out << "o object_"<< g << endl;
			//}

			colorRGB = hsv2rgb(colorHSV);
			
			out << "v" << " " << m_vertices[i].x() << " " << m_vertices[i].y() << " " << m_vertices[i].z() << " " << colorRGB.r << " " << colorRGB.g << " " << colorRGB.b << endl;
        }


        for (int i = 0; i < m_normals.size(); i++)
        {
            out << "vn" << " " << m_normals[i].x() << " " << m_normals[i].y() << " " << m_normals[i].z() << endl;
        }



        //out << "s 1" << endl;
        for (int i = 0; i < m_indices.size(); i+=3)
        {
			out << "f" << " " << m_indices[i] + 1 << "//" << m_indices[i] + 1 << " " << m_indices[i + 1] + 1 << "//" << m_indices[i + 1] + 1 << " " << m_indices[i + 2] + 1 << "//" << m_indices[i + 2] + 1 << endl;
        }
        file.close();
    }

    qDebug() << "OBJ Export Complete";
}

void Model::draw(QOpenGLShaderProgram* program, int wireframe, bool bbox)
{
    // Bind the program
    program->bind();

    // Bind the models VAO
    QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);
    
    // Tell OpenGL which VBOs to use
    m_vertexBuf.bind();
    m_normalBuf.bind();
    m_colorBuf.bind();
    m_indexBuf.bind();

    if (wireframe == Qt::Unchecked)
    {
        program->setUniformValue("isWireframe", (GLfloat)0.0);
        //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        // Draw cube geometry using indices from VBO 3
        glDrawElements(GL_TRIANGLES, m_indices.count(), GL_UNSIGNED_INT, 0);
    }
    else if (wireframe == Qt::PartiallyChecked)
    {
        program->setUniformValue("isWireframe", (GLfloat)0.0);
        //glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        // Draw cube geometry using indices from VBO 3
        glDrawElements(GL_TRIANGLES, m_indices.count(), GL_UNSIGNED_INT, 0);
        program->setUniformValue("isWireframe", (GLfloat)1.0);
        //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        // Draw cube geometry using indices from VBO 3
        glDrawElements(GL_TRIANGLES, m_indices.count(), GL_UNSIGNED_INT, 0);
    }
    else if (wireframe == Qt::Checked)
    {
        program->setUniformValue("isWireframe", (GLfloat)0.0);
        //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        // Draw cube geometry using indices from VBO 3
        glDrawElements(GL_TRIANGLES, m_indices.count(), GL_UNSIGNED_INT, 0);
    }

    // Release the buffers
    m_vertexBuf.release();
    m_normalBuf.release();
    m_colorBuf.release();
    m_indexBuf.release();

    if (bbox)
    {
        // Bind the Bbox VAO
        QOpenGLVertexArrayObject::Binder vaoBinder(&m_bboxvao);

        // Tell OpenGL which VBOs to use
        m_bboxBuf.bind();
        m_bboxindexBuf.bind();

        glDrawElements(GL_LINES, m_bboxindices.count(), GL_UNSIGNED_INT, 0);

        m_bboxBuf.release();
        m_bboxindexBuf.release();
    }
    // Release the program
    program->release();
}

/*
void Model::populateColorsPerVertex(int targetPointID)
{
    if (!m_colors.empty()) m_colors.clear();
    //qDebug() << "Colors per vertex populating.." << m_colors.size();
    // Give color based on y value in modelspace
    for (int i = 0; i < m_vertices.size(); i++)
    {
        m_colors.append(QVector3D(m_vertices[targetPointID].distanceToPoint(m_vertices[i])/100., 1, 0.5));
		
    }
    m_colors[targetPointID].setY(0);
    m_colors[targetPointID].setZ(0);
    //qDebug() << "Colors populated sucessfully" << m_colors.count();
    updateColors();
}
*/

void Model::populateColorsPerVertexList(std::vector<int> list)
{
	for (int i = 0; i < list.size(); i++)
	{
		m_colors[list.at(i)] = QVector3D(0.2, 1, 0.5); 
	}
	updateColors();
}

void Model::populateColorsPerVertexWithSegmentProperty(std::vector<int> & segmentPropertyMapVertex)
{

	std::vector<int>::iterator resultMax, resultMin;
	resultMax = std::max_element(segmentPropertyMapVertex.begin(), segmentPropertyMapVertex.end());
	resultMin = std::min_element(segmentPropertyMapVertex.begin(), segmentPropertyMapVertex.end());
	int max = segmentPropertyMapVertex.at(std::distance(segmentPropertyMapVertex.begin(), resultMax));
	int min = segmentPropertyMapVertex.at(std::distance(segmentPropertyMapVertex.begin(), resultMin));

	std::vector<int> list;
	std::cout << "Min is  " << min << " Max is  " << max << std::endl;
	for (int j = 0; j < max; j++){
		if (list.size() > 0){ list.clear(); }
		for (int i = 0; i < segmentPropertyMapVertex.size(); i++){
			if (segmentPropertyMapVertex.at(i) == j){
				list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(0) - 1);
				list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(3) - 1);
				list.push_back(lungmodel.at(lungmodel.size() - 1).faces.at(i).at(6) - 1);
			}
		}
		this->populateColorsPerVertexList(list, j*(1.0 / (double)max));
	}
	return;
}





void Model::populateColorsPerVertexList(std::vector<int> list, double HSVValue)
{

	if (((HSVValue) > 1) || ((HSVValue) < 0)){HSVValue = 0.5;}

	if (m_colors.empty()){
		for (int i = 0; i < m_vertices.size(); i++)
		{
			m_colors.append(QVector3D(0.8, 1, 0.5));
		}
	}


	for (int i = 0; i < list.size(); i++)
	{
		m_colors[list.at(i)] = QVector3D(HSVValue, 1, 0.5);
	}

	updateColors();
}







void Model::populateColorsPerVertexUniform(double HSVValue)
{
	if (!m_colors.empty()) m_colors.clear();
	if (m_colors.empty()){
		for (int i = 0; i < m_vertices.size(); i++)
		{
			m_colors.append(QVector3D(HSVValue, 1, 0.5));
		}
	}

	updateColors();
}


void Model::populateColorsWithCorrelations(int targetPointID)
{
    if (!m_colors.empty()) m_colors.clear();

    qDebug() << "Colors per vertex populating.." << targetPointID;

    // Give color based on y value in modelspace
    int i = 0;
    for (; i < targetPointID; i++)
    {
        m_colors.append(QVector3D(m_correlations[targetPointID][i], 1, 0.5));
    }

    // Black color for the selected point
    m_colors.append(QVector3D(0, 0, 0));
    i++;

    for (; i < m_vertices.size(); i++)
    {
        m_colors.append(QVector3D(m_correlations[i][targetPointID], 1, 0.5));
    }
    qDebug() << "Colors populated sucessfully" << m_vertices.count() << m_colors.count();
    updateColors();
}

void Model::populateColorsWithBar()
{
    if (!m_colors.empty()) m_colors.clear();

    qDebug() << "Colors with Bar populating..";

    // Give color based on the y coordinate in modelspace
    for (int i = 0; i < m_vertices.size(); i++)
    {
        m_colors.append(QVector3D((m_vertices[i].y() + m_offset + 100.)/200., 1, 0.5));
    }
    qDebug() << "Colors populated sucessfully" << m_colors.count();
    updateColors();
}

void Model::loadOBJ()
{
    qDebug() << "Model " << m_objFilename << " is Loading...";
    QFile file(m_objFilename);
    if (file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        // Clear previous data
        if (m_vertices.size()) m_vertices.clear();
        if (m_normals.size()) m_normals.clear();
        if (m_colors.size()) m_colors.clear();
        if (m_indices.size()) m_indices.clear();


        //qDebug() << m_vertices.size() << m_normals.size() << m_colors.size() << m_indices.size()/3;
        QTextStream in(&file);
        in.setRealNumberPrecision(6);
        QVector<QVector3D> tempNormals;
        while (!in.atEnd())
        {
            QString line = in.readLine();
            // If line is empty go to the next one
            if (line.isEmpty()) continue;
            // Split line in tokens
            QStringList tokens = line.split(" ", QString::SkipEmptyParts);
            //*
            if (tokens[0] == "v")
            {
                m_vertices.append(QVector3D(tokens[1].toFloat(),
                    tokens[2].toFloat(),
                    tokens[3].toFloat()));

                if (tokens.size() == 7)
                    m_colors.append(QVector3D(tokens[4].toFloat(), tokens[5].toFloat(), tokens[6].toFloat()));
                else
                    m_colors.append(QVector3D(0.8f, 0.8f, 0.8f));
            }
            else if (tokens[0] == "vn")
            {
                tempNormals.append(QVector3D(tokens[1].toFloat(), tokens[2].toFloat(), tokens[3].toFloat()));
            }
            else if (tokens[0] == "f")
            {
                QStringList tok1, tok2, tok3;
                tok1 = tokens[1].split("/", QString::KeepEmptyParts);
                tok2 = tokens[2].split("/", QString::KeepEmptyParts);
                tok3 = tokens[3].split("/", QString::KeepEmptyParts);

                int v1 = tok1[0].toInt() - 1;
                int v2 = tok2[0].toInt() - 1;
                int v3 = tok3[0].toInt() - 1;

                m_indices.append(v1);
                m_indices.append(v2);
                m_indices.append(v3);
            }
        }
        file.close();
        qDebug() << m_vertices.size() << m_normals.size() << m_colors.size() << m_indices.size() / 3;
        if (m_vertices.size() == tempNormals.size()) m_normals = tempNormals;
        else calculateNormals();

        qDebug() << "Model " << m_objFilename << " Loaded successfully!!" << m_vertices.count() << m_normals.count() << m_colors.count() << m_indices.count() / 3;
    }
    else
    {
        qDebug() << "Unable to open file: " << m_objFilename;
        file.close();
    }
}



void Model::calculateNormals()
{
    m_normals.fill(QVector3D(), m_vertices.size());
    qDebug() << "Calculating Normals...";

    qDebug() << __LINE__ << m_vertices.size() << m_normals.size() << m_colors.size() << m_indices.size() / 3;

    for (int i = 0; i < m_indices.size() / 3; i++)
    {
        // Non normalized triangle normal
        QVector3D norm = QVector3D::crossProduct((m_vertices[m_indices[3*i + 1]] - m_vertices[m_indices[3*i]]),
                                                 (m_vertices[m_indices[3*i + 2]] - m_vertices[m_indices[3*i]]));
        m_normals[m_indices[3*i]] += norm;
        m_normals[m_indices[3*i + 1]] += norm;
        m_normals[m_indices[3*i + 2]] += norm;
    }
    // Normalize per vertex normals
    for (int i = 0; i < m_normals.size(); i++)
    {
        m_normals[i].normalize();
    }
    qDebug() << "Normals Calculated sucessfully" << m_normals.size();
}

void Model::calculateBbox()
{
    QVector3D max;
    QVector3D min(m_vertices[0]);
    for (int i = 0; i < m_vertices.size(); i++)
    {
        if (max.x() < m_vertices[i].x()) max.setX(m_vertices[i].x());
        if (max.y() < m_vertices[i].y()) max.setY(m_vertices[i].y());
        if (max.z() < m_vertices[i].z()) max.setZ(m_vertices[i].z());
        if (min.x() > m_vertices[i].x()) min.setX(m_vertices[i].x());
        if (min.y() > m_vertices[i].y()) min.setY(m_vertices[i].y());
        if (min.z() > m_vertices[i].z()) min.setZ(m_vertices[i].z());
    }
    QVector3D cen((max + min) / 2);
    QVector3D rad((max - min) / 2);
    // Convert to world space
    m_bbox.set(m_worldMatrix * cen, m_worldMatrix * rad);
    qDebug() << "BBox: Center =" << cen << "Max =" << rad;

    // Clear previous data
    if (m_bboxdata.size()) m_bboxdata.clear();
    if (m_bboxindices.size()) m_bboxindices.clear();

    // Set up the vertex data
    // 0
    m_bboxdata.append(cen - rad);
    m_bboxdata.append(QVector3D(1, 0, 1));
    // 1
    m_bboxdata.append(cen + QVector3D(rad.x(), -rad.y(), -rad.z()));
    m_bboxdata.append(QVector3D(1, 0, 1));
    // 2
    m_bboxdata.append(cen + QVector3D(rad.x(), -rad.y(), rad.z()));
    m_bboxdata.append(QVector3D(1, 0, 1));
    // 3
    m_bboxdata.append(cen + QVector3D(-rad.x(), -rad.y(), rad.z()));
    m_bboxdata.append(QVector3D(1, 0, 1));
    // 4
    m_bboxdata.append(cen + QVector3D(-rad.x(), rad.y(), -rad.z()));
    m_bboxdata.append(QVector3D(1, 0, 1));
    // 5
    m_bboxdata.append(cen + QVector3D(rad.x(), rad.y(), -rad.z()));
    m_bboxdata.append(QVector3D(1, 0, 1));
    // 6
    m_bboxdata.append(cen + rad);
    m_bboxdata.append(QVector3D(1, 0, 1));
    // 7
    m_bboxdata.append(cen + QVector3D(-rad.x(), rad.y(), rad.z()));
    m_bboxdata.append(QVector3D(1, 0, 1));
   
    m_bboxindices.append(0);
    m_bboxindices.append(1);
    // 1-2
    m_bboxindices.append(1);
    m_bboxindices.append(2);
    // 2-3
    m_bboxindices.append(2);
    m_bboxindices.append(3);
    // 3-0
    m_bboxindices.append(3);
    m_bboxindices.append(0);
    // 4-5
    m_bboxindices.append(4);
    m_bboxindices.append(5);
    // 5-6
    m_bboxindices.append(5);
    m_bboxindices.append(6);
    // 6-7
    m_bboxindices.append(6);
    m_bboxindices.append(7);
    // 7-4
    m_bboxindices.append(7);
    m_bboxindices.append(4);
    // 0-4
    m_bboxindices.append(0);
    m_bboxindices.append(4);
    // 1-5
    m_bboxindices.append(1);
    m_bboxindices.append(5);
    // 2-6
    m_bboxindices.append(2);
    m_bboxindices.append(6);
    // 3-7
    m_bboxindices.append(3);
    m_bboxindices.append(7);
}

void Model::setWorldMatrix(VIEW type)
{
    // Set the world transformation matrix
    m_worldMatrix.setToIdentity();
    switch (type)
    {
    case VIEW::FRONT:
        break;
    case VIEW::BACK:
        m_worldMatrix.rotate(180, 1, 0, 0);
        break;
    case VIEW::RIGHT:
        m_worldMatrix.rotate(-90, 0, 1, 0);
        break;
    case VIEW::LEFT:
        m_worldMatrix.rotate(90, 0, 1, 0);
        break;
    case VIEW::UP:
        m_worldMatrix.rotate(90, 1, 0, 0);
        break;
    case VIEW::DOWN:
        m_worldMatrix.rotate(-90, 0, 1, 0);
        m_worldMatrix.rotate(-90, 1, 0, 0);
        break;
    }
}

QMatrix4x4& Model::worldMatrix()
{ 
    return m_worldMatrix;
}

void Model::intersectWithRay(Ray& r, int& closest, QList<int>& triList)
{
    QList<float> tauList;
    QList<int> closestList;

    float tau;
    int closestPoinID;
    for (int i = 0; i < m_indices.size() / 3; i++)
    {
        if (intersectTriWithRay(r, i, tau, closestPoinID))
        {
            triList.append(i);
            tauList.append(tau);
            closestList.append(closestPoinID);
        }
    }

    // Sort the List using quick bubblesort
    int n = triList.size();
    if (n != 0)
    {
        do
        {
            int newn = 0;
            for (int i = 1; i < n; i++)
            {
                if (tauList[i - 1] > tauList[i])
                {
                    tauList.swap(i - 1, i);
                    triList.swap(i - 1, i);
                    closestList.swap(i - 1, i);
                    newn = i;
                }
            }
            n = newn;
        } while (n != 0);
        closest = closestList[0];
    }
}

bool Model::intersectTriWithRay(Ray& r, int& triNum, float& tau, int& closest)
{
    const QVector3D& e = r.start();
    const QVector3D& d = r.dir();

    QVector3D &a = m_vertices[m_indices[3 * triNum]];
    QVector3D &b = m_vertices[m_indices[3 * triNum + 1]];
    QVector3D &c = m_vertices[m_indices[3 * triNum + 2]];

    // Solve the 3x3 equation e + t*d = a + gi*(b - a) + vi*(c - a)
    float xaxb = a.x() - b.x();
    float yayb = a.y() - b.y();
    float zazb = a.z() - b.z();
    float xaxc = a.x() - c.x();
    float yayc = a.y() - c.y();
    float zazc = a.z() - c.z();
    float xd = d.x();
    float yd = d.y();
    float zd = d.z();
    float xaxe = a.x() - e.x();
    float yaye = a.y() - e.y();
    float zaze = a.z() - e.z();

    float D = xaxb * (yayc * zd - yd   * zazc)
        - xaxc * (yayb * zd - yd   * zazb)
        + xd * (yayb * zazc - yayc * zazb);

    float Dg = xaxb * (yaye * zd - yd   * zaze)
        - xaxe * (yayb * zd - yd   * zazb)
        + xd   * (yayb * zaze - yaye * zazb);
    float vi = Dg / D;
    if (vi < 0 || vi > 1) return false;

    float Db = xaxe * (yayc * zd - yd   * zazc)
        - xaxc * (yaye * zd - yd   * zaze)
        + xd   * (yaye * zazc - yayc * zaze);
    float gi = Db / D;
    if (gi < 0 || gi > 1 - vi) return false;

    // There are three medians described in barycentric coordinates
    // gi = vi, gi = -2 * vi + 1, gi = -0.5 * vi + 0.5
    if (gi < vi)
    {
        if (gi < -2 * vi + 1) closest = m_indices[3 * triNum];
        else closest = m_indices[3 * triNum + 2];
    }
    else
    {
        if (gi < -0.5 * vi + 0.5) closest = m_indices[3 * triNum];
        else closest = m_indices[3 * triNum + 1];
    }

    float Dt = xaxb * (yayc * zaze - yaye * zazc)
        - xaxc * (yayb * zaze - yaye * zazb)
        + xaxe * (yayb * zazc - yayc * zazb);

    tau = Dt / D;


    return true;
}

void Model::setupVertexAttribs()
{
    // Transfer vertex data to VBO 0
    m_vertexBuf.bind();
    m_vertexBuf.allocate(m_vertices.constData(), m_vertices.count() * sizeof(QVector3D));
    m_vertexBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    m_vertexBuf.release();

    // Transfer vertex data to VBO 1
    m_normalBuf.bind();
    m_normalBuf.allocate(m_normals.constData(), m_normals.count() * sizeof(QVector3D));
    m_normalBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    m_normalBuf.release();

    // Transfer vertex data to VBO 2
    m_colorBuf.bind();
    m_colorBuf.allocate(m_colors.constData(), m_colors.count() * sizeof(QVector3D));
    m_colorBuf.release();

    // Transfer index data to VBO 3
    m_indexBuf.bind();
    m_indexBuf.allocate(m_indices.constData(), m_indices.count()* sizeof(GLuint));
    m_indexBuf.release();
}

void Model::setupBboxVertexAttribs()
{
    // Transfer vertex data to VBO 0
    m_bboxBuf.bind();
    m_bboxBuf.allocate(m_bboxdata.constData(), m_bboxdata.count() * sizeof(QVector3D));
    m_bboxBuf.setUsagePattern(QOpenGLBuffer::StaticDraw);
    m_bboxBuf.release();
    //*/

    // Transfer index data to VBO 1
    m_bboxindexBuf.bind();
    m_bboxindexBuf.allocate(m_bboxindices.constData(), m_bboxindices.count()* sizeof(GLuint));
    m_bboxindexBuf.release();
}

BBox& Model::bbox()
{ 
    return m_bbox;
}

void Model::updateVertices()
{
    QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);
    m_vertexBuf.bind();
    m_vertexBuf.write(0, m_vertices.constData(), m_vertices.count() * sizeof(QVector3D));
    m_vertexBuf.release();
    emit verticesChanged();
}

void Model::updateNormals()
{
    QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);
    m_normalBuf.bind();
    m_normalBuf.write(0, m_normals.constData(), m_normals.count() * sizeof(QVector3D));
    m_normalBuf.release();
    emit normalsChanged();
}

void Model::updateColors()
{
    QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);
    m_colorBuf.bind();
    m_colorBuf.write(0, m_colors.constData(), m_colors.count() * sizeof(QVector3D));
    m_colorBuf.release();
    emit colorsChanged();
}

void Model::updateAll()
{
    QOpenGLVertexArrayObject::Binder vaoBinder(&m_vao);
    // Update the Model
    m_vertexBuf.bind();
    m_vertexBuf.write(0, m_vertices.constData(), m_vertices.count() * sizeof(QVector3D));
    m_vertexBuf.release();

    m_normalBuf.bind();
    m_normalBuf.write(0, m_normals.constData(), m_normals.count() * sizeof(QVector3D));
    m_normalBuf.release();

    m_colorBuf.bind();
    m_colorBuf.write(0, m_colors.constData(), m_colors.count() * sizeof(QVector3D));
    m_colorBuf.release();

    // allocate the correct size
    m_indexBuf.bind();
    m_indexBuf.allocate(m_indices.constData(), m_indices.count()* sizeof(GLuint));
    m_indexBuf.release();

    QOpenGLVertexArrayObject::Binder vaoBboxBinder(&m_bboxvao);
    // Update the Bbox
    m_bboxBuf.bind();
    m_bboxBuf.write(0, m_bboxdata.constData(), m_bboxdata.count() * sizeof(QVector3D));
    m_bboxBuf.release();
    // No need to update bbox indices size it is always the same
}

void Model::setRandomCorrelations()
{
    if (m_correlations.size()) m_correlations.clear();
    m_correlations.reserve(m_vertices.size());

    qDebug() << "Initializing Correlation Matrix..";
    for (int i = m_vertices.size(); i > 0; i--)
    {
        QVector<float> temp(i);
        for (int j = 0; j < i; j++)
        {
            temp[j] = m_vertices[i-1].distanceToPoint(m_vertices[j]) / 100.;
            //qDebug() << __LINE__;
        }
        //qDebug() << temp.size();
        m_correlations.prepend(temp);
    }
    //for (int i = 0; i < m_correlations.size(); i++)
    //{
    //    qDebug() << m_correlations[i].size();
    //}
    qDebug() << "Correlations created" << m_correlations.size();
}


void Model::loadOBJ(dotObj * input)
{
	qDebug() << "Model " << m_objFilename << " is Loading...";
	QFile file(m_objFilename);
	QVector<QVector3D> tempNormals;
	
		// Clear previous data
		if (m_vertices.size()) m_vertices.clear();
		if (m_normals.size()) m_normals.clear();
		if (m_colors.size()) m_colors.clear();
		if (m_indices.size()) m_indices.clear();


		for (unsigned int i = 0; i < input->vertices.size(); i++){

			m_vertices.append(QVector3D(input->vertices.at(i).at(0),input->vertices.at(i).at(1),input->vertices.at(i).at(2)));
			m_colors.append(QVector3D(0.8f, 0.8f, 0.8f));
		}
		for (unsigned int i = 0; i < input->faces.size(); i++){

			m_indices.append(input->faces.at(i).at(0)-1);
			m_indices.append(input->faces.at(i).at(3)-1);
			m_indices.append(input->faces.at(i).at(6)-1);
		}
		for (unsigned int i = 0; i < input->normals.size(); i++){
			tempNormals.append(QVector3D(input->normals.at(i).at(0), input->normals.at(i).at(1), input->normals.at(i).at(2)));
		}

		


		if (m_vertices.size() == tempNormals.size()) 
		{ 
			m_normals = tempNormals; 
		}
		else
		{
		calculateNormals();
		}
		
		
		qDebug() << "Model " << m_objFilename << " Loaded successfully!!" << m_vertices.count() << m_normals.count() << m_colors.count() << m_indices.count() / 3;
	
		return;
}






hsv rgb2hsv(rgb in)
{
	hsv         out;
	double      min, max, delta;

	min = in.r < in.g ? in.r : in.g;
	min = min  < in.b ? min : in.b;

	max = in.r > in.g ? in.r : in.g;
	max = max  > in.b ? max : in.b;

	out.v = max;                                // v
	delta = max - min;
	if (delta < 0.00001)
	{
		out.s = 0;
		out.h = 0; // undefined, maybe nan?
		return out;
	}
	if (max > 0.0) { // NOTE: if Max is == 0, this divide would cause a crash
		out.s = (delta / max);                  // s
	}
	else {
		// if max is 0, then r = g = b = 0              
		// s = 0, v is undefined
		out.s = 0.0;
		out.h = NAN;                            // its now undefined
		return out;
	}
	if (in.r >= max)                           // > is bogus, just keeps compilor happy
		out.h = (in.g - in.b) / delta;        // between yellow & magenta
	else
		if (in.g >= max)
			out.h = 2.0 + (in.b - in.r) / delta;  // between cyan & yellow
		else
			out.h = 4.0 + (in.r - in.g) / delta;  // between magenta & cyan

	out.h *= 60.0;                              // degrees

	if (out.h < 0.0)
		out.h += 360.0;

	return out;
}


rgb hsv2rgb(hsv in)
{
	double      hh, p, q, t, ff;
	long        i;
	rgb         out;

	if (in.s <= 0.0) {       // < is bogus, just shuts up warnings
		out.r = in.v;
		out.g = in.v;
		out.b = in.v;
		return out;
	}
	hh = in.h;
	if (hh >= 360.0) hh = 0.0;
	hh /= 60.0;
	i = (long)hh;
	ff = hh - i;
	p = in.v * (1.0 - in.s);
	q = in.v * (1.0 - (in.s * ff));
	t = in.v * (1.0 - (in.s * (1.0 - ff)));

	switch (i) {
	case 0:
		out.r = in.v;
		out.g = t;
		out.b = p;
		break;
	case 1:
		out.r = q;
		out.g = in.v;
		out.b = p;
		break;
	case 2:
		out.r = p;
		out.g = in.v;
		out.b = t;
		break;

	case 3:
		out.r = p;
		out.g = q;
		out.b = in.v;
		break;
	case 4:
		out.r = t;
		out.g = p;
		out.b = in.v;
		break;
	case 5:
	default:
		out.r = in.v;
		out.g = p;
		out.b = q;
		break;
	}
	return out;
}
