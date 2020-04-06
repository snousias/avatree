#include "lungModelling.h"


void cellArray::initialize(std::vector<std::vector<float>> &verts,float cellsize){

	
	float minx = verts.at(0).at(0);
	float maxx = verts.at(0).at(0);

	float miny = verts.at(0).at(1);
	float maxy = verts.at(0).at(1);

	float minz = verts.at(0).at(2);
	float maxz = verts.at(0).at(2);

	for (int k = 0; k < verts.size(); k++){
		if (minx > verts.at(k).at(0)){ minx = verts.at(k).at(0); }
		if (maxx<verts.at(k).at(0)){ maxx = verts.at(k).at(0); }

		if (miny>verts.at(k).at(1)){ miny = verts.at(k).at(1); }
		if (maxy<verts.at(k).at(1)){ maxy = verts.at(k).at(1); }

		if (minz>verts.at(k).at(2)){ minz = verts.at(k).at(2); }
		if (maxz<verts.at(k).at(2)){ maxz = verts.at(k).at(2); }

	}
	//maxx += 0.02;
	//maxy += 0.02;
   ///maxz += 0.02;

	int xlen, ylen, zlen;
	float dx, dy, dz;

	/*
	xlen = 100;
	ylen = 100;
	zlen = 100;
	dx = (maxx - minx ) / (float)xlen;
	dy = (maxy - miny)  / (float)ylen;
	dz = (maxz - minz)  / (float)zlen;
	*/

	
	dx = cellsize;
	dy = cellsize;
	dz = cellsize;
	xlen = (int)((maxx - minx) / dx);
	ylen = (int)((maxy - miny) / dy);
	zlen = (int)((maxz - minz) / dz);
	
	
	//Init vertex per cell
	this->VertexPerCell.resize(verts.size());
	for (int i = 0; i < verts.size(); ++i) {
		this->VertexPerCell[i].resize(3);
	}



	//Init cell array
	std::vector<std::vector<std::vector<std::vector<int>> > > nullArray3D;
	this->cells = nullArray3D;
	this->cells.resize(xlen);
	for (int i = 0; i < xlen; ++i) {
		this->cells[i].resize(ylen);
		for (int j = 0; j < ylen; ++j)
			this->cells[i][j].resize(zlen);
	}


	int u, v, w;
	for (int k = 0; k < verts.size(); k++){
		u = (int)((verts.at(k).at(0) - minx) / dx);
		v = (int)((verts.at(k).at(1) - miny) / dy);
		w = (int)((verts.at(k).at(2) - minz) / dz);
		
		if (u == xlen){ u = u - 1; }
		if (v == ylen){ v = v - 1; }
		if (w == zlen){ w = w - 1; }
		//std::cout << u << "-" << v << "-" << w << std::endl;

		//Cell with coordinates u,v,w includes the vertex with index k
		this->cells[u][v][w].push_back(k);

		//Vertex with index k is inside the cell with coordinates u,v,w
		this->VertexPerCell.at(k).at(0) = u;
		this->VertexPerCell.at(k).at(1) = v;
		this->VertexPerCell.at(k).at(2) = w;

	}

	return;
}





void dotObj::initCells(float dx){

	this->cells.initialize(this->vertices,dx); //Distrubute points to a grid




	return;
}