#include "lungModelling.h"

#ifndef _CELLARRAY_
#define _CELLARRAY_



class cellArray{
public:
	std::vector<std::vector<std::vector<std::vector<int>>>> cells;
	std::vector<std::vector<int>> VertexPerCell;
	void initialize(std::vector<std::vector<float>> &verts, float dx);
};










#endif _CELLARRAY_