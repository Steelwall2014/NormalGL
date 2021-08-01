#pragma once

#include "Graphic.h"
class Raster
{
public:
	static void textOut( const char* pText , int x, int y , Color color );
	static void drawCell(int row, int col , unsigned color );
	static void drawGrid();
};

