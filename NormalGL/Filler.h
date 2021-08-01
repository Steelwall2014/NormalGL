#pragma once
#include "Graphic.h"
class Filler
{
private:
	static void __floodfill(int x, int y, Color oldcolor, Color newcolor);
public:
	// ·ºÌî³ä£¨Õ»£©
	static void StackFloodFill(Point2D seed, Color fill_color);
	// ·ºÌî³ä£¨µÝ¹é£©
	static void RecurseFloodFill(Point2D seed, Color fill_color);
	// É¨ÃèÏßÖÖ×Óµã
	static void ScanlineFill(Point2D seed, Color fill_color);
};