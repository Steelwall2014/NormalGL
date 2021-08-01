#pragma once
#ifndef _centercut_H
#define _centercut_H
#include <string>
#include "Graphic.h"
#include "GeoDefine.h"

class Painter;
class Clipper
{
private:
	vector<Point2D> boundary;
	bool boundary_is_rect = false;
public:
	Clipper() {}
	Clipper(Point2D window_lower_left, Point2D window_upper_right);
	Clipper(vector<Point2D> &boundary);
	bool CenterLineClip(Painter *painter, vector<Painter *> &container);	// 中点裁剪
	void SutherlandHodgmanClip(Painter *painter, Painter *new_painter);
	int WeilerAthertonClip(Painter *painter, vector<Painter *> &container);
	// 设置裁剪多边形
	void SetBoundary(Point2D window_lower_left, Point2D window_upper_right);
	void SetBoundary(vector<Point2D> &boundary);
	void Reset();
	vector<Point2D> GetBoundary() { return boundary; }
};


#endif