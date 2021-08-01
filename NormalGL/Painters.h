#pragma once
#include <string>

#include "ogr_geometry.h"
//#include "..\\3rdParty\\gdal\\include\\ogr_geometry.h"
#include "Graphic.h"
#include "Utils.h"
#include "GeoDefine.h"
#include "Matrix.h"
#include "Tin.h"
using namespace std;


enum DisplayType;
class Model;
class Painter
{
/// <summary>
/// 图形绘制的基类
/// </summary>
/// 作者：张径舟
protected:
	string __shape = "BASE";
	bool __is_closed_geom = false;

	virtual void __paintedge(PixelPoint* points, int point_count, Color edge_color) { return; };	// 在各个图形中实现
	virtual void __paintedge(int x1, int y1, int x2, int y2, Color color) { return; };	// 在各个图形中实现
public:
	Color edge_color = _RGB(0, 0, 0), fill_color = _RGB(100, 100, 100);
	Point2D *points;
	int point_count;
	bool visible = true;
	//Filler *filler = nullptr;
	Matrix transformation = Matrix::CreateIdentityMatrix(3);
	virtual bool Paint() { return false; }
	virtual bool Paint(Matrix &trans) { return false; }
	virtual Point2D *AdjustToDisplayType(DisplayType dtype, Point2D *points, int point_count);
	//virtual bool Fill(vector<SeedPoint> seed_points) { return false; };
	//virtual bool Fill() { return false; };
	string GetPaintShape() { return this->__shape; }
	static Painter *CreatePainter();
	virtual double GetX(int i) { return 0.0; }
	virtual double GetY(int i) { return 0.0; }
	virtual void SetX(int i, double x) { return; }
	virtual void SetY(int i, double y) { return; }
	virtual int GetPointCount() { return 0; }
	virtual Point2D *GetPoints() { return nullptr; }
	virtual void ClearData() { delete[] points; }
	//bool IsClosedGeom() { return filler? true: false; }
};



class RasterPainter
{
/// <summary>
///	栅格绘制类。（存在bug，颜色显示不对劲，暂时只支持BSQ格式）
/// </summary>
/// 作者：张径舟
/// 20210121：已弃用
private:
	Image* image;
public:
	void AddImage(Image* img);
	void AddImage(string filename);
	bool Paint(int bands[]);
};



class GridPainter : public Painter
{
/// <summary>
///	格网类。这个类比较特殊，只在切换到格网模式时调用Paint
/// </summary>
/// 作者：张径舟
public:
	GridPainter()
	{
		__shape = "GRID";
	}
	bool Paint();
};



class LinePainter : public Painter
{
/// <summary>
/// 线的绘制类
/// </summary>
/// 作者：汤子豪，张径舟
/// 整合：张径舟
private:
	//Point2D points[2];
	void __paintline(int x, int y, int offset_x, int offset_y, Color color, float k);
protected:
	void __paint(int x1, int y1, int x2, int y2, Color color);
public:
	LinePainter(const LinePainter &other)
	{
		__shape = "LINE";
		this->points = new Point2D[2];
		this->point_count = 2;
		this->points[0] = other.points[0];
		this->points[1] = other.points[1];
		this->transformation = other.transformation;
	}
	LinePainter()
	{
		__shape = "LINE";
		//double x1, y1, x2, y2
		this->points = new Point2D[2];
		this->point_count = 2;
		getRubberPoints(points[0].x, points[0].y, points[1].x, points[1].y);
	}
	LinePainter(double x1, double y1, double x2, double y2)
	{
		__shape = "LINE";
		this->points = new Point2D[2];
		this->point_count = 2;
		points[0].x = x1; points[0].y = y1; points[1].x = x2; points[1].y = y2;
	}
	double GetX(int i) { return points[i].x; }
	double GetY(int i) { return points[i].y; }
	void SetX(int i, double x) { points[i].x = x; }
	void SetY(int i, double y) { points[i].y = y; }
	int GetPointCount() { return 2; }
	bool Paint();
	bool Paint(Matrix &);
	Point2D *GetPoints();
	//void ClearData();
};

class PolylinePainter :public LinePainter
{
/// <summary>
/// 折线的绘制
/// </summary>
/// 作者：张径舟
//private:
//	int point_count;
//	Point2D* points;
public:
	PolylinePainter(const PolylinePainter &other)
	{
		__shape = "POLYLINE";
		this->point_count = other.point_count;
		this->points = new Point2D[other.point_count];
		for (int i = 0; i < point_count; i++)
			this->points[i] = other.points[i];
		this->transformation = other.transformation;
	}
	PolylinePainter()
	{
		__shape = "POLYLINE";
		point_count = getRubberPointCount();
		//PixelPoint *pts = new PixelPoint[point_count];
		points = new Point2D[point_count];
		getRubberPoints(points);
		//for (int i = 0; i < point_count; i++)
		//	points[i] = Point2D{ (double)pts[i].x, (double)pts[i].y };
		//delete [point_count]pts;
		fill_color = get_random_color();
	}
	PolylinePainter(Point2D *points, int point_count)
	{
		__shape = "POLYLINE";
		this->points = new Point2D[point_count];
		for (int i = 0; i < point_count; i++)
			this->points[i] = points[i];
		this->point_count = point_count;
		fill_color = get_random_color();
	}
	~PolylinePainter()
	{
		delete[] points;
		//delete filler;
	}
	double GetX(int i) { return points[i].x; }
	double GetY(int i) { return points[i].y; }
	void SetX(int i, double x) { points[i].x = x; }
	void SetY(int i, double y) { points[i].y = y; }
	int GetPointCount() { return point_count; }
	bool Paint();
	bool Paint(Matrix &);
	Point2D *GetPoints();
};





class CirclePainter : public Painter
{
/// <summary>
/// 圆的绘制类
/// </summary>
/// 作者：汤子豪
/// 整合：张径舟
private:
	//Point2D points[2];
	void __paintedge(int x1, int y1, int x2, int y2, Color color);
	void drawthu(int ox, int oy, int x, int y, Color color);
public:
	CirclePainter()
	{
		this->points = new Point2D[2];
		this->point_count = 2;
		this->transformation = Matrix::CreateIdentityMatrix(3);
		getRubberPoints(points[0].x, points[0].y, points[1].x, points[1].y);//获取橡皮线点坐标
		//filler = new Filler(this->edge_color);
		__shape = "CIRCLE";
	}
	CirclePainter(double x1, double y1, double x2, double y2)
	{
		this->points = new Point2D[2];
		this->point_count = 2;
		this->transformation = Matrix::CreateIdentityMatrix(3);
		points[0].x = x1; points[0].y = y1; points[1].x = x2; points[1].y = y2;
		__shape = "CIRCLE";
	}
	bool Paint();
	bool Paint(Matrix &);
};

	
class EllipsePainter : public Painter
{
/// <summary>
/// 椭圆的绘制类
/// </summary>
/// 作者：刘航
/// 整合：张径舟
private:
	//Point2D points[2];
	void __paintedge(int x1, int y1, int x2, int y2, Color color);
	void drawthu(int x, int y, int a, int b, Color color);
public:
	EllipsePainter()
	{
		this->points = new Point2D[2];
		this->point_count = 2;
		getRubberPoints(points[0].x, points[0].y, points[1].x, points[1].y);
		//filler = new Filler(this->edge_color);
		__shape = "ELLIPSE";
	}
	bool Paint();
};


class PolygonPainter : public PolylinePainter
{
/// <summary>
/// 多边形的绘制类
/// </summary>
/// 作者：张径舟
protected:
	//Point2D* points;
	//int point_count;
	void __paintedge(Point2D* points, int point_count, Color edge_color);
public:
	PolygonPainter()
	{	
		__shape = "POLYGON";
		__is_closed_geom = true;
		//filler = new Filler(this->edge_color);
		//point_count = 5;
		//points = new Point2D[6];
		//points[0] = Point2D{ 170, 456 };
		//points[1] = Point2D{ 169, 266 };
		//points[2] = Point2D{ 317, 215 };
		//points[3] = Point2D{ 268, 346 };
		//points[4] = Point2D{ 225, 299 };
		//points[5] = Point2D{ 170, 456 };

		point_count = getRubberPointCount() + 1;
		points = new Point2D[point_count];
		getRubberPoints(points);
		points[point_count-1] = points[0];
	}
	PolygonPainter(Point2D *points, int point_count)
	{
		__shape = "POLYGON";
		__is_closed_geom = true;
		//filler = new Filler(this->edge_color);
		this->points = new Point2D[point_count];
		for (int i = 0; i < point_count; i++)
			this->points[i] = points[i];
		this->point_count = point_count;
	}
	~PolygonPainter()
	{
		delete[] points;
	}
	bool Paint();
	bool Paint(Matrix &);
	Point2D *GetPoints();
};

class RectanglePainter : public PolygonPainter
{
private:
	//Point2D points[5];
	double min_x, max_x, min_y, max_y;
public:
	RectanglePainter()
	{
		this->points = new Point2D[5];
		getRubberPoints(points[0].x, points[0].y, points[1].x, points[1].y);//获取橡皮线点坐标
		__shape = "RECTANGLE";
		min_x = min(points[0].x, points[1].x);
		max_x = max(points[0].x, points[1].x);
		min_y = min(points[0].y, points[1].y);
		max_y = max(points[0].y, points[1].y);
		point_count = 5;
		points[0].x = min_x; points[0].y = min_y;
		points[1].x = min_x; points[1].y = max_y;
		points[2].x = max_x; points[2].y = max_y;
		points[3].x = max_x; points[3].y = min_y;
		points[4].x = min_x; points[4].y = min_y;
	}
	bool Paint();
	bool Paint(Matrix &);
	Point2D *GetPoints();
};

class TrianglePainter :public Painter
{
public:
	int point_count = 3;
	int trishape = 0;			//平底为1，平顶为2，都不是为0
	//int spvertex = 0;			//三角形为平底或者平顶时，存储不是底边的顶点索引；当其不为时，存储中间顶点（用于分割三角形）的顶点索引
	//int flatline[2] = { 0 };	//平顶或平底存储平边的两个顶点索引，拼接的为上下两个顶点
	Point3D spvertex;
	Point3D flatline[2];
	Model *model = nullptr;
	//Point3D operator[](int index);
	TrianglePainter() {}
	//TrianglePainter(Tin *Tin_grid, int pt1_index, int pt2_index, int pt3_index);
	//TrianglePainter(Tin *Tin_grid, int pts_index[3]);
	//TrianglePainter(Point3D p1, Point3D p2, Point3D p3, Model *model);
	TrianglePainter(Point3D p1, Point3D p2, Point3D p3);
	TrianglePainter(Point3D spvertex, Point3D flatline1, Point3D flatline2, int trishape);
	void DetermineTriShape(Point3D p1, Point3D p2, Point3D p3);
	void PaintflatTri();
	void PaintflattopTri();
	void PaintflatbottomTri();
	void PaintJointTri();
	bool Paint();
};


class Line3DPainter :public Painter
{
public:
	int point_count = 2;
	Point3D pt1, pt2;
	Line3DPainter() {}
	//TrianglePainter(Tin *Tin_grid, int pt1_index, int pt2_index, int pt3_index);
	//TrianglePainter(Tin *Tin_grid, int pts_index[3]);
	Line3DPainter(Point3D p1, Point3D p2)
	{
		if (p1.y < p2.y)
		{
			this->pt1 = p1;
			this->pt2 = p2;
		}
		else
		{
			this->pt2 = p1;
			this->pt1 = p2;
		}
	}
	bool Paint();
};