#include <string>
#include <stack>
#include <list>
#include <float.h>

#include "Graphic.h"
#include "Status.h"
#include "Utils.h"
#include "Painters.h"
#include "3DModel.h"

//using namespace Painters;
using namespace std;



Point2D *Painter::AdjustToDisplayType(DisplayType dtype, Point2D *points, int point_count)
{
	if (dtype == dtPixelMode)
		return nullptr;
	Point2D *new_points = new Point2D[point_count];
	PixelToGrid(points, new_points, point_count);
	return new_points;
}
Painter *Painter::CreatePainter()
{
	Painter *painter = nullptr;
	switch (status.GetPaintType())
	{
	case ptDrawCircle:
		painter = new CirclePainter();
		break;
	case ptDrawEllipse:
		painter = new EllipsePainter();
		break;
	case ptDrawLine:
		painter = new LinePainter();
		break;
	case ptDrawPolygon:
		painter = new PolygonPainter();
		break;
	case ptDrawPolyline:
		painter = new PolylinePainter();
		break;
	case ptDrawRect:
		painter = new RectanglePainter();
		break;
	}
	return painter;
}



void RasterPainter::AddImage(Image *img)
{
	/*
	* 向RasterPainter中添加一张影像
	* params:
	*	img: Image
	*/
	image = img;
}
void RasterPainter::AddImage(string filename)
{
	/*
	* 向RasterPainter中添加一张影像
	* params:
	*	filename: string，影像的文件位置
	*/
	Image *img = ReadData(filename);
	image = img;
}
bool RasterPainter::Paint(int bands[])
{
	/*
	* 绘制影像
	* params:
	*	bands: 数组，进行合成的波段数
	*/
	int offset = image->filehead->nRow * image->filehead->nCol;	// 后面的每一个波段的偏移量
	int r, g, b, pos;
	for (int row = 0; row < image->filehead->nRow; row++)
	{
		pos = row * image->filehead->nRow;
		for (int col = 0; col < image->filehead->nCol; col++)
		{
			r = image->array[pos + col + offset * (bands[0] - 1)];
			g = image->array[pos + col + offset * (bands[1] - 1)];
			b = image->array[pos + col + offset * (bands[2] - 1)];
			setPixel(col, image->filehead->nRow - row, _RGB(r, g, b));
		}
	}
	return true;
}



bool GridPainter::Paint()
{
	/*
	* 绘制格网
	*/
	winHeight = getWindowHeight();
	winWidth = getWindowWidth();
	int grid_x = 0, grid_y = 0;
	for (int x = 0; x <= winWidth; x += status.GetCellSize())
	{
		for (int y = 0; y <= winHeight; y += status.GetCellSize())
		{
			(grid_x + grid_y) % 2 == 0 ? set_unit(grid_x, grid_y, _RGB(255, 255, 255)) : set_unit(grid_x, grid_y, status.GetGridColor());
			grid_y++;
		}
		grid_x++;
		grid_y = 0;
	}
	return true;
}



void LinePainter::__paintline(int x, int y, int offset_x, int offset_y, Color color, float k)
{
	/*
	* 绘制线的实现
	* params:
	*	x: 开始点的坐标
	*	y: 开始点的坐标
	*	offset_x: 第二个点相对第一个点的偏移
	*	offset_y: 第二个点相对第一个点的偏移
	*	color: 颜色
	*	k: 直线斜率
	* 作者：汤子豪，张径舟
	*/
	int dx, dy, e, x1 = 0, y1 = 0;
	dx = x;
	dy = y;
	e = -dx;
	while (x1 <= x)
	{
		// 这些if else不太好看，但是还没想好怎么处理
		if (0.0 <= k && k <= 1.0)
		{
			set_unit(x1 + offset_x, y1 + offset_y, color);
		}
		else if (k > 1.0)
		{
			set_unit(y1 + offset_x, x1 + offset_y, color);
		}
		else if (-1.0 <= k && k < 0.0)
		{
			set_unit(x1 + offset_x, -y1 + offset_y, color);
		}
		else if (k < -1.0)
		{
			set_unit(y1 + offset_x, -x1 + offset_y, color);
		}
		e += 2 * dy;
		if (e > 0)
		{
			x1++;
			y1++;
			e -= 2 * dx;
		}
		else
		{
			x1++;
		}
	}
}
void LinePainter::__paint(int x1, int y1, int x2, int y2, Color color)
{
	/*
	* 绘制一条线段，线段以开始点的坐标和结束点的坐标表示
	* x1, y1: 开始点的坐标
	* x2, y2: 结束点的坐标
	* color: 颜色
	* 作者: 汤子豪，张径舟
	*/
	int xmin = min(x1, x2), ymin = min(y1, y2), xmax = max(x1, x2), ymax = max(y1, y2);
	if (xmin > winWidth || ymin > winHeight || xmax < 0 || ymax < 0)
		return;
	if (x1 > x2)	// 如果是从右向左绘制的，则交换两点坐标
	{
		swap(x1, x2);
		swap(y1, y2);
	}
	float k, x, y;
	y = y2 - y1;
	x = x2 - x1;
	k = y / x;		// 除零错误编译器帮我搞定了，直接是-inf
	if (0.0 <= k && k <= 1.0)
	{
		__paintline(x, y, x1, y1, color, k);
	}
	else if (k > 1.0)
	{
		__paintline(y, x, x1, y1, color, k);	// k>1的时候需要互换x,y
	}
	else if (-1.0 <= k && k < 0.0)
	{
		__paintline(x, -y, x1, y1, color, k);
	}
	else if (k < -1.0)
	{
		__paintline(-y, x, x1, y1, color, k);
	}
}
bool LinePainter::Paint(Matrix &trans)
{
	/*
	* 绘制线
	* 作者：汤子豪，张径舟
	*/
	Point2D  *temp_points;
	temp_points = WindowCoordinateTransform(this->points, 2, trans);
	//int _new_x1, _new_y1, _new_x2, _new_y2;
	//Matrix p = Matrix::CreatePointMatrix(x1, y1, x2, y2) * trans;
	//_new_x1 = p.data[0][0]; _new_y1 = p.data[0][1]; _new_x2 = p.data[1][0]; _new_y2 = p.data[1][1];
	//double xmin, ymin;
	//xmin = min(temp_points[0].x, temp_points[1].x)
	if (status.GetDisplayType() != dtPixelMode)
		PixelToGrid(temp_points, temp_points, 2);
	__paint(temp_points[0].x, temp_points[0].y, temp_points[1].x, temp_points[1].y, edge_color);
	delete[] temp_points;
	return true;
}
bool LinePainter::Paint()
{
	return this->Paint(this->transformation);
}
Point2D *LinePainter::GetPoints()
{
	Point2D *temp_points = new Point2D[2];
	for (int i = 0; i < 2; i++)
		temp_points[i] = points[i];
	return temp_points;
}
//void LinePainter::ClearData()
//{
//	points
//}



bool PolylinePainter::Paint(Matrix &trans)
{
	//Matrix p = Matrix::CreatePointMatrix(points, point_count) * trans;
	Point2D *temp_points = WindowCoordinateTransform(points, point_count, trans);
	//Point2D *temp_points2 = AdjustToDisplayType(status.GetDisplayType(), temp_points, this->point_count + 1);
	if (status.GetDisplayType() != dtPixelMode)
		PixelToGrid(temp_points, temp_points, point_count);
	for (int i = 0; i < point_count - 1; i++)
	{
		__paint(temp_points[i].x, temp_points[i].y, temp_points[i + 1].x, temp_points[i + 1].y, edge_color);	// 直接调用父类Line的__paint方法
	}
	delete[] temp_points;
	//delete[] temp_points2;
	return true;
}
bool PolylinePainter::Paint()
{
	return this->Paint(this->transformation);
}
Point2D *PolylinePainter::GetPoints()
{
	Point2D *temp_points = new Point2D[point_count];
	for (int i = 0; i < point_count; i++)
		temp_points[i] = points[i];
	return temp_points;
}



bool CirclePainter::Paint()
{
	//if (status.GetDisplayType() != dtPixelMode)
	//{
	//	Point2D temp_points[2];
	//	PixelToGrid(points, temp_points, 2);
	//	__paintedge(temp_points[0].x, temp_points[0].y, temp_points[1].x, temp_points[1].y, edge_color);
	//}
	//else
	//{
	//	__paintedge(points[0].x, points[0].y, points[1].x, points[1].y, edge_color);
	//}
	Paint(this->transformation);
	return true;
}
bool CirclePainter::Paint(Matrix &trans)
{
	Point2D *temp_points = WindowCoordinateTransform(points, 2, trans);
	if (status.GetDisplayType() != dtPixelMode)
		PixelToGrid(temp_points, temp_points, 2);
	__paintedge(temp_points[0].x, temp_points[0].y, temp_points[1].x, temp_points[1].y, edge_color);
	//filler->__paintfill_AET(temp_points, point_count, this->fill_color);
	delete[] temp_points;
	return true;
}
void CirclePainter::__paintedge(int x1, int y1, int x2, int y2, Color color)
{
	/*
	* 绘制圆
	* params:
	*	x1, y1, x2, y2: 圆的直径的两端坐标
	*	color: 颜色
	* 作者: 汤子豪
	*/
	int ox = (x1 + x2) / 2, oy = (y1 + y2) / 2;
	int radius = round(GetLength(x1, y1, x2, y2)) / 2;
	int d = 1 - radius;
	int x = 0, y = radius;
	set_unit(ox, oy + radius, color);
	int dr = radius / 1.414;
	while (x != dr)
	{
		if (d < -0.25)
		{
			d = d + 2 * x + 3;
		}
		else
		{
			y = y - 1;
			d = d + 2 * (x - y) + 5;
		}
		drawthu(ox, oy, x, y, color);
		drawthu(ox, oy, y, x, color);
		++x;
	}
	drawthu(ox, oy, dr, dr, color);
}
void CirclePainter::drawthu(int ox, int oy, int x, int y, Color color)
{
	set_unit(ox + x, oy + y, color);
	set_unit(ox + x, oy - y, color);
	set_unit(ox - x, oy + y, color);
	set_unit(ox - x, oy - y, color);
}



bool EllipsePainter::Paint()
{
	if (status.GetDisplayType() != dtPixelMode)
	{
		Point2D temp_points[2];
		PixelToGrid(points, temp_points, 2);
		__paintedge(temp_points[0].x, temp_points[0].y, temp_points[1].x, temp_points[1].y, edge_color);
	}
	else
	{
		__paintedge(points[0].x, points[0].y, points[1].x, points[1].y, edge_color);
	}
	return true;
}
void EllipsePainter::__paintedge(int x1, int y1, int x2, int y2, Color color)
{
	double a, b;
	a = abs((x1 - x2) / 2);
	b = abs((y1 - y2) / 2);
	int x = 0, y = b;
	double D;
	drawthu(x, y, ((x1 + x2) / 2), ((y1 + y2) / 2), color);
	D = b * b + a * a * (-b + 0.25);
	while (b * b * (x + 1) < a * a * (y - 0.5)) {

		if (D < 0) {
			x = x + 1;
			D = D + b * b * (2 * x + 3);
			drawthu(x, y, ((x1 + x2) / 2), ((y1 + y2) / 2), color);

		}
		else {
			x = x + 1;
			y = y - 1;
			D = D + b * b * (2 * x + 3) + a * a * (2 - 2 * y);
			drawthu(x, y, ((x1 + x2) / 2), ((y1 + y2) / 2), color);
		}
	}
	D = b * b * (x + 0.5) * (x + 0.5) + a * a * (y - 1) * (y - 1) - a * a * b * b;
	while (x <= a && y >= 0)
	{
		if (D < 0)
		{
			x = x + 1;
			y = y - 1;
			D = D + 2 * b * b * (x + 1) + a * a * (3 - y);
			drawthu(x, y, ((x1 + x2) / 2), ((y1 + y2) / 2), color);
		}
		else {
			y = y - 1;
			D = D + a * a * (3 - 2 * y);
			drawthu(x, y, ((x1 + x2) / 2), ((y1 + y2) / 2), color);
		}
	}
}
void EllipsePainter::drawthu(int x, int y, int a, int b, Color color)
{
	set_unit(x + a, y + b, color);
	set_unit(x + a, -y + b, color);
	set_unit(-x + a, y + b, color);
	set_unit(-x + a, -y + b, color);
}



bool PolygonPainter::Paint(Matrix &trans)
{
	//if (status.GetDisplayType() != dtPixelMode)
	//{
	//	Point2D *grid_points = new Point2D[point_count];
	//	PixelToGrid(points, grid_points, point_count);
	//	__paintedge(grid_points, point_count, edge_color);
	//}
	//else
	Point2D *temp_points = WindowCoordinateTransform(points, point_count, trans);
	if (status.GetDisplayType() != dtPixelMode)
		PixelToGrid(temp_points, temp_points, point_count);
	//Matrix p = Matrix::CreatePointMatrix(points, point_count+1) * trans;

	Bucket bucket = Bucket(temp_points, point_count - 1);
	if (bucket.ymin > winHeight + 1 || bucket.ymax < -1 || bucket.scanline_count < 5)
		return true;
	bucket.Correct();
	bucket.Sort();
	bucket.Draw(fill_color);

	__paintedge(temp_points, point_count, edge_color);
	//filler->__paintfill_AET(temp_points, point_count, this->fill_color);
	delete[] temp_points;
	return true;
}
bool PolygonPainter::Paint()
{
	//if (status.GetDisplayType() != dtPixelMode)
	//{
	//	Point2D *grid_points = new Point2D[point_count];
	//	PixelToGrid(points, grid_points, point_count);
	//	__paintedge(grid_points, point_count, edge_color);
	//}
	//else
	this->Paint(this->transformation);
	return true;
}
void PolygonPainter::__paintedge(Point2D *points, int point_count, Color edge_color)
{
	for (int i = 0; i < point_count - 1; i++)
	{
		__paint(points[i].x, points[i].y, points[i + 1].x, points[i + 1].y, edge_color);
	}
}
Point2D *PolygonPainter::GetPoints()
{
	Point2D *temp_points = new Point2D[point_count];
	for (int i = 0; i < point_count; i++)
		temp_points[i] = points[i];
	return temp_points;
}


bool RectanglePainter::Paint()
{
	this->Paint(this->transformation);
	return true;
}
bool RectanglePainter::Paint(Matrix &trans)
{
	Point2D *temp_points = WindowCoordinateTransform(points, point_count, trans);
	if (status.GetDisplayType() != dtPixelMode)
		PixelToGrid(temp_points, temp_points, point_count);
	__paintedge(temp_points, point_count, edge_color);
	delete[] temp_points;
	return true;
}
Point2D *RectanglePainter::GetPoints()
{
	return PolygonPainter::GetPoints();
	//Point2D *temp_points = new Point2D[point_count + 1];
	//for (int i = 0; i < point_count + 1; i++)
	//	temp_points[i] = points[i];
	//return temp_points;
}


//Point3D TrianglePainter::operator[](int index)
//{
//	switch (index)
//	{
//	case 0:
//		return spvertex;
//	case 1:
//		return flatline[0];
//	case 2:
//		return flatline[1];
//	default:
//		return Point3D{-1, -1, -1};
//	}
//}
//TrianglePainter::TrianglePainter(Tin *Tin_grid, int pt1_index, int pt2_index, int pt3_index)
//{
//	this->point_index[0] = pt1_index;
//	this->point_index[1] = pt2_index;
//	this->point_index[2] = pt3_index;
//	this->Tin_grid = Tin_grid;
//	DetermineTriShape();
//}
//
//TrianglePainter::TrianglePainter(Tin *Tin_grid, int pts_index[3])
//{
//	this->point_index[0] = pts_index[0];
//	this->point_index[1] = pts_index[1];
//	this->point_index[2] = pts_index[2];
//	this->Tin_grid = Tin_grid;
//	DetermineTriShape();
//}

//TrianglePainter::TrianglePainter(Point3D p1, Point3D p2, Point3D p3, Model *model)
//{
//	DetermineTriShape(p1, p2, p3);
//	this->model = model;
//}

TrianglePainter::TrianglePainter(Point3D p1, Point3D p2, Point3D p3)
{
	DetermineTriShape(p1, p2, p3);
}

TrianglePainter::TrianglePainter(Point3D spvertex, Point3D flatline1, Point3D flatline2, int trishape)
{
	this->spvertex = spvertex;
	this->flatline[0] = flatline1;
	this->flatline[1] = flatline2;
	this->trishape = trishape;
}

//TrianglePainter::TrianglePainter(Face &other)
//{
//	spvertex = other[other.spvertex];
//	flatline[0] = other[other.flatline[0]];
//	flatline[1] = other[other.flatline[1]];
//	this->trishape = other.trishape;
//
//}

/// <summary>
/// 判断光栅化三角形形状以及顶点状态
/// </summary>
//void TrianglePainter::DetermineTriShape()
//{
//	if ((*Tin_grid)[this->point_index[0]].y == (*Tin_grid)[this->point_index[1]].y)
//	{
//		this->spvertex = this->point_index[2];
//		this->flatline[0] = this->point_index[0];
//		this->flatline[1] = this->point_index[1];
//		if ((*Tin_grid)[this->point_index[2]].y > (*Tin_grid)[this->point_index[1]].y)
//		{
//			this->trishape = 1;//平底
//		}
//		else this->trishape = 2;//平顶
//	}
//	else if ((*Tin_grid)[this->point_index[1]].y == (*Tin_grid)[this->point_index[2]].y)
//	{
//		this->spvertex = this->point_index[0];
//		this->flatline[0] = this->point_index[1];
//		this->flatline[1] = this->point_index[2];
//		if ((*Tin_grid)[this->point_index[0]].y > (*Tin_grid)[this->point_index[1]].y)
//		{
//			this->trishape = 1;//平底
//		}
//		else this->trishape = 2;//平顶
//	}
//	else if ((*Tin_grid)[this->point_index[0]].y == (*Tin_grid)[this->point_index[2]].y)
//	{
//		this->spvertex = this->point_index[1];
//		this->flatline[0] = this->point_index[2];
//		this->flatline[1] = this->point_index[0];
//		if ((*Tin_grid)[this->point_index[1]].y > (*Tin_grid)[this->point_index[1]].y)
//		{
//			this->trishape = 1;//平底
//		}
//		else this->trishape = 2;//平顶
//	}
//	else//不平底或平顶
//	{
//		this->trishape = 0;
//		if ((*Tin_grid)[this->point_index[0]].y < (*Tin_grid)[this->point_index[1]].y)
//		{
//			if ((*Tin_grid)[this->point_index[0]].y > (*Tin_grid)[this->point_index[2]].y)//2在最下面
//			{
//				this->spvertex = this->point_index[0];
//			}
//			else if ((*Tin_grid)[this->point_index[2]].y > (*Tin_grid)[this->point_index[1]].y)//2在最上面
//			{
//				this->spvertex = this->point_index[1];
//			}
//			else this->spvertex = this->point_index[2];//2在中间
//		}
//		else
//		{
//			if ((*Tin_grid)[this->point_index[1]].y > (*Tin_grid)[this->point_index[2]].y)//2在最下面
//			{
//				this->spvertex = this->point_index[1];
//			}
//			else if ((*Tin_grid)[this->point_index[2]].y > (*Tin_grid)[this->point_index[0]].y)//2在最上面
//			{
//				this->spvertex = this->point_index[0];
//			}
//			else this->spvertex = this->point_index[2];//2在中间
//		}
//	}
//}
void TrianglePainter::PaintflatTri()
{
	int flag = 1;				// 用来在后面判断扫描线从下往上还是从上往下
	double bottomy, topy;
	bottomy = this->flatline[0].y;
	topy = this->spvertex.y;
	if (bottomy > topy)
	{
		flag = 0;
		swap(bottomy, topy);
	}

	double reciprocal_k1 = (this->spvertex.x - this->flatline[0].x) / (this->spvertex.y - this->flatline[0].y);
	double reciprocal_k2 = (this->spvertex.x - this->flatline[1].x) / (this->spvertex.y - this->flatline[1].y);
	double reciprocal_kz1 = (this->spvertex.z - this->flatline[0].z) / (this->spvertex.y - this->flatline[0].y);
	double reciprocal_kz2 = (this->spvertex.z - this->flatline[1].z) / (this->spvertex.y - this->flatline[1].y);
	double x1 = this->flatline[0].x;
	double x2 = this->flatline[1].x;
	double left_z = flatline[0].z, right_z = flatline[1].z;

	// 这个函数里带有“I”的都是光强计算相关的
	double left_IR, right_IR;
	double left_IG, right_IG;
	double left_IB, right_IB;
	double left_rk_IR, right_rk_IR;
	double left_rk_IG, right_rk_IG;
	double left_rk_IB, right_rk_IB;
	if (model != nullptr)
	{
		left_IR = model->vertex_IR[flatline[0].id]; right_IR = model->vertex_IR[flatline[1].id];
		left_IG = model->vertex_IG[flatline[0].id]; right_IG = model->vertex_IG[flatline[1].id];
		left_IB = model->vertex_IB[flatline[0].id]; right_IB = model->vertex_IB[flatline[1].id];

		left_rk_IR = (model->vertex_IR[spvertex.id] - model->vertex_IR[flatline[0].id]) / (spvertex.y - flatline[0].y);
		left_rk_IG = (model->vertex_IG[spvertex.id] - model->vertex_IG[flatline[0].id]) / (spvertex.y - flatline[0].y);
		left_rk_IB = (model->vertex_IB[spvertex.id] - model->vertex_IB[flatline[0].id]) / (spvertex.y - flatline[0].y);

		right_rk_IR = (model->vertex_IR[spvertex.id] - model->vertex_IR[flatline[1].id]) / (spvertex.y - flatline[1].y);
		right_rk_IG = (model->vertex_IG[spvertex.id] - model->vertex_IG[flatline[1].id]) / (spvertex.y - flatline[1].y);
		right_rk_IB = (model->vertex_IB[spvertex.id] - model->vertex_IB[flatline[1].id]) / (spvertex.y - flatline[1].y);
	}


	if (x1 > x2)//如果x1大于x2，将两者互换
	{
		swap(x1, x2);
		swap(reciprocal_k1, reciprocal_k2);
		swap(left_z, right_z);
		swap(reciprocal_kz1, reciprocal_kz2);

		if (model != nullptr)
		{
			swap(left_IR, right_IR);
			swap(left_rk_IR, right_rk_IR);
			swap(left_IG, right_IG);
			swap(left_rk_IG, right_rk_IG);
			swap(left_IB, right_IB);
			swap(left_rk_IB, right_rk_IB);
		}
	}


	if (flag)
	{
		for (double y = bottomy; y <= topy; ++y)	//从下往上
		{
			double reciprocal_hori_kz = (right_z - left_z) / (x2 - x1);
			double temp_left_z = left_z;

			double temp_left_IR, temp_left_IG, temp_left_IB;
			double r_hori_k_IR, r_hori_k_IG, r_hori_k_IB;
			if (model != nullptr)
			{
				temp_left_IR = left_IR;
				temp_left_IG = left_IG;
				temp_left_IB = left_IB;
				r_hori_k_IR = (right_IR - left_IR) / (x2 - x1);
				r_hori_k_IG = (right_IG - left_IG) / (x2 - x1);
				r_hori_k_IB = (right_IB - left_IB) / (x2 - x1);
			}

			for (int x = x1-1; x <= x2; ++x)
			{
				Color color;
				if (model != nullptr)
				{
					int R = temp_left_IR > 1.0 ? 255 : temp_left_IR * 255;
					int G = temp_left_IG > 1.0 ? 255 : temp_left_IG * 255;
					int B = temp_left_IB > 1.0 ? 255 : temp_left_IB * 255;
					color = _RGB(R, G, B);
				}
				else
					color = fill_color;
				setPixel(x, y, temp_left_z, color);
				//setPixel(x, y, _RGB(0, 0, 0));
				temp_left_z += reciprocal_hori_kz;
				if (model != nullptr)
				{
					temp_left_IR += r_hori_k_IR;
					temp_left_IG += r_hori_k_IG;
					temp_left_IB += r_hori_k_IB;
				}
			}
			x1 += reciprocal_k1;
			x2 += reciprocal_k2;
			left_z += reciprocal_kz1;
			right_z += reciprocal_kz2;
			if (model != nullptr)
			{
				left_IR += left_rk_IR;
				left_IG += left_rk_IG;
				left_IB += left_rk_IB;
				right_IR += right_rk_IR;
				right_IG += right_rk_IG;
				right_IB += right_rk_IB;
			}
		}
	}
	else
	{
		for (double y = topy; y >= bottomy; --y)	//从上往下
		{
			double reciprocal_hori_kz = (right_z - left_z) / (x2 - x1);
			double temp_left_z = left_z;

			double temp_left_IR, temp_left_IG, temp_left_IB;
			double r_hori_k_IR, r_hori_k_IG, r_hori_k_IB;
			if (model != nullptr)
			{
				temp_left_IR = left_IR;
				temp_left_IG = left_IG;
				temp_left_IB = left_IB;
				r_hori_k_IR = (right_IR - left_IR) / (x2 - x1);
				r_hori_k_IG = (right_IG - left_IG) / (x2 - x1);
				r_hori_k_IB = (right_IB - left_IB) / (x2 - x1);
			}

			for (int x = x1-1; x <= x2; ++x)
			{
				Color color;
				if (model != nullptr)
				{
					int R = temp_left_IR > 1.0 ? 255 : temp_left_IR * 255;
					int G = temp_left_IG > 1.0 ? 255 : temp_left_IG * 255;
					int B = temp_left_IB > 1.0 ? 255 : temp_left_IB * 255;
					color = _RGB(R, G, B);
				}
				else
					color = fill_color;
				setPixel(x, y, temp_left_z, color);
				//setPixel(x, y, _RGB(0, 0, 0));
				temp_left_z += reciprocal_hori_kz;
				if (model != nullptr)
				{
					temp_left_IR += r_hori_k_IR;
					temp_left_IG += r_hori_k_IG;
					temp_left_IB += r_hori_k_IB;
				}
			}
			x1 -= reciprocal_k1;
			x2 -= reciprocal_k2;
			left_z -= reciprocal_kz1;
			right_z -= reciprocal_kz2;
			if (model != nullptr)
			{
				left_IR -= left_rk_IR;
				left_IG -= left_rk_IG;
				left_IB -= left_rk_IB;
				right_IR -= right_rk_IR;
				right_IG -= right_rk_IG;
				right_IB -= right_rk_IB;
			}
		}
	}
}

// 绘制平顶，已弃用
void TrianglePainter::PaintflattopTri()
{
	/*
		*********
		 *     *
		  *   *
		   * *
			* 
	*/
	double bottomy = spvertex.y;//最低y
	double topy = flatline[0].y;//最高y
	double left_z = flatline[0].z, right_z = flatline[1].z;
	double left_IR, right_IR;
	double left_IG, right_IG;
	double left_IB, right_IB;
	double left_rk_IR, right_rk_IR;
	double left_rk_IG, right_rk_IG;
	double left_rk_IB, right_rk_IB;

	if (model != nullptr)
	{
		left_IR = model->vertex_IR[flatline[0].id]; right_IR = model->vertex_IR[flatline[1].id];
		left_IG = model->vertex_IG[flatline[0].id]; right_IG = model->vertex_IG[flatline[1].id];
		left_IB = model->vertex_IB[flatline[0].id]; right_IB = model->vertex_IB[flatline[1].id];

		left_rk_IR = (model->vertex_IR[spvertex.id] - model->vertex_IR[flatline[0].id]) / (spvertex.y - flatline[0].y);
		left_rk_IG = (model->vertex_IG[spvertex.id] - model->vertex_IG[flatline[0].id]) / (spvertex.y - flatline[0].y);
		left_rk_IB = (model->vertex_IB[spvertex.id] - model->vertex_IB[flatline[0].id]) / (spvertex.y - flatline[0].y);

		right_rk_IR = (model->vertex_IR[spvertex.id] - model->vertex_IR[flatline[1].id]) / (spvertex.y - flatline[1].y);
		right_rk_IG = (model->vertex_IG[spvertex.id] - model->vertex_IG[flatline[1].id]) / (spvertex.y - flatline[1].y);
		right_rk_IB = (model->vertex_IB[spvertex.id] - model->vertex_IB[flatline[1].id]) / (spvertex.y - flatline[1].y);
	}

	//k倒数
	double reciprocal_k1 = (this->spvertex.x - this->flatline[0].x) / (this->spvertex.y - this->flatline[0].y);
	double reciprocal_k2 = (this->spvertex.x - this->flatline[1].x) / (this->spvertex.y - this->flatline[1].y);
	double reciprocal_kz1 = (this->spvertex.z - this->flatline[0].z) / (this->spvertex.y - this->flatline[0].y);
	double reciprocal_kz2 = (this->spvertex.z - this->flatline[1].z) / (this->spvertex.y - this->flatline[1].y);

	double x1 = this->flatline[0].x;
	double x2 = this->flatline[1].x;

	if (x1 > x2)//如果x1大于x2，将两者互换
	{
		swap(x1, x2);
		swap(reciprocal_k1, reciprocal_k2);
		swap(left_z, right_z);
		swap(reciprocal_kz1, reciprocal_kz2);
		if (model != nullptr)
		{
			swap(left_IR, right_IR);
			swap(left_rk_IR, right_rk_IR);
			swap(left_IG, right_IG);
			swap(left_rk_IG, right_rk_IG);
			swap(left_IB, right_IB);
			swap(left_rk_IB, right_rk_IB);
		}
	}

	int i = 0;
	for (double y = topy; y >= bottomy; --y)//循环画水平扫描线
	{
		double reciprocal_hori_kz = (right_z - left_z) / (x2 - x1);
		double temp_left_z = left_z;

		double temp_left_IR, temp_left_IG, temp_left_IB;
		if (model != nullptr)
		{
			temp_left_IR = left_IR;
			temp_left_IG = left_IG;
			temp_left_IB = left_IB;
		}


		double r_hori_k_IR, r_hori_k_IG, r_hori_k_IB;
		if (model != nullptr)
		{
			r_hori_k_IR = (right_IR - left_IR) / (x2 - x1);
			r_hori_k_IG = (right_IG - left_IG) / (x2 - x1);
			r_hori_k_IB = (right_IB - left_IB) / (x2 - x1);
		}
		for (int x = x1; x <= x2; ++x)
		{
			Color color;
			if (model != nullptr)
				color = _RGB(int(temp_left_IR * 255), int(temp_left_IG * 255), int(temp_left_IB * 255));
			else
				color = this->fill_color;
			setPixel(x, y, temp_left_z, color);
			//setPixel(x, y, _RGB(0, 0, 0));
			temp_left_z += reciprocal_hori_kz;
			if (model != nullptr)
			{
				temp_left_IR += r_hori_k_IR;
				temp_left_IG += r_hori_k_IG;
				temp_left_IB += r_hori_k_IB;
			}

		}
		//set_unit(x1, y, _RGB(0, 0, 0));
		//set_unit(x2, y, _RGB(0, 0, 0));
		x1 -= reciprocal_k1;
		x2 -= reciprocal_k2;
		left_z -= reciprocal_kz1;
		right_z -= reciprocal_kz2;
		if (model != nullptr)
		{
			left_IR -= left_rk_IR;
			left_IG -= left_rk_IG;
			left_IB -= left_rk_IB;
			right_IR -= right_rk_IR;
			right_IG -= right_rk_IG;
			right_IB -= right_rk_IB;
		}

	}
}
// 绘制平底，已弃用
void TrianglePainter::PaintflatbottomTri()
{
	/*
			*
		   * *
	      *   *
	     *     *
		*********
	*/
	double bottomy = this->flatline[0].y;
	double topy = this->spvertex.y;
	double left_z = flatline[0].z, right_z = flatline[1].z;

	//k倒数
	double reciprocal_k1 = (double)(this->spvertex.x - this->flatline[0].x) / (double)(this->spvertex.y - this->flatline[0].y);
	double reciprocal_k2 = (double)(this->spvertex.x - this->flatline[1].x) / (double)(this->spvertex.y - this->flatline[1].y);
	double reciprocal_kz1 = (this->spvertex.z - this->flatline[0].z) / (this->spvertex.y - this->flatline[0].y);
	double reciprocal_kz2 = (this->spvertex.z - this->flatline[1].z) / (this->spvertex.y - this->flatline[1].y);
	double x1 = this->flatline[0].x;
	double x2 = this->flatline[1].x;

	if (x1 > x2)//如果x1大于x2，将两者互换
	{
		swap(x1, x2);
		swap(reciprocal_k1, reciprocal_k2);
		swap(left_z, right_z);
		swap(reciprocal_kz1, reciprocal_kz2);
	}

	int i = 0;
	for (double y = bottomy; y <= topy; ++y)//循环画水平扫描线
	{
		double reciprocal_hori_kz = (right_z - left_z) / (x2 - x1);
		double temp_left_z = left_z;
		for (int x = x1; x <= x2; ++x)
		{
			setPixel(x, y, temp_left_z, this->fill_color);
			//setPixel(x, y, _RGB(0, 0, 0));
			temp_left_z += reciprocal_hori_kz;
		}
		//for (int x = x1; x <= x2; ++x)
		//{
		//	setPixel(x, y, _RGB(0, 0, 0));
		//}
		//set_unit(x1, y, _RGB(0, 0, 0));
		//set_unit(x2, y, _RGB(0, 0, 0));
		/*for (int x = x1; x <= x2; ++x)
		{
			set_unit(x, y, _RGB(0, 0, 255));
		}*/
		x1 += reciprocal_k1;
		x2 += reciprocal_k2;
		left_z += reciprocal_kz1;
		right_z += reciprocal_kz2;
	}
}
// 绘制既不是平顶也不是平底
void TrianglePainter::PaintJointTri()
{
	double midy = this->spvertex.y;//中间顶点的y

	//被分割边的k
	double reciprocal_k = (double)(this->flatline[1].x - this->flatline[0].x) / (double)(this->flatline[1].y - this->flatline[0].y);

	//一个顶点
	double y1 = this->flatline[0].y;
	double x1 = this->flatline[0].x;

	double crossx = (midy - y1) * reciprocal_k + x1;//切割交点x
	double mid_depth = (midy - this->flatline[0].y) / (this->flatline[1].y - this->flatline[0].y) * (this->flatline[1].z - this->flatline[0].z) + this->flatline[0].z;

	Point3D crosspt(crossx, midy, mid_depth);
	double mid_IR, mid_IG, mid_IB;
	if (model != nullptr)
	{
		crosspt.id = model->vertex_IR.size();
		mid_IR = (midy - this->flatline[0].y) / (this->flatline[1].y - this->flatline[0].y) * (model->vertex_IR[flatline[1].id] - model->vertex_IR[flatline[0].id]) + model->vertex_IR[flatline[0].id];
		mid_IG = (midy - this->flatline[0].y) / (this->flatline[1].y - this->flatline[0].y) * (model->vertex_IG[flatline[1].id] - model->vertex_IG[flatline[0].id]) + model->vertex_IG[flatline[0].id];
		mid_IB = (midy - this->flatline[0].y) / (this->flatline[1].y - this->flatline[0].y) * (model->vertex_IB[flatline[1].id] - model->vertex_IB[flatline[0].id]) + model->vertex_IB[flatline[0].id];
		model->vertex_IR.push_back(mid_IR);
		model->vertex_IG.push_back(mid_IG);
		model->vertex_IB.push_back(mid_IB);
	}

	//TrianglePainter tri1(crosspt, this->spvertex, this->flatline[0]);
	//tri1.spvertex = this->flatline[0];
	//tri1.flatline[0] = newindex; tri1.flatline[1] = this->spvertex;

	//TrianglePainter tri2(Tin_grid, newindex, this->spvertex, this->flatline[1]);
	//tri2.spvertex = this->flatline[1];
	//tri2.flatline[0] = newindex; tri1.flatline[1] = this->spvertex;

	if (this->flatline[0].y > this->flatline[1].y)
	{
		TrianglePainter tri1 = TrianglePainter(this->flatline[0], crosspt, this->spvertex, 1);
		tri1.fill_color = this->fill_color;
		tri1.model = this->model;
		//tri1.PaintflatbottomTri();
		tri1.PaintflatTri();

		TrianglePainter tri2 = TrianglePainter(this->flatline[1], crosspt, this->spvertex, 2);
		tri2.fill_color = this->fill_color;
		tri2.model = this->model;
		//tri2.PaintflattopTri();
		tri2.PaintflatTri();
	}
	else
	{
		TrianglePainter tri1 = TrianglePainter(this->flatline[1], crosspt, this->spvertex, 1);
		tri1.fill_color = this->fill_color;
		tri1.model = this->model;
		//tri1.PaintflatbottomTri();
		tri1.PaintflatTri();

		TrianglePainter tri2 = TrianglePainter(this->flatline[0], crosspt, this->spvertex, 2);
		tri2.fill_color = this->fill_color;
		tri2.model = this->model;
		//tri2.PaintflattopTri();
		tri2.PaintflatTri();
	}
	if (model != nullptr)
	{
		model->vertex_IR.pop_back();
		model->vertex_IG.pop_back();
		model->vertex_IB.pop_back();
	}
}

bool TrianglePainter::Paint()
{
	//if (this->trishape == 1)//平底
	//{
	//	this->PaintflatbottomTri();
	//}

	//else if (this->trishape == 2)//平顶
	//{
	//	this->PaintflattopTri();
	//}
	if (this->trishape == 1 || this->trishape == 2)
		this->PaintflatTri();
	else this->PaintJointTri();//拼接

	return true;
}

/// <summary>
/// 判断光栅化三角形形状以及顶点状态
/// </summary>
void TrianglePainter::DetermineTriShape(Point3D p1, Point3D p2, Point3D p3)
{
	if (abs(p1.y - p2.y) < 1e-6)
	{
		this->spvertex = p3;
		this->flatline[0] = p1;
		this->flatline[1] = p2;
		if (p3.y > p2.y)
		{
			this->trishape = 1;//平底
		}
		else this->trishape = 2;//平顶
	}
	else if (abs(p2.y - p3.y) < 1e-6)
	{
		this->spvertex = p1;
		this->flatline[0] = p2;
		this->flatline[1] = p3;
		if (p1.y > p2.y)
		{
			this->trishape = 1;//平底
		}
		else this->trishape = 2;//平顶
	}
	else if (abs(p1.y - p3.y) < 1e-6)
	{
		this->spvertex = p2;
		this->flatline[0] = p3;
		this->flatline[1] = p1;
		if (p2.y > p1.y)
		{
			this->trishape = 1;//平底
		}
		else this->trishape = 2;//平顶
	}
	else//不平底或平顶
	{
		this->trishape = 0;
		if (p1.y < p2.y)
		{
			if (p1.y > p3.y)//2在最下面
			{
				this->spvertex = p1;
				this->flatline[0] = p2;
				this->flatline[1] = p3;
			}
			else if (p3.y > p2.y)//2在最上面
			{
				this->spvertex = p2;
				this->flatline[0] = p1;
				this->flatline[1] = p3;
			}
			else
			{
				this->spvertex = p3;//2在中间
				this->flatline[0] = p1;
				this->flatline[1] = p2;
			}
				
		}
		else
		{
			if (p2.y > p3.y)//2在最下面
			{
				this->spvertex = p2;
				this->flatline[0] = p1;
				this->flatline[1] = p3;
			}
			else if (p3.y > p1.y)//2在最上面
			{
				this->spvertex = p1;
				this->flatline[0] = p2;
				this->flatline[1] = p3;
			}
			else
			{
				this->spvertex = p3;//2在中间
				this->flatline[0] = p1;
				this->flatline[1] = p2;
			}
		}
	}
}

bool Line3DPainter::Paint()
{
	double rk = (pt2.x - pt1.x) / (pt2.y - pt1.y);
	double rkz = (pt2.z - pt1.z) / (pt2.y - pt1.y);
	int temp_x = (int)pt1.x;
	double temp_z = pt1.z;
	for (int temp_y = (int)pt1.y; temp_y <= pt2.y; temp_y++)
	{
		setPixel(temp_x, temp_y, temp_z, edge_color);
		temp_x += rk;
		temp_z += rkz;
	}
	return true;
}
