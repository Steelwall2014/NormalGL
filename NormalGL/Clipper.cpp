#include "Clipper.h"
#include <windows.h>
#include <math.h>
#include <vector>
#include <assert.h>
#include <locale.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "Painters.h"
//#include "centercut.h"
using namespace std;
/// <summary>
/// 出入点类，继承Point2D，？？？类直接继承了结构体也不知道有没有问题
/// </summary>
class CrossPoint :public Point2D {
public:
	double x, y;
	int pindex[2] = { 0,0 };//对应的多边形边,对应的裁剪多边形的边
	int flag = 0;//flag出点为0，入点为1
	CrossPoint()
	{
		this->x = 0;
		this->y = 0;
		this->pindex[0] = 0;
		this->pindex[1] = 0;
		this->flag = 0;
	}
	CrossPoint(double x, double y, int flag, int pindex, int cindex)
	{
		this->flag = flag;
		this->x = x;
		this->y = y;
		this->pindex[0] = pindex;
		this->pindex[1] = cindex;
	}

	CrossPoint(Point2D pt, int flag, int pindex, int cindex)
	{
		this->x = pt.x;
		this->y = pt.y;
		this->flag = flag;
		this->pindex[0] = pindex;
		this->pindex[1] = cindex;
	}

	CrossPoint(Point2D pt)//用于存储非交点，利用flag为-1区分
	{
		this->x = pt.x;
		this->y = pt.y;
		this->flag = -1;
	}
	bool operator==(CrossPoint &p1)
	{
		if (fabs(this->x - p1.x) < 0.001 && fabs(this->y - p1.y) < 0.001) return true;
		else return false;
	}
	bool operator==(Point2D &p1)
	{
		if (fabs(this->x - p1.x) < 0.001 && fabs(this->y - p1.y) < 0.001) return true;
		else return false;
	}
};
bool centercut(PixelPoint *p, int x1, int x2, int y1, int y2, int threshold);
void ChangeClockwise(Point2D *poly, int num);
int RTInside(Point2D s, Point2D e, Point2D p);
vector<Point2D> subSH_Clip(vector<Point2D> pg_vec, vector<Point2D> clip_vec);
Point2D GetCrossPoint(Point2D p1, Point2D p2, Point2D s1, Point2D s2);
/// <summary>
/// 将在同一条裁剪窗口边上的点进行顺时针排列
/// </summary>
/// <param name="cross">交点</param>
/// <param name="clip">裁剪点</param>
void GetRightOrder(vector<CrossPoint> &cross, vector<Point2D> clip)
{
	for (int i = 0; i < cross.size() - 1; ++i)
	{
		for (int j = i + 1; j < cross.size(); j++)
		{
			if (cross[i].pindex[1] == cross[j].pindex[1])
			{
				int index = cross[i].pindex[1];
				int disa = abs(cross[i].x - clip[index].x);
				int disb = abs(cross[j].x - clip[index].x);
				if (disa > disb)
				{
					CrossPoint temp = cross[i];
					cross[i] = cross[j];
					cross[j] = temp;
				}
				else if (disa == disb)
				{
					int disya = abs(cross[i].y - clip[index].y);
					int disyb = abs(cross[j].y - clip[index].y);
					if (disya > disyb)
					{
						CrossPoint temp = cross[i];
						cross[i] = cross[j];
						cross[j] = temp;
					}
				}
			}
		}
	}
}
bool is_clockwise(Point2D *points, int point_count)
{
	// points需要首尾相连
	double s = 0.0;
	for (int i = 0; i < point_count-1; i++)
	{
		s += (points[i].y + points[i+1].y) * (points[i+1].x - points[i].x);
	}
	return s > 0.0;
}
void to_clockwise(Point2D *points, int point_count)
{
	if (is_clockwise(points, point_count))
		return;
	else
	{
		Point2D temp;
		//Point2D *newpoly = new Point2D[num];
		for (int i = 0; i < point_count / 2; i++)
		{
			temp = points[i];
			points[i] = points[point_count - 1 - i];
			points[point_count - 1 - i] = temp;
			//newpoly[i] = poly[num - 1 - i];
		}
	}
}

void to_anticlockwise(Point2D *points, int point_count)
{
	if (!is_clockwise(points, point_count))
		return;
	else
	{
		Point2D temp;
		for (int i = 0; i < point_count / 2; i++)
		{
			temp = points[i];
			points[i] = points[point_count - 1 - i];
			points[point_count - 1 - i] = temp;
		}
	}
}

int FindXmaxPt(Point2D *polygon, int nodenum);
Point2D *JudgeClockWise(Point2D *polygon, int nodenum, int pos);

int GetIntersectx(int y, Point2D p1, Point2D p2);
int ptInPolyByCorner(Point2D pt, vector<Point2D> polygon);
int ptInPolyByCorner(CrossPoint pt, vector<CrossPoint> polygon);
vector<CrossPoint> returnSortClipCross(vector<CrossPoint> crossP);
vector<CrossPoint> GetIndex(vector<CrossPoint> crossP, vector<Point2D> clip);
vector<CrossPoint> Merge(vector<CrossPoint> cross, vector<Point2D> unmerge, int flag, vector<Point2D> clip);
vector<CrossPoint> GetCrossPts(vector<Point2D> &pg_vec, vector<Point2D> &clip_vec);
vector<vector<CrossPoint>> WA_Clip(vector<Point2D> pg_vec, vector<Point2D> clip_vec);





Clipper::Clipper(vector<Point2D> &boundary)
{
	SetBoundary(boundary);
}
Clipper::Clipper(Point2D window_lower_left, Point2D window_upper_right)
{
	SetBoundary(window_lower_left, window_upper_right);
}
void Clipper::Reset()
{
	boundary.clear();
}
void Clipper::SetBoundary(Point2D window_lower_left, Point2D window_upper_right)
{
	boundary_is_rect = true;
	this->boundary.clear();
	this->boundary.push_back(window_lower_left);
	this->boundary.push_back(Point2D(window_lower_left.x, window_upper_right.y));
	this->boundary.push_back(window_upper_right);
	this->boundary.push_back(Point2D(window_upper_right.x, window_lower_left.y));
	this->boundary.push_back(window_lower_left);
}
void Clipper::SetBoundary(vector<Point2D> &boundary)
{
	this->boundary.clear();
	for (int i = 0; i < boundary.size(); i++)
		this->boundary.push_back(boundary[i]);
	this->boundary.push_back(boundary[0]);
}
bool Clipper::CenterLineClip(Painter *painter, vector<Painter *> &container)
{
	if (!boundary_is_rect)
		throw("The boundary must be a rectangle when using CenterLineClip");
	if (painter->GetPaintShape() != "POLYLINE" && painter->GetPaintShape() != "LINE")
		throw("This function(CenterLineClip) only applies to lines");
	bool flag = false; // 用来标注有没有剪裁成功
	for (int i = 0; i < painter->GetPointCount() - 1; i++)
	{
		PixelPoint p[2];
		p[0] = PixelPoint{ (int)painter->GetX(i), (int)painter->GetY(i) };
		p[1] = PixelPoint{ (int)painter->GetX(i + 1), (int)painter->GetY(i + 1) };
		bool temp_flag = centercut(p, boundary[0].x, boundary[3].x, boundary[0].y, boundary[1].y, 2);	// 如果是false，那么说明这个线段裁剪后已经没了
		flag |= temp_flag;
		if (temp_flag)
		{
			LinePainter *new_painter = new LinePainter(p[0].x, p[0].y, p[1].x, p[1].y);
			container.push_back(new_painter);
		}

		//new_painter->SetX(0, p[0].x); new_painter->SetY(0, p[0].y);
		//new_painter->SetX(1, p[1].x); new_painter->SetY(1, p[1].y);
	}
	return flag;
	//return painters;
	//else
	//{
	//	//new_painter = (PolylinePainter *)new_painter;
	//	vector<Painter *> *painters = new vector<Painter *>;
	//	for (int i = 0; i < painter->GetPointCount() - 1; i++)
	//	{
	//		PixelPoint p[2];
	//		p[0] = PixelPoint{ (int)painter->GetX(i), (int)painter->GetY(i) };
	//		p[1] = PixelPoint{ (int)painter->GetX(i+1), (int)painter->GetY(i+1) };
	//		centercut(p, boundary[0].x, boundary[3].x, boundary[0].y, boundary[1].y, 2);
	//		LinePainter *new_painter = new LinePainter(p[0].x, p[0].y, p[1].x, p[1].y);
	//		painters->push_back(new_painter);
	//		////if (!is_in_threshold(PixelPoint{ (int)painter->GetX(i), (int)painter->GetY(i) }, boundary[0].x, boundary[3].x, boundary[0].y, boundary[1].y, 2))
	//		//	new_points.push_back(p[0]);
	//		////if (!is_in_threshold(PixelPoint{ (int)painter->GetX(i + 1), (int)painter->GetY(i + 1) }, boundary[0].x, boundary[3].x, boundary[0].y, boundary[1].y, 2))
	//		//	new_points.push_back(p[1]);
	//	}
	//	//new_painter->ClearData();
	//	//new_painter->points = new Point2D[new_points.size()];
	//	//for (int i = 0; i < new_points.size(); i++)
	//	//{
	//	//	new_painter->SetX(i, new_points[i].x); 
	//	//	new_painter->SetY(i, new_points[i].y);
	//	//}
	//		
	//		//new_painter->SetX(i + 1, p2.x); new_painter->SetY(1, p2.y);
	//}
}
void Clipper::SutherlandHodgmanClip(Painter *painter, Painter *new_painter)
{
	Point2D *pg = painter->GetPoints(), *clip;
	int pgnum = painter->GetPointCount(), clipnum = this->boundary.size();
	clip = new Point2D[clipnum];
	for (int i = 0; i < clipnum; i++)
		clip[i] = this->boundary[i];
	to_anticlockwise(pg, pgnum);			// 这里实现的sh算法其实是不需要多边形首尾相连的，但是这个函数需要，所以前面分配了多一个空间的内存，这边传入的是点的数量+1
	to_anticlockwise(clip, clipnum);		// 不过我打算这次汇报之后把点的数量统一成带上最后一个点的数量，而不是真实点的数量
	vector<Point2D> pgvec(pg, pg + pgnum - 1);
	vector<Point2D> clipvec(clip, clip + clipnum - 1);

	vector<Point2D> clipped;				// 这个结果是首尾相连的
	clipped = subSH_Clip(pgvec, clipvec);
	delete[] new_painter->points;
	new_painter->points = new Point2D[clipped.size()];
	new_painter->point_count = clipped.size();
	for (int i = 0; i < clipped.size(); i++)
		new_painter->points[i] = clipped[i];
	delete[] pg;
	delete[] clip;
	//return clipped;
}
int Clipper::WeilerAthertonClip(Painter *painter, vector<Painter *> &container)
{
	Point2D *pg = painter->GetPoints(), *clip;
	int pgnum = painter->GetPointCount(), clipnum = this->boundary.size();
	clip = new Point2D[clipnum];
	for (int i = 0; i < clipnum; i++)
		clip[i] = this->boundary[i];
	//pg = JudgeClockWise(pg, pgnum, FindXmaxPt(pg, pgnum));
	//clip = JudgeClockWise(clip, clipnum, FindXmaxPt(clip, clipnum));
	to_clockwise(pg, pgnum);
	to_clockwise(clip, clipnum);
	vector<Point2D> pg_vec(pg, pg + pgnum);
	vector<Point2D> clip_vec(clip, clip + clipnum);
	int i = 0;
	for (; i < pg_vec.size(); i++)
	{
		if (!ptInPolyByCorner(pg_vec[i], clip_vec))
			break;
	}
	if (i == pg_vec.size())
	{
		return 0;			// 如果被剪裁多边形完全在剪裁多边形内就直接返回，用0代表这种情况
	}
	for (i = 0; i < clip_vec.size(); i++)
	{
		if (!ptInPolyByCorner(clip_vec[i], pg_vec))
			break;
	}
	if (i == clip_vec.size())
	{
		PolygonPainter *new_painter = new PolygonPainter(clip_vec.data(), clip_vec.size());
		new_painter->fill_color = painter->fill_color;
		container.push_back(new_painter);
		return 1;			// 如果剪裁多边形完全在被剪裁多边形内就直接返回，用1代表这种情况
	}


	vector<vector<CrossPoint>> clippedvector = WA_Clip(pg_vec, clip_vec);	// 结果是首尾相连的
	for (auto &clip_result : clippedvector)
	{
		Point2D *points = new Point2D[clip_result.size()];
		for (int i = 0; i < clip_result.size(); i++)
			points[i] = Point2D{ clip_result[i].x, clip_result[i].y };
		PolygonPainter *new_painter = new PolygonPainter(points, clip_result.size());
		new_painter->fill_color = painter->fill_color;
		container.push_back(new_painter);
		delete[] points;
	}
	delete[] pg;
	delete[] clip;
	return 2;
}



/*点所在区域编号
点：p
x1：左
x2：右
y1：下
y2：上*/
int encode(PixelPoint p, int  x1, int x2, int y1, int y2)
{
	int a = 0;
	int xmin, xmax, ymin, ymax;
	xmin = min(x1, x2);
	xmax = max(x1, x2);
	ymin = min(y1, y2);
	ymax = max(y1, y2);
	if (p.x < xmin)
		a |= 1;	//0001
	else if (p.x > xmax)
		a |= 2;	//0010
	if (p.y < ymin)
		a |= 4; //0100
	else if (p.y > ymax)
		a |= 8; //1000

	//if (p.x >= x1)
	//{
	//	if (p.x <= x2)
	//	{
	//		if (p.y >= y1)
	//		{
	//			if (p.y <= y2)
	//			{
	//				a = 0;
	//			}
	//			else { a = 1; }
	//		}
	//		else { a = 2; }
	//	}
	//	else {
	//		if (p.y >= y1)
	//		{
	//			if (p.y <= y2)
	//			{
	//				a = 3;
	//			}
	//			else { a = 13; }
	//		}
	//		else { a = 23; }
	//	}
	//}
	//else {
	//	if (p.y >= y1)
	//	{
	//		if (p.y <= y2)
	//		{
	//			a = 4;
	//		}
	//		else { a = 14; }
	//	}
	//	else { a = 24; }
	//}
	return a;
}
/*
* 找中点
*/
PixelPoint centerpoint(PixelPoint p0, PixelPoint p1)
{
	PixelPoint center;
	center.x = (p0.x + p1.x) >> 1;
	center.y = (p0.y + p1.y) >> 1;
	return center;
}
/*
* 点是否已经符合范围
* threshold：阈值
*/
bool is_in_threshold(PixelPoint center, int x1, int x2, int y1, int y2, int threshold)
{
	//if ((x1- threshold <= center.x  && center.x <= x1 + threshold) || 
	//	(x2 - threshold <= center.x && center.x <= x2 + threshold))
	//	if ((y1 - threshold <= center.y && center.y <= y1 + threshold) ||
	//		(y2 - threshold <= center.y && center.y <= y2 + threshold))
	//		return true;
	if ((y1 - threshold <= center.y && center.y <= y1 + threshold) ||
		(y2 - threshold <= center.y && center.y <= y2 + threshold))
	{
		if ((x1 - threshold <= center.x && center.x <= x2 + threshold))
			return true;
	}
	else if ((x1 - threshold <= center.x && center.x <= x1 + threshold) ||
		(x2 - threshold <= center.x && center.x <= x2 + threshold))
	{
		if ((y1 - threshold <= center.y && center.y <= y2 + threshold))
			return true;
	}

	//if (center.x >= x2 - threshold && center.x <= x2) return true;
	//if (center.y >= y1 && center.y <= y1 + threshold) return true;
	//if (center.y >= y2 - threshold && center.y <= y2) return true;
	return false;
}
PixelPoint get_farthest_visible_point(PixelPoint p0, PixelPoint p1, int x1, int x2, int y1, int y2, int threshold)
{
	// 寻找p0在p1方向上在范围内的最远可见点
	int code0 = encode(p0, x1, x2, y1, y2);
	int code1 = encode(p1, x1, x2, y1, y2);
	if (code1 == 0)							// p1在窗口内
		return p1;
	else if (code0 & code1)					// p0和p1都在同侧的外面
	{
		return PixelPoint{ -999, -999 };
	}
	else
	{
		PixelPoint temp_p0 = p0, temp_p1 = p1, pm = centerpoint(p0, p1);	// 分别用来表示寻找过程中的临时p0、临时p1和中点
		while (true)
		{
			if (GetLength(temp_p0.x, temp_p0.y, temp_p1.x, temp_p1.y) < threshold)
			{
				if (is_in_threshold(pm, x1, x2, y1, y2, threshold/2))
					return pm;
				else
					return PixelPoint{ -999, -999 };
			}

			if (encode(pm, x1, x2, y1, y2) == 0)				// 如果pm在窗口内
			{
				temp_p0 = pm;
			}
			else if (encode(pm, x1, x2, y1, y2) & encode(temp_p1, x1, x2, y1, y2))	// 如果pm和temp_p1处于同一侧
			{
				temp_p1 = pm;
			}
			pm = centerpoint(temp_p0, temp_p1);
		}
	}
}
////递归
//PixelPoint cut(PixelPoint *p0, PixelPoint *p1, int x1, int x2, int y1, int y2, int threshold)
//{
//	//int code0 = encode(*p0, x1, x2, y1, y2), code1 = encode(*p1, x1, x2, y1, y2);
//	//if (code0)
//	//PixelPoint center;
//	//centerpoint(&center, p0, p1);
//	//int a2 = encode(*p1, x1, x2, y1, y2);
//	//int temp = encode(center, x1, x2, y1, y2);
//	//if (is_in_threshold(*p0, x1, x2, y1, y2, threshold)) {
//	//	p1 = p0;
//	//	return *p1;
//	//}
//	//if (temp == 0)
//	//{
//	//	p0 = &center;
//	//	cut(p0, p1, x1, x2, y1, y2, threshold);
//	//}
//	//else {
//	//	if (is_in_same_outer_zone(temp, a2))
//	//	{
//	//		p1 = &center;
//	//	}
//	//	else {
//	//		p0 = &center;
//	//		
//	//	}
//	//	cut(p0, p1, x1, x2, y1, y2, threshold);
//	//}
//}
/*p为输入的需要裁剪的线段的两端点
* x1左
* x2右
* y1下
* y2上
* threshold阈值
*/
bool centercut(PixelPoint *p, int x1, int x2, int y1, int y2, int threshold)
{
	PixelPoint save1, save2;
	save1 = get_farthest_visible_point(p[0], p[1], x1, x2, y1, y2, threshold);
	save2 = get_farthest_visible_point(p[1], p[0], x1, x2, y1, y2, threshold);
	// 如果返回的结果是-999，说明寻找到的最远点不在窗口内，说明线段不经过窗口
	bool flag = true;
	if (save1.x == -999 && save1.y == -999 || save2.x == -999 && save2.y == -999)
	{
		flag = false;
	}
	else
	{
		p[1] = save1;
		p[0] = save2;
	}

	return flag;



}


bool is_on_segment(Point2D p, Point2D l1, Point2D l2)
{
	if (fabs((p.x - l1.x) * (l2.y - l1.y) - (l2.x - l1.x) * (p.y - l1.y)) < 0.001  //叉乘

		&& (min(l1.x, l2.x)-0.001 <= p.x && p.x <= max(l1.x, l2.x)+0.001)
		//保证p点坐标在l1l2之间
		&& min(l1.y, l2.y) - 0.001 <= p.y && p.y <= max(l1.y, l2.y) + 0.001)
		return true;
	else
		return false;
}

//判断交点是否在线段上
//int is_on_segment(Point2D p, Point2D l1, Point2D l2)
//{
//	
//	if (fabs((p.x - l1.x) * (l2.y - l1.y) - (l2.x - l1.x) * (p.y - l1.y)) < 0.001  //叉乘
//
//		&& (min(l1.x, l2.x) <= p.x && p.x <= max(l1.x, l2.x))
//		//保证p点坐标在l1l2之间
//		&& min(l1.y, l2.y) <= p.y && p.y <= max(l1.y, l2.y))
//		return false;
//	else
//		return true;
//}
//判断两直线是否相交，不相交返回1
int Crosses(Point2D p1, Point2D p2, Point2D s1, Point2D s2)
{
	if ((s1.y - s2.y) * (p1.x - p2.x) - (p1.y - p2.y) * (s1.x - s2.x) != 0) return 1;
	else return 0;
}
/// <summary>
/// 找出多边形顶点中x最大的点，即必定为凸的点
/// </summary>
/// <param name="polygon">多边形顶点集合数组</param>
/// <param name="nodenum">多边形顶点数量</param>
/// <returns>最大值点所在下标</returns>
int FindXmaxPt(Point2D *polygon, int nodenum)
{
	Point2D Xmaxpt;
	Xmaxpt.x = 0; Xmaxpt.y = 0;
	int pos = 0;
	for (int i = 0; i < nodenum; ++i)
	{
		if (Xmaxpt.x < polygon[i].x)
		{
			Xmaxpt.x = polygon[i].x;
			pos = i;
		}
	}
	return pos;
}
/// <summary>
/// 通过确定的凸点下标判断多边形是否为顺时针，不是则逆置输出
/// </summary>
/// <param name="polygon">多边形顶点数组</param>
/// <param name="nodenum">顶点个数</param>
/// <param name="pos">x最大点下标</param>
/// <returns>顺时针顶点数组</returns>
Point2D *JudgeClockWise(Point2D *polygon, int nodenum, int pos)
{
	Point2D xmaxP = polygon[pos];

	Point2D backP;
	if (pos == 0) backP = polygon[nodenum - 1];//特殊情况：多边形顶点数组末尾为第一个顶点特殊考虑；
	else backP = polygon[pos - 1];//最大不可能在末尾的第一个顶点（因为第一个就是它）所以不需要考虑这种情况
	Point2D frontP = polygon[pos + 1];

	Point2D vec1, vec2;
	vec1.x = xmaxP.x - backP.x; vec1.y = xmaxP.y - backP.y;//向量1
	vec2.x = frontP.x - xmaxP.x; vec2.y = frontP.y - xmaxP.y;//向量2
	int cross = vec1.x * vec2.y - vec1.y * vec2.x;

	Point2D *newpoly = new Point2D[nodenum];
	if (cross < 0) return polygon;//叉积为负则为顺时针
	else//叉积为正将其变为顺时针
	{

		for (int i = 0; i < nodenum; ++i)
		{
			newpoly[i] = polygon[nodenum - i - 1];
		}
	}
	return newpoly;
}
/// <summary>
/// 获取水平射线与两点确定线段的交点（用于转角法）
/// </summary>
/// <param name="y">坐标y值</param>
/// <param name="p1">线段端点1</param>
/// <param name="p2">线段端点2</param>
/// <returns>有交点返回坐标x，没有则返回0</returns>
double GetIntersectx(double y, Point2D p1, Point2D p2) {
	double x = ((p2.x - p1.x) * (y - p1.y) / (p2.y - p1.y)) + p1.x;
	if (((x - p1.x) * (x - p2.x)) <= 0) 
		return x;
	else 
		return 0;
}
/// <summary>
/// 转角法判断点是否在多边形内部，返回值不为0则在内部
/// </summary>
/// <param name="pt">点</param>
/// <param name="polygon">多边形点集</param>
/// <returns>返回值</returns>
int ptInPolyByCorner(Point2D pt, vector<Point2D> polygon) {//转角算法判断点是否在多边形内部
	int wn = 0;
	for (int i = 0; i < polygon.size() - 1; ++i) {
		if (polygon[i].y == polygon[i + 1].y) continue;
		if ((pt.y - polygon[i].y) * (pt.y - polygon[i + 1].y) > 0) continue;
		double intersectx = GetIntersectx(pt.y, polygon[i], polygon[i + 1]);

		if (polygon[i + 1].y < polygon[i].y) {
			if (intersectx > pt.x) ++wn;
		}
		else if (polygon[i + 1].y > polygon[i].y)
		{
			if (intersectx > pt.x) --wn;
		}
	}
	return wn;
}
/// <summary>
/// 转角法判断点是否在多边形内部，返回值不为0则在内部
/// </summary>
/// <param name="pt">点</param>
/// <param name="polygon">多边形点集</param>
/// <returns>返回值</returns>
int ptInPolyByCorner(CrossPoint pt, vector<CrossPoint> polygon) 
{//转角算法判断点是否在多边形内部
	int wn = 0;
	for (int i = 0; i < polygon.size() - 1; ++i) {
		if (fabs(polygon[i].y - polygon[i + 1].y) < 0.001) continue;
		if ((pt.y - polygon[i].y) * (pt.y - polygon[i + 1].y) > 0) continue;
		double intersectx = GetIntersectx(pt.y, Point2D{ polygon[i].x, polygon[i].y }, Point2D{ polygon[i + 1].x, polygon[i + 1].y });

		if (polygon[i + 1].y < polygon[i].y) {
			if (intersectx > pt.x) ++wn;
		}
		else if (polygon[i + 1].y > polygon[i].y)
		{
			if (intersectx > pt.x) --wn;
		}
	}
	return wn;
}
/*
当时比较懵该用什么方法
就都写出来了
*/
bool SortIndexp(CrossPoint a, CrossPoint b)
{
	return a.pindex[0] < b.pindex[0];
}
bool SortIndexc(CrossPoint a, CrossPoint b)
{
	return a.pindex[1] < b.pindex[1];
}
void SortClipCrossc(vector<CrossPoint> &crossP)
{
	sort(crossP.begin(), crossP.end(), SortIndexc);
}
void SortClipCrossp(vector<CrossPoint> &crossP)
{
	sort(crossP.begin(), crossP.end(), SortIndexp);
}
vector<CrossPoint> returnSortClipCross(vector<CrossPoint> crossP)
{
	sort(crossP.begin(), crossP.end(), SortIndexc);
	return crossP;
}
/// <summary>
/// 获取交点的pindex在对应融合交点vector中的索引位置
/// </summary>
/// <param name="crossP">交点</param>
/// <returns>更新索引的交点</returns>
/// <summary>
/// 获取交点的pindex在对应融合交点vector中的索引位置
/// </summary>
/// <param name="crossP">交点</param>
/// <returns>更新索引的交点</returns>
vector<CrossPoint> GetIndex(vector<CrossPoint> crossP, vector<Point2D> clip)
{
	int nump = 1;
	int numc = 1;
	for (int i = 0; i < crossP.size(); ++i)
	{
		crossP[i].pindex[0] += nump;
		nump++;
	}
	GetRightOrder(crossP, clip);
	SortClipCrossc(crossP);
	for (int j = 0; j < crossP.size(); j++)
	{
		crossP[j].pindex[1] += numc;
		numc++;
	}
	return crossP;
}
/// <summary>
/// 按照GetIndex处理后得到的交点pindex获得融合交点顶点vector
/// </summary>
/// <param name="cross">交点</param>
/// <param name="unmerge">顶点</param>
/// <param name="flag">0为融合多边形，1为融合窗口</param>
/// <returns>融合的点</returns>
vector<CrossPoint> Merge(vector<CrossPoint> cross, vector<Point2D> unmerge, int flag, vector<Point2D> clip)
{
	cross = GetIndex(cross, clip);
	vector<CrossPoint> merged;
	int i = 0, index = 0, num = 0;

	if (flag == 0) SortClipCrossp(cross);
	for (i; i < cross.size(); ++i)
	{
		int cur = cross[i].pindex[flag];
		for (index; index < cur - num; index++)
		{
			merged.push_back(unmerge[index]);
		}
		merged.push_back(cross[i]);
		++num;
	}
	for (index; index < unmerge.size(); index++)
	{
		merged.push_back(unmerge[index]);
	}
	return merged;
}

vector<CrossPoint> Merge(vector<CrossPoint> &cross, vector<Point2D> &unmerge)
{
	// 这里的unmerged是首尾相连的
	vector<CrossPoint> merged = vector<CrossPoint>(unmerge.begin(), unmerge.end());
	int i = 0, index = 0, num = 0;
	for (i; i < cross.size(); ++i)
	{
		vector<CrossPoint>::iterator insert_iter = merged.begin(), end_iter = merged.end();
		end_iter--; 
		CrossPoint cross_point = cross[i];
		Point2D cross_point_point2d = Point2D(cross_point.x, cross_point.y);
		for (; insert_iter != end_iter; insert_iter++)
		{
			Point2D l1 = Point2D((*insert_iter).x, (*insert_iter).y);
			insert_iter++;
			Point2D l2 = Point2D((*insert_iter).x, (*insert_iter).y);
			insert_iter--;
			if (is_on_segment(cross_point_point2d, l1, l2))
			{
				insert_iter++;
				break;
			}
			// 暂时没有考虑自交然后有一个交点正好处于自交点的情况
		}
		merged.insert(insert_iter, cross_point);
	}
	return merged;
}
vector<vector<CrossPoint>> UpdateIndex(vector<CrossPoint> &pg_vec, vector<CrossPoint> &clip_vec, vector<CrossPoint> cross_points)
{
	vector<vector<CrossPoint>> result;
	vector<CrossPoint> pg_with_cross, clip_with_cross;
	int flag = 0;
	if (ptInPolyByCorner(pg_vec[0], clip_vec) == 0) flag = 1;	//如果第一个点在外面那么后面如果有交点第一个必定是入点所以flag为1
	for (auto &poly_point : pg_vec)
	{
		if (poly_point.flag == -1)
			continue;
		poly_point.flag = flag;
		flag = !flag;
		for (int j = 0; j < clip_vec.size(); j++)//auto &clip_point : clip_vec)
		{
			CrossPoint clip_point = clip_vec[j];
			if (poly_point == clip_point && clip_point.flag != -1)
			{
				poly_point.pindex[1] = j;
				break;
			}
		}
	}
	for (auto &clip_point : clip_vec)
	{
		if (clip_point.flag == -1)
			continue;
		for (int j = 0; j < pg_vec.size(); j++)//auto &clip_point : clip_vec)
		{
			CrossPoint poly_point = pg_vec[j];
			if (poly_point == clip_point)
			{
				clip_point.pindex[0] = j;
				clip_point.flag = poly_point.flag;
				break;
			}
		}
	}
	return result;
}
/// <summary>
/// 获取裁剪窗口与多边形的交点及其出入状态和位置信息
/// </summary>
/// <param name="pg_vec">多边形点</param>
/// <param name="clip_vec">窗口点</param>
/// <returns>交点</returns>
vector<CrossPoint> GetCrossPts(vector<Point2D> &pg_vec, vector<Point2D> &clip_vec)
{
	int pgnum = pg_vec.size();
	int clipnum = clip_vec.size();
	int flag = 0;
	vector<CrossPoint> crossPts;

	//if (ptInPolyByCorner(pg_vec[0], clip_vec) == 0) flag = 1;//如果第一个点在外面那么后面如果有交点第一个必定是入点所以flag为1

	for (int i = 0; i < pgnum - 1; ++i)
	{
		for (int j = 0; j < clipnum - 1; ++j)
		{
			if (!Crosses(pg_vec[i], pg_vec[i + 1], clip_vec[j], clip_vec[j + 1])) continue;
			Point2D across = GetCrossPoint(pg_vec[i], pg_vec[i + 1], clip_vec[j], clip_vec[j + 1]);
			if (!is_on_segment(across, pg_vec[i], pg_vec[i + 1])) continue;
			if (!is_on_segment(across, clip_vec[j], clip_vec[j + 1])) continue;
			//if ((across.x - pg_vec[i].x) * (across.x - pg_vec[i + 1].x) >= 0) continue;
			//if ((across.x - clip_vec[j].x) * (across.x - clip_vec[j + 1].x) > 0) continue;
			CrossPoint acrp(across, flag, -1, -1);//创建一个交点
			crossPts.push_back(acrp);
			//flag = !flag;//下一个交点一定与本次交点为反向
			// 这边flag大部分时候没问题，但是有的时候会出错
			// 比如http://hi.csdn.net/attachment/201005/24/0_1274687832F3b5.gif这张图片中的ef的flag就会翻转
		}
	}
	return crossPts;
}

int find_in_point_in_clip(vector<CrossPoint> &clipandcross, int start_index, vector<CrossPoint> &aclipped)
{
	//for (int i = start_index; i < clipandcross.size(); i++)
	int i = start_index % clipandcross.size();
	while (true)
	{
		if (clipandcross[i].flag != 1)
		{
			aclipped.push_back(clipandcross[i]);
		}
		else
		{
			return clipandcross[i].pindex[0];
		}
		i++;
		if (i == clipandcross.size())
			i %= clipandcross.size();
	}
}

int find_out_point_in_poly(vector<CrossPoint> &pgandcross, int start_index, vector<CrossPoint> &aclipped, CrossPoint poly_start)
{
	//for (int i = start_index; i < pgandcross.size(); i++)
	int i = start_index % pgandcross.size();
	while (true)
	{
		if (pgandcross[i] == poly_start)
		{
			aclipped.push_back(pgandcross[i]);
			return -1;
		}
		if (pgandcross[i].flag != 0)
		{
			pgandcross[i].flag = -1;
			aclipped.push_back(pgandcross[i]);
		}
		else
		{
			return pgandcross[i].pindex[1];
		}
		i++;
		if (i == pgandcross.size())
			i %= pgandcross.size();
	}
}


//vector<CrossPoint> func(vector<CrossPoint> &pgandcross, vector<CrossPoint> &clipandcross)
//{
//	vector<CrossPoint> aclipped;
//	CrossPoint stpt;
//	vector<int> tracked_indices;
//	for (int index = 0; index < pgandcross.size(); index++)
//	{
//		if (pgandcross[index].flag != 1)	// 寻找入点
//			continue;						// 如果这个点不是那就到下一个
//		// 找到一个入点了
//		//in_found = true;
//		stpt = pgandcross[index];			// 设为起点
//		tracked_indices.push_back(index);
//		//aclipped.push_back(stpt);			// 加入第一个结果，下面开始寻找出点
//		for (int outindex = index; outindex < pgandcross.size(); outindex++)
//		{
//			if (pgandcross[outindex].flag != 0)
//			{
//				aclipped.push_back(pgandcross[outindex]);	// 如果不是出点，那么就加到结果里面
//				tracked_indices.push_back(outindex);
//				continue;									// 到下一个点
//			}
//			// 找到一个出点了
//			//aclipped.push_back(pgandcross[outindex]);		// 出点加到结果里
//			int inindex_in_clip = pgandcross[outindex].pindex[1];	// 这个出点在裁剪窗口中的index，下面开始寻找入点
//			for (; inindex_in_clip < clipandcross.size(); inindex_in_clip++)
//			{
//				
//			}
//		}
//	}
//	return aclipped;
//}
/// <summary>
/// 应用WA算法进行裁剪
/// </summary>
/// <param name="pg_vec">顺时针多边形顶点</param>
/// <param name="clip_vec">顺时针裁剪窗口顶点</param>
/// <returns>装有得到裁剪多边形顶点vector的vector</returns>
vector<vector<CrossPoint>> WA_Clip(vector<Point2D> pg_vec, vector<Point2D> clip_vec)
{
	//pg_vec.clear();
	//clip_vec.clear();
	//pg_vec.push_back(Point2D(297.74444444444441, 431.00000000000000));
	//pg_vec.push_back(Point2D(553.31081081081061, 431.00000000000000));
	//pg_vec.push_back(Point2D(598.91891891891896, 306.00000000000000));
	//pg_vec.push_back(Point2D(248.43888888888890, 306.00000000000000));
	//pg_vec.push_back(Point2D(297.74444444444441, 431.00000000000000));
	//clip_vec.push_back(Point2D(476.00000000000000, 256.00000000000000));
	//clip_vec.push_back(Point2D(440.00000000000000, 471.00000000000000));
	//clip_vec.push_back(Point2D(674.00000000000000, 458.00000000000000));
	//clip_vec.push_back(Point2D(728.00000000000000, 269.00000000000000));
	//clip_vec.push_back(Point2D(476.00000000000000, 256.00000000000000));
	//clip_vec.push_back(Point2D(840.00000000000000, 488.00000000000000));
	//clip_vec.push_back(Point2D(849.00000000000000, 334.00000000000000));
	//{x = 297.74444444444441, y = 431.00000000000000 }
	//{x = 553.31081081081061, y = 431.00000000000000 }
	//{x = 598.91891891891896, y = 306.00000000000000 }
	//{x = 248.43888888888890, y = 306.00000000000000 }
	//{x = 297.74444444444441, y = 431.00000000000000 }

	//{x = 476.00000000000000 y = 256.00000000000000 }
	//{x = 440.00000000000000 y = 471.00000000000000 }
	//{x = 674.00000000000000 y = 458.00000000000000 }
	//{x = 728.00000000000000 y = 269.00000000000000 }
	//{x = 476.00000000000000 y = 256.00000000000000 }
	vector<vector<CrossPoint>> output;//输出
	vector<CrossPoint> crossPts = GetCrossPts(pg_vec, clip_vec);//获得交点
	vector<CrossPoint> pgandcross = Merge(crossPts, pg_vec);
	vector<CrossPoint> clipandcross = Merge(crossPts, clip_vec);
	UpdateIndex(pgandcross, clipandcross, crossPts);
	//vector<CrossPoint> pgandcross = Merge(crossPts, pg_vec, 0, clip_vec);//融合交点和窗口
	//vector<CrossPoint> clipandcross = Merge(crossPts, clip_vec, 1, clip_vec);//融合交点和多边形顶点
	// e和f对调了，导致报错
	pgandcross.pop_back();
	clipandcross.pop_back();	// 这个实现需要首尾不相连
	vector<CrossPoint> aclipped;//一个被裁剪多边形
	CrossPoint stpt;//起始点

	int pos = 0;
	//pgandcross.clear();
	//clipandcross.clear();
	//pgandcross.push_back(CrossPoint(745.52277362338543, 445.89123045547245, -1, 0, 0));
	//pgandcross.push_back(CrossPoint(786.00000000000000, 443.00000000000000, -1, 0, 0));
	//pgandcross.push_back(CrossPoint(816.90446773992267, 336.68863097466618, 1, -1, 1));
	//pgandcross.push_back(CrossPoint(836.00000000000000, 271.00000000000000, -1, 0, 0));
	//pgandcross.push_back(CrossPoint(720.19100499722379, 270.25763464741812, -1, 0, 0));
	//pgandcross.push_back(CrossPoint(730.81257461201756, 343.90051730998806, 0, -1, 2));	// 这两个多边形的0和1反了
	//clipandcross.push_back(CrossPoint(849.00000000000000, 334.00000000000000, -1, 0, 0));
	//clipandcross.push_back(CrossPoint(816.90446773992267, 336.68863097466618, 1, 2, -1));
	//clipandcross.push_back(CrossPoint(730.81257461201756, 343.90051730998806, 0, 5, -1));
	//clipandcross.push_back(CrossPoint(658.00000000000000, 350.00000000000000, -1, 0, 0));
	//clipandcross.push_back(CrossPoint(719.00000000000000, 483.00000000000000, -1, 0, 0));
	//clipandcross.push_back(CrossPoint(840.00000000000000, 488.00000000000000, -1, 0, 0));
	for (int index = 0; index < pgandcross.size(); index++)
	{
		if (pgandcross[index].flag != 1)	// 寻找入点
			continue;						// 如果这个点不是那就到下一个
		// 找到一个入点了
		stpt = pgandcross[index];			// 设为起点
		aclipped.push_back(stpt);
		int start_index = index + 1;
		while (true)
		{
			start_index = find_out_point_in_poly(pgandcross, start_index, aclipped, stpt);
			if (start_index == -1)	// 回到起点了
			{
				output.push_back(aclipped);
				aclipped.clear();
				break;
			}
			start_index = find_in_point_in_clip(clipandcross, start_index, aclipped);
		}
		////aclipped.push_back(stpt);			// 加入第一个结果，下面开始寻找出点
		//for (int outindex = index; outindex < pgandcross.size(); outindex++)
		//{
		//	if (pgandcross[outindex].flag != 0)
		//	{
		//		aclipped.push_back(pgandcross[outindex]);	// 如果不是出点，那么就加到结果里面
		//		pgandcross[outindex].flag = -1;
		//		continue;									// 到下一个点
		//	}
		//	// 找到一个出点了
		//	//aclipped.push_back(pgandcross[outindex]);		// 出点加到结果里
		//	int inindex_in_clip = pgandcross[outindex].pindex[1];	// 这个出点在裁剪窗口中的index，下面开始寻找入点
		//	for (; inindex_in_clip < clipandcross.size(); inindex_in_clip++)
		//	{
		//		if (pgandcross[inindex_in_clip].flag != 1)
		//		{
		//			aclipped.push_back(pgandcross[inindex_in_clip]);	// 如果不是入点，那么就加到结果里面
		//			pgandcross[inindex_in_clip].flag = -1;
		//			continue;									// 到下一个点
		//		}
		//		// 找到一个入点了
		//		if (clipandcross[inindex_in_clip] == stpt)		//如果为起点则将当前裁剪加入输出，寻找是否有其他被裁开的多边形
		//		{
		//			output.push_back(aclipped);
		//			break;
		//		}
		//	}
		//}
	}
	// todo: 置空可能会出问题



	//for (index; index < pgandcross.size(); ++index)
	//{
	//	if (pgandcross[index].flag == 1)//寻找一个入点
	//	{
	//		stpt = pgandcross[index];//设为起点
	//		aclipped.push_back(stpt);//加入第一个结果

	//		int outindex = index + 1;//存储找到出点的索引
	//		for (outindex; outindex < pgandcross.size(); ++outindex)//寻找出点
	//		{
	//			aclipped.push_back(pgandcross[outindex]);
	//			if (pgandcross[outindex].flag == 0)
	//			{
	//				outindex = pgandcross[outindex].pindex[1];//获得其在裁剪窗口中的索引位置
	//				int inindex = outindex + 1;//存储找到入点的索引
	//				for (inindex; inindex < clipandcross.size(); ++inindex)//寻找入点
	//				{
	//					aclipped.push_back(clipandcross[inindex]);
	//					if (clipandcross[inindex].flag == 1) break;
	//				}

	//				if (clipandcross[inindex] == stpt)//如果为起点则将当前裁剪加入输出，寻找是否有其他被裁开的多边形
	//				{
	//					output.push_back(aclipped);
	//					break;
	//				}
	//			}
	//		}
	//	}
	//}
	return output;
}







/// <summary>
/// 判断多边形顶点是否为逆时针若不是则改为逆时针
/// </summary>
/// <param name="poly">多边形顶点数组</param>
/// <param name="num">多边形顶点个数</param>
/// <returns>多边形逆时针顶点数组</returns>
void ChangeClockwise(Point2D *poly, int num)
{
	int cross = (poly[0].x - poly[2].x) * (poly[1].y - poly[2].y) - (poly[1].x - poly[2].x) * (poly[0].y - poly[2].y);
	int flag = (cross > 0) ? 0 : 1;
	
	if (flag)
	{
		Point2D temp;
		//Point2D *newpoly = new Point2D[num];
		for (int i = 0; i < (num+1)/2; i++)
		{
			temp = poly[i];
			poly[i] = poly[num - 1 - i];
			poly[num - 1 - i] = temp;
			//newpoly[i] = poly[num - 1 - i];
		}
	}
}
/// <summary>
/// 通过求叉积的方式判断点在向量的内侧或外侧（顺序为逆时针方向）
/// </summary>
/// <param name="s">起始点</param>
/// <param name="e">终止点</param>
/// <param name="p">需要判断的点</param>
/// <returns>外侧为0，内侧为1</returns>
int RTInside(Point2D s, Point2D e, Point2D p)//叉积判断，外侧为0，内侧为1
{
	int cross = (e.x - s.x) * (p.y - s.y) - (p.x - s.x) * (e.y - s.y);
	int flag = (cross > 0) ? 1 : 0;
	return flag;
}
/// <summary>
/// 通过四个点得到所在两条直线的交点
/// </summary>
/// <param name="p1">第一条直线的一个点</param>
/// <param name="p2">第一条直线的另一个点</param>
/// <param name="s1">第二条直线的一个点</param>
/// <param name="s2">第二条直线的另一个点</param>
/// <returns>交点</returns>
Point2D GetCrossPoint(Point2D p1, Point2D p2, Point2D s1, Point2D s2)
{
	double par3 = (s1.y - s2.y) * (p1.x - p2.x) - (p1.y - p2.y) * (s1.x - s2.x);
	double par2 = p1.y * p2.x - p1.x * p2.y;
	double par1 = s1.y * s2.x - s1.x * s2.y;
	double x = ((p1.x - p2.x) * par1 - (s1.x - s2.x) * par2) / par3;
	double y = ((p1.y - p2.y) * par1 - (s1.y - s2.y) * par2) / par3;
	Point2D crosspoint;
	crosspoint.x = x;
	crosspoint.y = y;
	return crosspoint;
}
/// <summary>
/// SH裁剪算法
/// </summary>
/// <param name="polygon">待裁剪多边形顶点vector</param>
/// <param name="clipbound">裁剪边框顶点vector</param>
/// <returns>裁剪完毕的多边形vector</returns>
vector<Point2D> subSH_Clip(vector<Point2D> pg_vec, vector<Point2D> clip_vec)
{
	if (pg_vec.empty())
		return pg_vec;
	int i = 0, flag = 0;
	Point2D start, stop;//被剪裁多边形的边向量起点和终点
	Point2D sp, ep;//剪裁窗口边界向量的起点和终点
	vector<Point2D> temp;

	sp = clip_vec[clip_vec.size() - 1];

	for (i; i < clip_vec.size(); i++)
	{
		ep = clip_vec[i];
		start = pg_vec[pg_vec.size() - 1];
		flag = RTInside(sp, ep, start);

		for (int j = 0; j < pg_vec.size(); j++)
		{
			stop = pg_vec[j];
			if (RTInside(sp, ep, stop))//当前第i个顶点是否在边界内侧
			{
				if (flag == 0)//前一个点是否在外侧
				{
					flag = 1;//由外到内的情况					
					temp.push_back(GetCrossPoint(sp, ep, start, stop));//存入交点
				}
				temp.push_back(stop);
			}
			else
			{
				if (flag)
				{
					flag = 0;//由内到外
					temp.push_back(GetCrossPoint(sp, ep, start, stop));//存入交点
				}
			}
			start = stop;
		}
		sp = ep;
		if (temp.size() == 0) return temp;
		pg_vec = temp;//将经过一条边裁剪的vector继续裁剪
		temp.clear();//清空中间变量用于下一次裁剪
	}
	pg_vec.push_back(pg_vec[0]);
	return pg_vec;
}

