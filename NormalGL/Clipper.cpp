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
/// ������࣬�̳�Point2D����������ֱ�Ӽ̳��˽ṹ��Ҳ��֪����û������
/// </summary>
class CrossPoint :public Point2D {
public:
	double x, y;
	int pindex[2] = { 0,0 };//��Ӧ�Ķ���α�,��Ӧ�Ĳü�����εı�
	int flag = 0;//flag����Ϊ0�����Ϊ1
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

	CrossPoint(Point2D pt)//���ڴ洢�ǽ��㣬����flagΪ-1����
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
/// ����ͬһ���ü����ڱ��ϵĵ����˳ʱ������
/// </summary>
/// <param name="cross">����</param>
/// <param name="clip">�ü���</param>
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
	// points��Ҫ��β����
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
	bool flag = false; // ������ע��û�м��óɹ�
	for (int i = 0; i < painter->GetPointCount() - 1; i++)
	{
		PixelPoint p[2];
		p[0] = PixelPoint{ (int)painter->GetX(i), (int)painter->GetY(i) };
		p[1] = PixelPoint{ (int)painter->GetX(i + 1), (int)painter->GetY(i + 1) };
		bool temp_flag = centercut(p, boundary[0].x, boundary[3].x, boundary[0].y, boundary[1].y, 2);	// �����false����ô˵������߶βü����Ѿ�û��
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
	to_anticlockwise(pg, pgnum);			// ����ʵ�ֵ�sh�㷨��ʵ�ǲ���Ҫ�������β�����ģ��������������Ҫ������ǰ������˶�һ���ռ���ڴ棬��ߴ�����ǵ������+1
	to_anticlockwise(clip, clipnum);		// �����Ҵ�����λ㱨֮��ѵ������ͳһ�ɴ������һ�������������������ʵ�������
	vector<Point2D> pgvec(pg, pg + pgnum - 1);
	vector<Point2D> clipvec(clip, clip + clipnum - 1);

	vector<Point2D> clipped;				// ����������β������
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
		return 0;			// ��������ö������ȫ�ڼ��ö�����ھ�ֱ�ӷ��أ���0�����������
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
		return 1;			// ������ö������ȫ�ڱ����ö�����ھ�ֱ�ӷ��أ���1�����������
	}


	vector<vector<CrossPoint>> clippedvector = WA_Clip(pg_vec, clip_vec);	// �������β������
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



/*������������
�㣺p
x1����
x2����
y1����
y2����*/
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
* ���е�
*/
PixelPoint centerpoint(PixelPoint p0, PixelPoint p1)
{
	PixelPoint center;
	center.x = (p0.x + p1.x) >> 1;
	center.y = (p0.y + p1.y) >> 1;
	return center;
}
/*
* ���Ƿ��Ѿ����Ϸ�Χ
* threshold����ֵ
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
	// Ѱ��p0��p1�������ڷ�Χ�ڵ���Զ�ɼ���
	int code0 = encode(p0, x1, x2, y1, y2);
	int code1 = encode(p1, x1, x2, y1, y2);
	if (code1 == 0)							// p1�ڴ�����
		return p1;
	else if (code0 & code1)					// p0��p1����ͬ�������
	{
		return PixelPoint{ -999, -999 };
	}
	else
	{
		PixelPoint temp_p0 = p0, temp_p1 = p1, pm = centerpoint(p0, p1);	// �ֱ�������ʾѰ�ҹ����е���ʱp0����ʱp1���е�
		while (true)
		{
			if (GetLength(temp_p0.x, temp_p0.y, temp_p1.x, temp_p1.y) < threshold)
			{
				if (is_in_threshold(pm, x1, x2, y1, y2, threshold/2))
					return pm;
				else
					return PixelPoint{ -999, -999 };
			}

			if (encode(pm, x1, x2, y1, y2) == 0)				// ���pm�ڴ�����
			{
				temp_p0 = pm;
			}
			else if (encode(pm, x1, x2, y1, y2) & encode(temp_p1, x1, x2, y1, y2))	// ���pm��temp_p1����ͬһ��
			{
				temp_p1 = pm;
			}
			pm = centerpoint(temp_p0, temp_p1);
		}
	}
}
////�ݹ�
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
/*pΪ�������Ҫ�ü����߶ε����˵�
* x1��
* x2��
* y1��
* y2��
* threshold��ֵ
*/
bool centercut(PixelPoint *p, int x1, int x2, int y1, int y2, int threshold)
{
	PixelPoint save1, save2;
	save1 = get_farthest_visible_point(p[0], p[1], x1, x2, y1, y2, threshold);
	save2 = get_farthest_visible_point(p[1], p[0], x1, x2, y1, y2, threshold);
	// ������صĽ����-999��˵��Ѱ�ҵ�����Զ�㲻�ڴ����ڣ�˵���߶β���������
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
	if (fabs((p.x - l1.x) * (l2.y - l1.y) - (l2.x - l1.x) * (p.y - l1.y)) < 0.001  //���

		&& (min(l1.x, l2.x)-0.001 <= p.x && p.x <= max(l1.x, l2.x)+0.001)
		//��֤p��������l1l2֮��
		&& min(l1.y, l2.y) - 0.001 <= p.y && p.y <= max(l1.y, l2.y) + 0.001)
		return true;
	else
		return false;
}

//�жϽ����Ƿ����߶���
//int is_on_segment(Point2D p, Point2D l1, Point2D l2)
//{
//	
//	if (fabs((p.x - l1.x) * (l2.y - l1.y) - (l2.x - l1.x) * (p.y - l1.y)) < 0.001  //���
//
//		&& (min(l1.x, l2.x) <= p.x && p.x <= max(l1.x, l2.x))
//		//��֤p��������l1l2֮��
//		&& min(l1.y, l2.y) <= p.y && p.y <= max(l1.y, l2.y))
//		return false;
//	else
//		return true;
//}
//�ж���ֱ���Ƿ��ཻ�����ཻ����1
int Crosses(Point2D p1, Point2D p2, Point2D s1, Point2D s2)
{
	if ((s1.y - s2.y) * (p1.x - p2.x) - (p1.y - p2.y) * (s1.x - s2.x) != 0) return 1;
	else return 0;
}
/// <summary>
/// �ҳ�����ζ�����x���ĵ㣬���ض�Ϊ͹�ĵ�
/// </summary>
/// <param name="polygon">����ζ��㼯������</param>
/// <param name="nodenum">����ζ�������</param>
/// <returns>���ֵ�������±�</returns>
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
/// ͨ��ȷ����͹���±��ж϶�����Ƿ�Ϊ˳ʱ�룬�������������
/// </summary>
/// <param name="polygon">����ζ�������</param>
/// <param name="nodenum">�������</param>
/// <param name="pos">x�����±�</param>
/// <returns>˳ʱ�붥������</returns>
Point2D *JudgeClockWise(Point2D *polygon, int nodenum, int pos)
{
	Point2D xmaxP = polygon[pos];

	Point2D backP;
	if (pos == 0) backP = polygon[nodenum - 1];//�������������ζ�������ĩβΪ��һ���������⿼�ǣ�
	else backP = polygon[pos - 1];//��󲻿�����ĩβ�ĵ�һ�����㣨��Ϊ��һ�������������Բ���Ҫ�����������
	Point2D frontP = polygon[pos + 1];

	Point2D vec1, vec2;
	vec1.x = xmaxP.x - backP.x; vec1.y = xmaxP.y - backP.y;//����1
	vec2.x = frontP.x - xmaxP.x; vec2.y = frontP.y - xmaxP.y;//����2
	int cross = vec1.x * vec2.y - vec1.y * vec2.x;

	Point2D *newpoly = new Point2D[nodenum];
	if (cross < 0) return polygon;//���Ϊ����Ϊ˳ʱ��
	else//���Ϊ�������Ϊ˳ʱ��
	{

		for (int i = 0; i < nodenum; ++i)
		{
			newpoly[i] = polygon[nodenum - i - 1];
		}
	}
	return newpoly;
}
/// <summary>
/// ��ȡˮƽ����������ȷ���߶εĽ��㣨����ת�Ƿ���
/// </summary>
/// <param name="y">����yֵ</param>
/// <param name="p1">�߶ζ˵�1</param>
/// <param name="p2">�߶ζ˵�2</param>
/// <returns>�н��㷵������x��û���򷵻�0</returns>
double GetIntersectx(double y, Point2D p1, Point2D p2) {
	double x = ((p2.x - p1.x) * (y - p1.y) / (p2.y - p1.y)) + p1.x;
	if (((x - p1.x) * (x - p2.x)) <= 0) 
		return x;
	else 
		return 0;
}
/// <summary>
/// ת�Ƿ��жϵ��Ƿ��ڶ�����ڲ�������ֵ��Ϊ0�����ڲ�
/// </summary>
/// <param name="pt">��</param>
/// <param name="polygon">����ε㼯</param>
/// <returns>����ֵ</returns>
int ptInPolyByCorner(Point2D pt, vector<Point2D> polygon) {//ת���㷨�жϵ��Ƿ��ڶ�����ڲ�
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
/// ת�Ƿ��жϵ��Ƿ��ڶ�����ڲ�������ֵ��Ϊ0�����ڲ�
/// </summary>
/// <param name="pt">��</param>
/// <param name="polygon">����ε㼯</param>
/// <returns>����ֵ</returns>
int ptInPolyByCorner(CrossPoint pt, vector<CrossPoint> polygon) 
{//ת���㷨�жϵ��Ƿ��ڶ�����ڲ�
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
��ʱ�Ƚ��¸���ʲô����
�Ͷ�д������
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
/// ��ȡ�����pindex�ڶ�Ӧ�ںϽ���vector�е�����λ��
/// </summary>
/// <param name="crossP">����</param>
/// <returns>���������Ľ���</returns>
/// <summary>
/// ��ȡ�����pindex�ڶ�Ӧ�ںϽ���vector�е�����λ��
/// </summary>
/// <param name="crossP">����</param>
/// <returns>���������Ľ���</returns>
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
/// ����GetIndex�����õ��Ľ���pindex����ںϽ��㶥��vector
/// </summary>
/// <param name="cross">����</param>
/// <param name="unmerge">����</param>
/// <param name="flag">0Ϊ�ں϶���Σ�1Ϊ�ںϴ���</param>
/// <returns>�ںϵĵ�</returns>
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
	// �����unmerged����β������
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
			// ��ʱû�п����Խ�Ȼ����һ���������ô����Խ�������
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
	if (ptInPolyByCorner(pg_vec[0], clip_vec) == 0) flag = 1;	//�����һ������������ô��������н����һ���ض����������flagΪ1
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
/// ��ȡ�ü����������εĽ��㼰�����״̬��λ����Ϣ
/// </summary>
/// <param name="pg_vec">����ε�</param>
/// <param name="clip_vec">���ڵ�</param>
/// <returns>����</returns>
vector<CrossPoint> GetCrossPts(vector<Point2D> &pg_vec, vector<Point2D> &clip_vec)
{
	int pgnum = pg_vec.size();
	int clipnum = clip_vec.size();
	int flag = 0;
	vector<CrossPoint> crossPts;

	//if (ptInPolyByCorner(pg_vec[0], clip_vec) == 0) flag = 1;//�����һ������������ô��������н����һ���ض����������flagΪ1

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
			CrossPoint acrp(across, flag, -1, -1);//����һ������
			crossPts.push_back(acrp);
			//flag = !flag;//��һ������һ���뱾�ν���Ϊ����
			// ���flag�󲿷�ʱ��û���⣬�����е�ʱ������
			// ����http://hi.csdn.net/attachment/201005/24/0_1274687832F3b5.gif����ͼƬ�е�ef��flag�ͻᷭת
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
//		if (pgandcross[index].flag != 1)	// Ѱ�����
//			continue;						// �������㲻���Ǿ͵���һ��
//		// �ҵ�һ�������
//		//in_found = true;
//		stpt = pgandcross[index];			// ��Ϊ���
//		tracked_indices.push_back(index);
//		//aclipped.push_back(stpt);			// �����һ����������濪ʼѰ�ҳ���
//		for (int outindex = index; outindex < pgandcross.size(); outindex++)
//		{
//			if (pgandcross[outindex].flag != 0)
//			{
//				aclipped.push_back(pgandcross[outindex]);	// ������ǳ��㣬��ô�ͼӵ��������
//				tracked_indices.push_back(outindex);
//				continue;									// ����һ����
//			}
//			// �ҵ�һ��������
//			//aclipped.push_back(pgandcross[outindex]);		// ����ӵ������
//			int inindex_in_clip = pgandcross[outindex].pindex[1];	// ��������ڲü������е�index�����濪ʼѰ�����
//			for (; inindex_in_clip < clipandcross.size(); inindex_in_clip++)
//			{
//				
//			}
//		}
//	}
//	return aclipped;
//}
/// <summary>
/// Ӧ��WA�㷨���вü�
/// </summary>
/// <param name="pg_vec">˳ʱ�����ζ���</param>
/// <param name="clip_vec">˳ʱ��ü����ڶ���</param>
/// <returns>װ�еõ��ü�����ζ���vector��vector</returns>
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
	vector<vector<CrossPoint>> output;//���
	vector<CrossPoint> crossPts = GetCrossPts(pg_vec, clip_vec);//��ý���
	vector<CrossPoint> pgandcross = Merge(crossPts, pg_vec);
	vector<CrossPoint> clipandcross = Merge(crossPts, clip_vec);
	UpdateIndex(pgandcross, clipandcross, crossPts);
	//vector<CrossPoint> pgandcross = Merge(crossPts, pg_vec, 0, clip_vec);//�ںϽ���ʹ���
	//vector<CrossPoint> clipandcross = Merge(crossPts, clip_vec, 1, clip_vec);//�ںϽ���Ͷ���ζ���
	// e��f�Ե��ˣ����±���
	pgandcross.pop_back();
	clipandcross.pop_back();	// ���ʵ����Ҫ��β������
	vector<CrossPoint> aclipped;//һ�����ü������
	CrossPoint stpt;//��ʼ��

	int pos = 0;
	//pgandcross.clear();
	//clipandcross.clear();
	//pgandcross.push_back(CrossPoint(745.52277362338543, 445.89123045547245, -1, 0, 0));
	//pgandcross.push_back(CrossPoint(786.00000000000000, 443.00000000000000, -1, 0, 0));
	//pgandcross.push_back(CrossPoint(816.90446773992267, 336.68863097466618, 1, -1, 1));
	//pgandcross.push_back(CrossPoint(836.00000000000000, 271.00000000000000, -1, 0, 0));
	//pgandcross.push_back(CrossPoint(720.19100499722379, 270.25763464741812, -1, 0, 0));
	//pgandcross.push_back(CrossPoint(730.81257461201756, 343.90051730998806, 0, -1, 2));	// ����������ε�0��1����
	//clipandcross.push_back(CrossPoint(849.00000000000000, 334.00000000000000, -1, 0, 0));
	//clipandcross.push_back(CrossPoint(816.90446773992267, 336.68863097466618, 1, 2, -1));
	//clipandcross.push_back(CrossPoint(730.81257461201756, 343.90051730998806, 0, 5, -1));
	//clipandcross.push_back(CrossPoint(658.00000000000000, 350.00000000000000, -1, 0, 0));
	//clipandcross.push_back(CrossPoint(719.00000000000000, 483.00000000000000, -1, 0, 0));
	//clipandcross.push_back(CrossPoint(840.00000000000000, 488.00000000000000, -1, 0, 0));
	for (int index = 0; index < pgandcross.size(); index++)
	{
		if (pgandcross[index].flag != 1)	// Ѱ�����
			continue;						// �������㲻���Ǿ͵���һ��
		// �ҵ�һ�������
		stpt = pgandcross[index];			// ��Ϊ���
		aclipped.push_back(stpt);
		int start_index = index + 1;
		while (true)
		{
			start_index = find_out_point_in_poly(pgandcross, start_index, aclipped, stpt);
			if (start_index == -1)	// �ص������
			{
				output.push_back(aclipped);
				aclipped.clear();
				break;
			}
			start_index = find_in_point_in_clip(clipandcross, start_index, aclipped);
		}
		////aclipped.push_back(stpt);			// �����һ����������濪ʼѰ�ҳ���
		//for (int outindex = index; outindex < pgandcross.size(); outindex++)
		//{
		//	if (pgandcross[outindex].flag != 0)
		//	{
		//		aclipped.push_back(pgandcross[outindex]);	// ������ǳ��㣬��ô�ͼӵ��������
		//		pgandcross[outindex].flag = -1;
		//		continue;									// ����һ����
		//	}
		//	// �ҵ�һ��������
		//	//aclipped.push_back(pgandcross[outindex]);		// ����ӵ������
		//	int inindex_in_clip = pgandcross[outindex].pindex[1];	// ��������ڲü������е�index�����濪ʼѰ�����
		//	for (; inindex_in_clip < clipandcross.size(); inindex_in_clip++)
		//	{
		//		if (pgandcross[inindex_in_clip].flag != 1)
		//		{
		//			aclipped.push_back(pgandcross[inindex_in_clip]);	// ���������㣬��ô�ͼӵ��������
		//			pgandcross[inindex_in_clip].flag = -1;
		//			continue;									// ����һ����
		//		}
		//		// �ҵ�һ�������
		//		if (clipandcross[inindex_in_clip] == stpt)		//���Ϊ����򽫵�ǰ�ü����������Ѱ���Ƿ����������ÿ��Ķ����
		//		{
		//			output.push_back(aclipped);
		//			break;
		//		}
		//	}
		//}
	}
	// todo: �ÿտ��ܻ������



	//for (index; index < pgandcross.size(); ++index)
	//{
	//	if (pgandcross[index].flag == 1)//Ѱ��һ�����
	//	{
	//		stpt = pgandcross[index];//��Ϊ���
	//		aclipped.push_back(stpt);//�����һ�����

	//		int outindex = index + 1;//�洢�ҵ����������
	//		for (outindex; outindex < pgandcross.size(); ++outindex)//Ѱ�ҳ���
	//		{
	//			aclipped.push_back(pgandcross[outindex]);
	//			if (pgandcross[outindex].flag == 0)
	//			{
	//				outindex = pgandcross[outindex].pindex[1];//������ڲü������е�����λ��
	//				int inindex = outindex + 1;//�洢�ҵ���������
	//				for (inindex; inindex < clipandcross.size(); ++inindex)//Ѱ�����
	//				{
	//					aclipped.push_back(clipandcross[inindex]);
	//					if (clipandcross[inindex].flag == 1) break;
	//				}

	//				if (clipandcross[inindex] == stpt)//���Ϊ����򽫵�ǰ�ü����������Ѱ���Ƿ����������ÿ��Ķ����
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
/// �ж϶���ζ����Ƿ�Ϊ��ʱ�����������Ϊ��ʱ��
/// </summary>
/// <param name="poly">����ζ�������</param>
/// <param name="num">����ζ������</param>
/// <returns>�������ʱ�붥������</returns>
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
/// ͨ�������ķ�ʽ�жϵ����������ڲ����ࣨ˳��Ϊ��ʱ�뷽��
/// </summary>
/// <param name="s">��ʼ��</param>
/// <param name="e">��ֹ��</param>
/// <param name="p">��Ҫ�жϵĵ�</param>
/// <returns>���Ϊ0���ڲ�Ϊ1</returns>
int RTInside(Point2D s, Point2D e, Point2D p)//����жϣ����Ϊ0���ڲ�Ϊ1
{
	int cross = (e.x - s.x) * (p.y - s.y) - (p.x - s.x) * (e.y - s.y);
	int flag = (cross > 0) ? 1 : 0;
	return flag;
}
/// <summary>
/// ͨ���ĸ���õ���������ֱ�ߵĽ���
/// </summary>
/// <param name="p1">��һ��ֱ�ߵ�һ����</param>
/// <param name="p2">��һ��ֱ�ߵ���һ����</param>
/// <param name="s1">�ڶ���ֱ�ߵ�һ����</param>
/// <param name="s2">�ڶ���ֱ�ߵ���һ����</param>
/// <returns>����</returns>
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
/// SH�ü��㷨
/// </summary>
/// <param name="polygon">���ü�����ζ���vector</param>
/// <param name="clipbound">�ü��߿򶥵�vector</param>
/// <returns>�ü���ϵĶ����vector</returns>
vector<Point2D> subSH_Clip(vector<Point2D> pg_vec, vector<Point2D> clip_vec)
{
	if (pg_vec.empty())
		return pg_vec;
	int i = 0, flag = 0;
	Point2D start, stop;//�����ö���εı����������յ�
	Point2D sp, ep;//���ô��ڱ߽������������յ�
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
			if (RTInside(sp, ep, stop))//��ǰ��i�������Ƿ��ڱ߽��ڲ�
			{
				if (flag == 0)//ǰһ�����Ƿ������
				{
					flag = 1;//���⵽�ڵ����					
					temp.push_back(GetCrossPoint(sp, ep, start, stop));//���뽻��
				}
				temp.push_back(stop);
			}
			else
			{
				if (flag)
				{
					flag = 0;//���ڵ���
					temp.push_back(GetCrossPoint(sp, ep, start, stop));//���뽻��
				}
			}
			start = stop;
		}
		sp = ep;
		if (temp.size() == 0) return temp;
		pg_vec = temp;//������һ���߲ü���vector�����ü�
		temp.clear();//����м����������һ�βü�
	}
	pg_vec.push_back(pg_vec[0]);
	return pg_vec;
}

