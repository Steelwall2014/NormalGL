#pragma once
#ifndef _tin_H
#define _tin_H

#include "Graphic.h"
//#include "GeoDefine.h"
class Point3D; 
class Face
{
public:
	int point_count = 3;
	int face_id;
	//int point_count = 3;
	//int trishape = 0;//平底为1，平顶为2，都不是为0
	//int spvertex = 0;//三角形为平底或者平顶时，存储不是底边的顶点索引；当其不为时，存储中间顶点（用于分割三角形）的顶点索引
	//int flatline[2] = { 0 };//平顶或平底存储平边的两个顶点索引，拼接的为上下两个顶点
	int point_index[3];
	Tin *Tin_grid;
	Point3D operator[](int index);
	Face() {}
	Face(Tin *Tin_grid, int pt1_index, int pt2_index, int pt3_index);
	Face(Tin *Tin_grid, int pts_index[3]);
	Face(const Face &other);
	//void DetermineTriShape();
};
class Tin
{
private:
	vector<Point3D> *points;
public:
	Tin() 
	{
		points = new vector<Point3D>();
	}
	~Tin()
	{
		delete points;
	}
	Tin(Point3D*, int);
	Tin(vector<Point3D> &);
	Point3D& operator[](int);
	bool push(Point3D);
	Point3D pop(Point3D);
	int size();
	// 用来设置新加进去的点的id（话说其实可以放在push方法里）
	void SetLastElemId();
};




#endif
