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
	//int trishape = 0;//ƽ��Ϊ1��ƽ��Ϊ2��������Ϊ0
	//int spvertex = 0;//������Ϊƽ�׻���ƽ��ʱ���洢���ǵױߵĶ������������䲻Ϊʱ���洢�м䶥�㣨���ڷָ������Σ��Ķ�������
	//int flatline[2] = { 0 };//ƽ����ƽ�״洢ƽ�ߵ���������������ƴ�ӵ�Ϊ������������
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
	// ���������¼ӽ�ȥ�ĵ��id����˵��ʵ���Է���push�����
	void SetLastElemId();
};




#endif
