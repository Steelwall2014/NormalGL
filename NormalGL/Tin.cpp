#include "Graphic.h"
#include "GeoDefine.h"
#include "Tin.h"
Point3D Face::operator[](int index)
{
	return (*this->Tin_grid)[point_index[index]];
}
Face::Face(Tin *Tin_grid, int pt1_index, int pt2_index, int pt3_index)
{
	this->point_index[0] = pt1_index;
	this->point_index[1] = pt2_index;
	this->point_index[2] = pt3_index;
	this->Tin_grid = Tin_grid;
	//DetermineTriShape();
}

Face::Face(Tin *Tin_grid, int pts_index[3])
{
	this->point_index[0] = pts_index[0];
	this->point_index[1] = pts_index[1];
	this->point_index[2] = pts_index[2];
	this->Tin_grid = Tin_grid;
	//DetermineTriShape();
}

Face::Face(const Face &other)
{
	this->point_index[0] = other.point_index[0];
	this->point_index[1] = other.point_index[1];
	this->point_index[2] = other.point_index[2];
	this->Tin_grid = other.Tin_grid;
	//DetermineTriShape();
}

Tin::Tin(Point3D* pts, int point_num)
{
	//Point3D apt;
	points = new vector<Point3D>();
	for (int i = 0; i < point_num; ++i)
	{
		//apt = pts[i];
		points->push_back(pts[i]);
	}
}

Tin::Tin(vector<Point3D> &pts)
{
	points = new vector<Point3D>(pts.size());
	for (int i = 0; i < pts.size(); i++)
		(*points)[i] = pts[i];
}

Point3D& Tin::operator[](int index)
{
	if (index >= points->size())
	{
		throw "Out of range";
	}
	return (*points)[index];
}

bool Tin::push(Point3D pt)
{
	//Point3D newpt = pt;
	points->push_back(pt);
	return true;
}

Point3D Tin::pop(Point3D pt)
{
	Point3D lastpt = (*points)[points->size()-1];
	points->pop_back();
	return lastpt;
}

int Tin::size()
{
	return points->size();
}

void Tin::SetLastElemId()
{
	if (points->size() != 0)
		(*points)[points->size() - 1].id = points->size() - 1;
}
