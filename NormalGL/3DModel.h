#pragma once
#include <string>
#include "Tin.h"
//#include "GeoDefine.h"
using namespace std;

//class Tin;

class Model
{
public:
	Tin *vertice;
	vector<Face> faces;
	vector<Point3D> vertex_norm_vecs;
	vector<Point3D> faces_norm_vecs;
	vector<double> vertex_IR;
	vector<double> vertex_IG;
	vector<double> vertex_IB;
	Model() {}
	Model(string filepath);
private:
	void _add_vertex_norm(int index, Matrix &norm);
	void _add_vertex_norm(int index, Point3D norm);
};