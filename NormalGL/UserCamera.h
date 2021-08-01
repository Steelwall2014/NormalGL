#pragma once
#include "GeoDefine.h"
//#include "Matrix.h"

#include <vector>
using namespace std;
class Matrix;
class Painter;
class Model;
class UserCamera
{
public:
	double eye_x, eye_y, eye_z;
	double look_x, look_y, look_z;
	double up_x, up_y, up_z;
	double z_near, z_far;
	double height, width;
	double aspect, fovy;
	Matrix u, v, n;
	Matrix viewMatrix;
	Matrix perspectMatrix;
public:
	UserCamera() {};
	UserCamera(double eye_x, double eye_y, double eye_z, 
			   double look_x, double look_y, double look_z,
		double up_x, double up_y, double up_z);

	void SetPerspective(double aspect, double fovy, double z_near, double z_far);
	void SetFrustum(double width , double height, double z_near, double z_far);
	void SetEyePos(double eye_x, double eye_y, double eye_z);
	void MoveFrontBack(double distance);
	//void MoveBackward(double distance);

	// Look函数的作用是一样的，给定一个modek或者face，输出绘制的图形
	Painter *Look(Face & triangle);
	Painter *Look(Face &triangle, Matrix &transform);
	vector<Painter *> Look(Model &model, Matrix &transform);
	Painter *Look(Line3D &line, Matrix &rotation);
	// 屏幕上的坐标转换到世界坐标，这个函数是一个简化版，不适用于任意位置的屏幕
	Point3D WindowToWorld(double x, double y);
};