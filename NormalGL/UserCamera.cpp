#include <vector>
#include <cmath>
#include "Matrix.h"
#include "Painters.h"
#include "UserCamera.h"
#include "3DModel.h"
#include "Status.h"

extern double Kd, Ks, Ka;
extern int shininess;
UserCamera::UserCamera(double eye_x, double eye_y, double eye_z, double look_x, double look_y, double look_z, double up_x, double up_y, double up_z)
{
	// n是看向后的向量，v是指向上方的，u是两个向量叉乘
	this->eye_x = eye_x;
	this->eye_y = eye_y;
	this->eye_z = eye_z;
	this->look_x = look_x;
	this->look_y = look_y;
	this->look_z = look_z;
	this->up_x = up_x;
	this->up_y = up_y;
	this->up_z = up_z;
	n = Matrix(1, 3);
	n.data[0][0] = eye_x - look_x; n.data[0][1] = eye_y - look_y; n.data[0][2] = eye_z - look_z;
	n.Normalize();
	u = Matrix(1, 3);
	u.data[0][0] = up_x; u.data[0][1] = up_y; u.data[0][2] = up_z;
	u = u.CrossProduct(n);
	u.Normalize();
	v = n.CrossProduct(u);
}

void UserCamera::SetPerspective(double aspect, double fovy, double z_near, double z_far)
{
	this->z_near = z_near;
	this->z_far = z_far;
	this->aspect = aspect;
	this->fovy = fovy;
	double dz = z_near - this->eye_z;
	this->height = 2 * dz * tan(fovy / 2);
	this->width = this->height * aspect;
	viewMatrix = Matrix::CreateViewMatrix(*this);
	perspectMatrix = Matrix::CreatePerspectiveMatrix(*this);
}

void UserCamera::SetFrustum(double width, double height, double z_near, double z_far)
{
	this->z_near = z_near;
	this->z_far = z_far;
	this->height = height;
	this->width = width;
	this->aspect = (double)width / (double)height;
	this->fovy = atan(height / 2.0 / z_near) * 2;
	viewMatrix = Matrix::CreateViewMatrix(*this);
	perspectMatrix = Matrix::CreatePerspectiveMatrix(*this);
}

void UserCamera::SetEyePos(double eye_x, double eye_y, double eye_z)
{
	this->eye_x = eye_x;
	this->eye_y = eye_y;
	this->eye_z = eye_z;
	
	n.data[0][0] = eye_x - look_x; n.data[0][1] = eye_y - look_y; n.data[0][2] = eye_z - look_z;
	n.Normalize();
	u.data[0][0] = up_x; u.data[0][1] = up_y; u.data[0][2] = up_z;
	u = u.CrossProduct(n);
	u.Normalize();
	v = n.CrossProduct(u);
}

void UserCamera::MoveFrontBack(double distance)
{
	eye_x += n.data[0][0] * (-distance);
	eye_y += n.data[0][1] * (-distance);
	eye_z += n.data[0][2] * (-distance);
	viewMatrix = Matrix::CreateViewMatrix(*this);
	perspectMatrix = Matrix::CreatePerspectiveMatrix(*this);
}

Painter * UserCamera::Look(Face & triangle)
{
	Point3D projected_points[3];
	for (int i = 0; i < 3; i++)
	{
		Point3D temp = triangle[i] * viewMatrix;
		temp = temp * perspectMatrix;
		if (temp.depth < 0)
			return nullptr;
		temp.x /= temp.depth; temp.y /= temp.depth; temp.z /= temp.depth;
		projected_points[i] = Point3D{ (temp.x - (-1)) * getWindowWidth() / 2.0, (temp.y - (-1)) * getWindowHeight() / 2.0, temp.z };
	}
	TrianglePainter *painter = new TrianglePainter(projected_points[0], projected_points[1], projected_points[2]);
	//projected_points[3] = projected_points[0];
	//TrianglePainter *painter = new TrianglePainter(projected_points);
	//PolygonPainter *painter = new PolygonPainter(projected_points, 4);
	return painter;
}


Painter *UserCamera::Look(Face & triangle, Matrix &transform)
{
	Point3D projected_points[3];
	for (int i = 0; i < 3; i++)
	{
		Point3D temp = triangle[i] * transform;
		temp = temp * viewMatrix;
		temp = temp * perspectMatrix;
		if (temp.depth < 0)
			return nullptr;
		temp.x /= temp.depth; temp.y /= temp.depth; temp.z /= temp.depth;
		projected_points[i] = Point3D{ (temp.x - (-1)) * getWindowWidth() / 2.0, (temp.y - (-1)) * getWindowHeight() / 2.0, temp.z };
	}
	TrianglePainter *painter = new TrianglePainter(projected_points[0], projected_points[1], projected_points[2]);
	//projected_points[3] = projected_points[0];
	//TrianglePainter *painter = new TrianglePainter(projected_points);
	//PolygonPainter *painter = new PolygonPainter(projected_points, 4);
	return painter;
}

vector<Painter *> UserCamera::Look(Model &model, Matrix &transform)
{
	vector<Painter *> painters;
	/**********************************真实感着色参数设置***********************************/
	// 各色的入射光强
	double Ir = 0, Ig = 1.5, Ib = 2.5;
	// 给定入射环境光光强
	double Ia = 0.5;
	// 给定光源位于相机处
	Matrix light_source_vec = Matrix(1, 3);
	light_source_vec.data[0][0] = -n.data[0][0];
	light_source_vec.data[0][1] = -n.data[0][1];
	light_source_vec.data[0][2] = -n.data[0][2];
	// 给定漫反射常数、镜面反射常数、环境光漫反射常数
	//double Kd = 0.4, Ks = 1, Ka = 0.2;
	// 给定高光指数
	//int shininess = 20;
	/**********************************真实感着色参数设置***********************************/

	for (int i = 0; i < model.vertice->size(); i++)
	{
		Point3D after_trans_norm_point3d = model.vertex_norm_vecs[i] * transform;//inv_transpose_transform;
		after_trans_norm_point3d.Normalize();
		//if (after_trans_norm_point3d.y < -0.98 && (*model.vertice)[i].y < 0)// && after_trans_norm_point3d.z < 0.01 && after_trans_norm_point3d.z > -0.01)
		//if ((*model.vertice)[i].y < -250 && (*model.vertice)[i].x <10 && (*model.vertice)[i].x > -10)
		//	double aaa = 0;
		Matrix after_trans_norm = Matrix(1, 3);
		after_trans_norm.data[0][0] = after_trans_norm_point3d.x;
		after_trans_norm.data[0][1] = after_trans_norm_point3d.y;
		after_trans_norm.data[0][2] = after_trans_norm_point3d.z;
		Matrix V = Matrix(1, 3);						// V是顶点到观察点的向量
		Point3D after_trans_point = (*model.vertice)[i] * transform;
		V.data[0][0] = eye_x - after_trans_point.x;
		V.data[0][1] = eye_y - after_trans_point.y;
		V.data[0][2] = eye_z - after_trans_point.z;
		V.Normalize();
		// R是镜面反射后的向量
		Matrix R = light_source_vec - after_trans_norm * (2.0 * (light_source_vec.DotProduct(after_trans_norm)));
		R.Normalize();
		double LN = after_trans_norm.DotProduct(-light_source_vec);
		double RV = R.DotProduct(V);
		// 漫反射光强
		double Idr = Ir * Kd * max(0.0, LN);
		double Idg = Ig * Kd * max(0.0, LN);
		double Idb = Ib * Kd * max(0.0, LN);
		// 镜面反射
		double Isr = Ir * Ks * pow(max(0, RV), shininess);
		double Isg = Ig * Ks * pow(max(0, RV), shininess);
		double Isb = Ib * Ks * pow(max(0, RV), shininess);
		// 环境光
		double Ie = Ia * Ka;

		double IR = Ie + Idr + Isr;
		double IG = Ie + Idg + Isg;
		double IB = Ie + Idb + Isb;
		model.vertex_IR[i] = IR;
		model.vertex_IG[i] = IG;
		model.vertex_IB[i] = IB;
	}
	for (int i = 0; i < model.faces.size(); i++)
	{
		Point3D projected_points[3];

		// 后向面消隐
		if (status.cull_face)
		{
			Matrix assumed_view_direction = Matrix(1, 3);
			assumed_view_direction.data[0][0] = 0;
			assumed_view_direction.data[0][1] = 1;
			assumed_view_direction.data[0][2] = 0;
			Point3D face_norm = model.faces_norm_vecs[model.faces[i].face_id];
			if (face_norm.DotProduct(assumed_view_direction) > 0)
				continue;
		}	 

		for (int j = 0; j < 3; j++)
		{
			Point3D temp = model.faces[i][j] * transform * viewMatrix * perspectMatrix;
			if (temp.depth < 0)
				break;
			temp.id = model.faces[i][j].id;
			temp.x /= temp.depth; temp.y /= temp.depth; temp.z /= temp.depth;
			temp.x = (temp.x - (-1)) * getWindowWidth() / 2.0;
			temp.y = (temp.y - (-1)) * getWindowHeight() / 2.0;
			projected_points[j] = temp;
		}
		TrianglePainter *painter = new TrianglePainter(projected_points[0], projected_points[1], projected_points[2]);
		painter->model = &model;
		painters.push_back(painter);
	}
	return painters;
}

Painter *UserCamera::Look(Line3D &line, Matrix &transform)
{
	Point2D projected_points[2];
	for (int i = 0; i < 2; i++)
	{
		Point3D temp = line.pts[i] * transform;
		temp = temp * viewMatrix;
		temp = temp * perspectMatrix;
		if (temp.depth < 0)
			return nullptr;
		temp.x /= temp.depth; temp.y /= temp.depth; temp.z /= temp.depth;
		projected_points[i] = Point2D{ (temp.x - (-1)) * getWindowWidth() / 2.0, (temp.y - (-1)) * getWindowHeight() / 2.0 };
	}
	LinePainter *painter = new LinePainter(projected_points[0].x, projected_points[0].y, projected_points[1].x, projected_points[1].y);
	//projected_points[3] = projected_points[0];
	//TrianglePainter *painter = new TrianglePainter(projected_points);
	//PolygonPainter *painter = new PolygonPainter(projected_points, 4);
	return painter;
}

Point3D UserCamera::WindowToWorld(double x, double y)
{
	// 这个仅仅适用于UserCamera参数为 0, -500, 0, 0, 0, 0, 0, 0, 1 的情况
	Point3D world_point;
	world_point.x = eye_x + (-n.data[0][0] * z_near);
	world_point.y = eye_y + (-n.data[0][1] * z_near);
	world_point.z = eye_z + (-n.data[0][2] * z_near);
	world_point.x -= width / 2 - x;
	world_point.z -= height / 2 - y;
	return world_point;
}
